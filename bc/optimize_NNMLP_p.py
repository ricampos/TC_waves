#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
optimize_NNMLP_p.py  (parallel version)

Hyperparameter search for the MLP-NN res_hs model, using Optuna (TPE
sampler) with storm-grouped cross-validation (GroupKFold on cid).

Parallelization
---------------
Every MLP fit (one trial x one CV fold x one seed) is independent, so they
are farmed out to N_JOBS worker processes with joblib. To keep all cores
busy, trials are evaluated in small batches via Optuna's ask/tell
interface: TRIALS_PER_BATCH trials are asked at once, all their
(fold, seed) fits are run together in one parallel pool, then the scores
are told back to the study. With 3 folds x 6 seeds = 18 fits per trial and
TRIALS_PER_BATCH = 3, that is 54 tasks per batch, ~2 waves on 28 cores.
TPE uses constant_liar=True so the trials asked within one batch spread out
instead of all landing on the same point.

BLAS/OpenMP threads are pinned to 1 per worker (env vars below) so that 28
workers do not each start 28 threads and thrash the CPU.

Scoring (changed vs. previous version)
--------------------------------------
- Scalers are fit on the TRAINING part of each fold only (no min/max leak
  from the validation storms).
- RMSE is computed in physical units of res_hs (after inverse-transforming
  the prediction), not in normalized space. In normalized space a [-1, 1]
  target has twice the range of a [0, 1] target, so its RMSE is inflated by
  ~2x and the search would be biased toward '0_1' regardless of real skill.

With solver='adam', scikit-learn ignores learning_rate (schedule), power_t,
and momentum/nesterovs_momentum, so they are not searched.

Output: best_mlp_params.json, optuna_trials.csv, optuna_history.png,
optuna_importance.png.
"""

# --- pin BLAS/OpenMP threads BEFORE numpy/sklearn are imported -------------
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_v] = "1"

import ast
import json
import time
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import netCDF4 as nc
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import GroupKFold
from sklearn.neural_network import MLPRegressor
from joblib import Parallel, delayed

import optuna
from optuna.samplers import TPESampler
optuna.logging.set_verbosity(optuna.logging.WARNING)


# ============================================================
# Data reading (same as featureSel_dataScreening.py)
# ============================================================

def read_tc_dataset(fname, OMVARS, LVARS):
    ds = nc.Dataset(fname, 'r')
    time_epoch = ds.variables['time'][:].astype('float64')
    time_dt = pd.to_datetime(time_epoch, unit='s')
    tlag = ds.variables['tlag'][:].astype('int32')

    raw = {'time_epoch': time_epoch, 'time_dt': time_dt, 'tlag': tlag}
    for v in OMVARS:
        arr = ds.variables[v][:]
        if np.ma.isMaskedArray(arr):
            arr = arr.filled(np.nan) if arr.dtype.kind == 'f' else arr.data
        raw[v] = arr
    for v in LVARS:
        arr = ds.variables[v][:, :]
        if np.ma.isMaskedArray(arr):
            arr = arr.filled(np.nan)
        raw[v] = arr

    df = pd.DataFrame({'time': time_dt})
    for v in OMVARS:
        vals = raw[v]
        if vals.dtype.kind in ('S', 'U', 'O'):
            vals = np.array([x.decode() if isinstance(x, bytes) else str(x).strip()
                              for x in vals])
        df[v] = vals

    for v in LVARS:
        for i, lag in enumerate(tlag):
            df[f'{v}_lag{int(lag)}'] = raw[v][i, :]

    ds.close()
    return raw, df


def read_feature_list(fname):
    with open(fname) as f:
        lines = [ln.strip() for ln in f if ln.strip()]
    if not lines:
        raise ValueError(f"{fname} is empty")
    first_cols = lines[0].split('\t')
    if first_cols[0].lower() == 'feature':
        features = [ln.split('\t')[0] for ln in lines[1:]]
    else:
        features = [ln.split('\t')[0] for ln in lines]
    return features


# ============================================================
# Config
# ============================================================

FNAME = "TC_ModelObs_CyclInfo.nc"
FEATURE_FILE = "FS_MLPNN_res.txt"
TARGET = "res_hs"

CV_SPLITS = 3                     # fewer folds during search (speed)
OPT_SEEDS = (1, 2, 3, 5, 10, 15)  # seeds averaged per trial (robustness, not cherry-picking)
MAX_ITER_CAP = 10000              # early_stopping cuts this short in practice
N_TRIALS = 150                    # total Optuna trials

N_JOBS = 28                       # worker processes
TRIALS_PER_BATCH = 3              # trials evaluated concurrently (3 x 18 fits = 54 tasks)

ARCH_OPTIONS = [(2,), (5,), (10,), (20,), (50,), (80,), (150,),
                (5, 5), (10, 10), (20, 20), (50, 50), (80, 80), (150, 150),
                (100, 100, 50), (150, 100, 50)]

OUT_PARAMS_JSON = "best_mlp_params.json"
OUT_TRIALS_CSV = "optuna_trials.csv"
OUT_HISTORY_PNG = "optuna_history.png"
OUT_IMPORTANCE_PNG = "optuna_importance.png"

RANGE_MAP = {'-1_1': (-1, 1), '0_1': (0, 1)}


# ============================================================
# Load data once (module-level)
# ============================================================

LVARS = [
    'clat', 'clon', 'vmx', 'rmw', 'vfm', 'cdir_sin', 'cdir_cos',
    'r34', 'r34_ne', 'r34_nw', 'r34_se', 'r34_sw',
    'r50', 'r50_ne', 'r50_nw', 'r50_se', 'r50_sw',
    'r64', 'r64_ne', 'r64_nw', 'r64_se', 'r64_sw']

OMVARS = [
    'jlday', 'obs_ID', 'cid', 'mlat', 'mlon', 'csec', 'ccdist', 'ccangle',
    'ohs', 'ownd', 'mhs', 'mtp', 'mwnd', 'hhs', 'hwnd',
    'ghs_c', 'ghs_em', 'gwnd_c', 'gwnd_em']

ds, df = read_tc_dataset(FNAME, OMVARS, LVARS)
df['res_hs'] = df['mhs'] - df['ohs']

selected_features = read_feature_list(FEATURE_FILE)
missing = [f for f in selected_features if f not in df.columns]
if missing:
    raise ValueError(f"Selected features not found in df: {missing}")

keep_cols = selected_features + [TARGET, 'cid']
mdf = df[keep_cols].dropna(subset=selected_features + [TARGET]).reset_index(drop=True)
print(f"{mdf.shape[0]}/{df.shape[0]} rows valid (no NaN in selected features/target)")

X_raw = mdf[selected_features].values.astype('float64')
y_raw = mdf[TARGET].values.reshape(-1, 1).astype('float64')
groups = mdf['cid'].values

gkf = GroupKFold(n_splits=CV_SPLITS)
# split once (same folds reused across all trials -- fair comparison between trials)
FOLD_INDICES = list(gkf.split(X_raw, y_raw, groups=groups))


# ============================================================
# One MLP fit (runs in a worker process)
# ============================================================

def fit_one(X, y, tr_idx, te_idx, p, seed):
    """Train one MLP on one fold with one seed; return validation RMSE in
    physical units of the target (res_hs)."""
    try:
        rng_range = RANGE_MAP[p['norm_range']]
        xs = MinMaxScaler(feature_range=rng_range).fit(X[tr_idx])
        ys = MinMaxScaler(feature_range=rng_range).fit(y[tr_idx])

        model = MLPRegressor(
            hidden_layer_sizes=p['arch'], activation=p['activation'], solver='adam',
            alpha=p['alpha'], batch_size=p['batch_size'],
            learning_rate_init=p['learning_rate_init'],
            max_iter=MAX_ITER_CAP, shuffle=True, random_state=seed,
            early_stopping=True, n_iter_no_change=p['n_iter_no_change'],
            validation_fraction=p['validation_fraction'],
            beta_1=0.9, beta_2=0.999, epsilon=1e-08,
        )
        model.fit(xs.transform(X[tr_idx]), ys.transform(y[tr_idx])[:, 0])

        pred_n = model.predict(xs.transform(X[te_idx]))
        pred = ys.inverse_transform(pred_n.reshape(-1, 1))[:, 0]
        return float(np.sqrt(np.mean((pred - y[te_idx, 0]) ** 2)))
    except Exception:
        return float('nan')


# ============================================================
# Search space
# ============================================================

def suggest_params(trial):
    arch_str = trial.suggest_categorical('hidden_layer_sizes', [str(a) for a in ARCH_OPTIONS])
    bs = trial.suggest_categorical('batch_size', ['32', '64', '128', '256', 'auto'])
    return {
        'arch': ast.literal_eval(arch_str),
        'activation': trial.suggest_categorical('activation', ['tanh', 'relu', 'logistic']),
        'norm_range': trial.suggest_categorical('norm_range', ['-1_1', '0_1']),
        'alpha': trial.suggest_float('alpha', 1e-6, 1e-1, log=True),
        'learning_rate_init': trial.suggest_float('learning_rate_init', 1e-5, 1e-1, log=True),
        'batch_size': bs if bs == 'auto' else int(bs),
        'n_iter_no_change': trial.suggest_int('n_iter_no_change', 10, 30),
        'validation_fraction': trial.suggest_float('validation_fraction', 0.1, 0.2),
    }


# ============================================================
# Run the study
# ============================================================

if __name__ == "__main__":

    sampler = TPESampler(seed=0, constant_liar=True)
    study = optuna.create_study(direction='minimize', sampler=sampler)

    fits_per_trial = len(FOLD_INDICES) * len(OPT_SEEDS)
    print(f"{N_TRIALS} trials x {fits_per_trial} fits/trial = "
          f"{N_TRIALS * fits_per_trial} MLP fits on {N_JOBS} workers "
          f"({TRIALS_PER_BATCH} trials per batch)")

    t0 = time.time()
    n_done = 0
    # one persistent pool for the whole search (workers are reused across batches)
    with Parallel(n_jobs=N_JOBS, backend='loky') as parallel:
        while n_done < N_TRIALS:
            k = min(TRIALS_PER_BATCH, N_TRIALS - n_done)
            trials = [study.ask() for _ in range(k)]
            plist = [suggest_params(t) for t in trials]

            tasks = [(ti, tr, te, seed)
                     for ti in range(k)
                     for (tr, te) in FOLD_INDICES
                     for seed in OPT_SEEDS]
            results = parallel(
                delayed(fit_one)(X_raw, y_raw, tr, te, plist[ti], seed)
                for (ti, tr, te, seed) in tasks)

            # regroup per trial and tell Optuna
            per_trial = [[] for _ in range(k)]
            for (ti, _, _, _), r in zip(tasks, results):
                per_trial[ti].append(r)

            for ti, trial in enumerate(trials):
                vals = np.array(per_trial[ti], dtype='float64')
                if np.isfinite(vals).all():
                    study.tell(trial, float(vals.mean()))
                else:
                    study.tell(trial, state=optuna.trial.TrialState.FAIL)

            n_done += k
            best = study.best_value if any(
                t.state == optuna.trial.TrialState.COMPLETE for t in study.trials) else float('nan')
            print(f"[{n_done:4d}/{N_TRIALS}] best RMSE so far = {best:.5f}  "
                  f"elapsed {(time.time() - t0) / 60:.1f} min", flush=True)

    print("\nBest trial:")
    print(f"  RMSE (res_hs units, m) = {study.best_value:.5f}")
    print(f"  params = {study.best_params}")

    best_params = dict(study.best_params)
    best_params['hidden_layer_sizes'] = list(ast.literal_eval(best_params['hidden_layer_sizes']))

    with open(OUT_PARAMS_JSON, 'w') as f:
        json.dump(best_params, f, indent=2)
    print(f"\nSaved best params -> {OUT_PARAMS_JSON}")

    trials_df = study.trials_dataframe()
    trials_df.to_csv(OUT_TRIALS_CSV, sep='\t', index=False, float_format='%.6f')
    print(f"Saved all trials  -> {OUT_TRIALS_CSV}")

    try:
        from optuna.visualization.matplotlib import (
            plot_optimization_history, plot_param_importances)

        plot_optimization_history(study)
        plt.tight_layout()
        plt.savefig(OUT_HISTORY_PNG, dpi=200, bbox_inches='tight')
        plt.close('all')

        plot_param_importances(study)
        plt.tight_layout()
        plt.savefig(OUT_IMPORTANCE_PNG, dpi=200, bbox_inches='tight')
        plt.close('all')

        print(f"Saved plots       -> {OUT_HISTORY_PNG}, {OUT_IMPORTANCE_PNG}")
    except Exception as e:
        print(f"Plotting skipped ({e}); best_mlp_params.json and "
              f"optuna_trials.csv are still saved.")
