#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
optimize_NNMLP.py

Optuna search for the MLP-NN res_hs model

Fixed (from the first run):  activation='tanh', batch_size=128,
                             n_iter_no_change=20, validation_fraction=0.15
Searched:                    hidden_layer_sizes, norm_range, alpha,
                             learning_rate_init (all in narrower ranges)

Differences vs. the first serial script
---------------------------------------
1. RMSE is computed in physical units of res_hs (after inverse-transform),
   and scalers are fit on each fold's TRAINING data only. In the first run
   the RMSE was measured in normalized space, where a [-1, 1] target has
   twice the range of a [0, 1] target, so '-1_1' trials had ~2x inflated
   scores and were unfairly penalized. norm_range is therefore re-opened
   here.
2. The study is stored in SQLite (optuna_mlp_narrow.db): if the run is
   interrupted, re-running the script resumes where it stopped, and nothing
   is lost. best_mlp_params.json is rewritten every time a new best trial
   appears.
3. The best configuration found in the first run is enqueued as the first
   two trials (once with each norm_range), giving a direct fair comparison
   of '0_1' vs '-1_1' on the same configuration.
4. Optional parallelism INSIDE each trial (N_JOBS): the 4 seeds x 3 folds
   of one trial are run concurrently. N_JOBS = 1 gives the plain serial
   behavior. Try N_JOBS = 12 for a speed test.

No execution here -- run on the server.
"""

# N_JOBS has to be known before numpy is imported (BLAS thread pinning)
N_JOBS = 12          # 1 = serial (BLAS uses all cores); >1 = parallel fits, 1 BLAS thread each

import os
if N_JOBS > 1:
    for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
               "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        os.environ[_v] = "1"

import ast
import json
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

CV_SPLITS = 3
OPT_SEEDS = (1, 2, 3, 5)     # seeds averaged per trial (was 6; 4 keeps the noise reasonable)
MAX_ITER_CAP = 3000          # early_stopping cuts this short in practice
N_TRIALS = 40                # total trials in the study (resumable)

# fixed from the first run
FIXED = dict(activation='tanh', batch_size=128,
             n_iter_no_change=20, validation_fraction=0.15)

ARCH_OPTIONS = [ (2,),(5,),(10,),(20,),(50,),(80,),(150,),
    (5,5),(10,10),(20,20),(50,50),(80,80),(150,150)]

ALPHA_RANGE = (1e-3, 2e-2)
LR_RANGE = (1e-4, 5e-3)

# best configuration of the first run (trial 62) -> enqueued as the first trials
BEST_FIRST_RUN = {'hidden_layer_sizes': '(20,)', 'alpha': 0.00423759670849982,
                  'learning_rate_init': 0.0005815387120251358}

STUDY_NAME = "mlp_res_hs_narrow"
STORAGE = "sqlite:///optuna_mlp_narrow.db"
OUT_PARAMS_JSON = "best_mlp_params.json"
OUT_TRIALS_CSV = "optuna_trials_narrow.csv"
OUT_HISTORY_PNG = "optuna_history_narrow.png"
OUT_IMPORTANCE_PNG = "optuna_importance_narrow.png"

RANGE_MAP = {'-1_1': (-1, 1), '0_1': (0, 1)}


# ============================================================
# Load data once
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
FOLD_INDICES = list(gkf.split(X_raw, y_raw, groups=groups))   # same folds for every trial


# ============================================================
# One MLP fit -> validation RMSE in physical units (m)
# ============================================================

def fit_one(X, y, tr_idx, te_idx, p, seed):
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
# Optuna objective
# ============================================================

def objective(trial):
    arch_str = trial.suggest_categorical('hidden_layer_sizes', [str(a) for a in ARCH_OPTIONS])
    p = dict(FIXED)
    p['arch'] = ast.literal_eval(arch_str)
    p['norm_range'] = trial.suggest_categorical('norm_range', ['0_1', '-1_1'])
    p['alpha'] = trial.suggest_float('alpha', *ALPHA_RANGE, log=True)
    p['learning_rate_init'] = trial.suggest_float('learning_rate_init', *LR_RANGE, log=True)

    results = Parallel(n_jobs=N_JOBS)(
        delayed(fit_one)(X_raw, y_raw, tr, te, p, seed)
        for (tr, te) in FOLD_INDICES for seed in OPT_SEEDS)

    vals = np.array(results, dtype='float64')
    if not np.isfinite(vals).all():
        return float('nan')
    trial.set_user_attr('rmse_std_over_fits', float(vals.std()))   # noise level of this score
    return float(vals.mean())


def save_best(study, trial):
    """Rewrite best_mlp_params.json whenever a new best trial appears."""
    if study.best_trial.number != trial.number:
        return
    bp = dict(study.best_params)
    bp['hidden_layer_sizes'] = list(ast.literal_eval(bp['hidden_layer_sizes']))
    bp.update(FIXED)
    bp['cv_rmse_res_hs_m'] = study.best_value
    with open(OUT_PARAMS_JSON, 'w') as f:
        json.dump(bp, f, indent=2)


# ============================================================
# Run
# ============================================================

if __name__ == "__main__":

    study = optuna.create_study(
        study_name=STUDY_NAME, storage=STORAGE, load_if_exists=True,
        direction='minimize', sampler=TPESampler(seed=0))

    if len(study.trials) == 0:
        # first-run best config, once per normalization -> fair '0_1' vs '-1_1' comparison
        for nr in ('0_1', '-1_1'):
            study.enqueue_trial({**BEST_FIRST_RUN, 'norm_range': nr})

    n_finished = len([t for t in study.trials if t.state.is_finished()])
    remaining = max(N_TRIALS - n_finished, 0)
    print(f"{n_finished} trials already in the study; running {remaining} more "
          f"(N_JOBS={N_JOBS}, {len(FOLD_INDICES) * len(OPT_SEEDS)} fits per trial)")

    study.optimize(objective, n_trials=remaining, callbacks=[save_best])

    print("\nBest trial:")
    print(f"  RMSE (res_hs units, m) = {study.best_value:.5f}")
    print(f"  params = {study.best_params}")

    study.trials_dataframe().to_csv(OUT_TRIALS_CSV, sep='\t', index=False, float_format='%.6f')
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
        print(f"Plotting skipped ({e})")
