#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
optimize_NNMLP.py

Hyperparameter search for the MLP-NN res_hs model, using Optuna (TPE
sampler) with storm-grouped cross-validation (GroupKFold on cid).

IMPORTANT: with solver='adam' (used throughout), scikit-learn's
MLPRegressor IGNORES `learning_rate` (constant/invscaling/adaptive),
`power_t`, and `momentum`/`nesterovs_momentum` -- those are sgd-only.
They are not included in the search space below because tuning them
would have no effect on the trained model.

What's searched instead (the parameters that actually matter for Adam):
  - hidden_layer_sizes  (architecture)
  - activation          (tanh / relu / logistic)
  - normalization range for X/y ([-1,1] vs [0,1])
  - alpha               (L2 regularization)
  - learning_rate_init
  - batch_size
  - n_iter_no_change, validation_fraction (early stopping controls)

Each trial's score is the MEAN validation RMSE (normalized space) across
GroupKFold folds AND across a small set of seeds -- averaging across
seeds, not taking the best one, so a hyperparameter set isn't rewarded
for getting lucky with one random initialization.

Output: best_mlp_params.json (feed this into the final training script),
optuna_trials.csv (every trial's params + score), and two diagnostic
plots (optimization history, parameter importance).

No execution here -- run this on the server.
"""

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

CV_SPLITS = 3           # fewer folds during search (speed); final training uses more
OPT_SEEDS = (1,2,3,5,10,15)      # seeds averaged per trial (robustness, not cherry-picking)
MAX_ITER_CAP = 10000     # early_stopping cuts this short in practice
N_TRIALS = 150          # Optuna trials -- adjust based on available time/compute

ARCH_OPTIONS = [ (2,),(5,),(10,),(20,),(50,),(80,),(150,),
    (5,5),(10,10),(20,20),(50,50),(80,80),(150,150),
    (100, 100, 50), (150, 100, 50),
]


OUT_PARAMS_JSON = "best_mlp_params.json"
OUT_TRIALS_CSV = "optuna_trials.csv"
OUT_HISTORY_PNG = "optuna_history.png"
OUT_IMPORTANCE_PNG = "optuna_importance.png"


# ============================================================
# Load data once (module-level, reused by every trial)
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

X_raw = mdf[selected_features].values
y_raw = mdf[TARGET].values.reshape(-1, 1)
groups = mdf['cid'].values

RANGE_MAP = {'-1_1': (-1, 1), '0_1': (0, 1)}
gkf = GroupKFold(n_splits=CV_SPLITS)
# split once (same folds reused across all trials -- fair comparison between trials)
FOLD_INDICES = list(gkf.split(X_raw, y_raw, groups=groups))


# ============================================================
# Optuna objective
# ============================================================

def objective(trial):
    hidden_layer_sizes = trial.suggest_categorical(
        'hidden_layer_sizes', [str(a) for a in ARCH_OPTIONS])  # Optuna wants hashable scalars
    arch = eval(hidden_layer_sizes)  # back to tuple, e.g. "(100, 100)" -> (100, 100)

    activation = trial.suggest_categorical('activation', ['tanh', 'relu', 'logistic'])
    norm_range = trial.suggest_categorical('norm_range', ['-1_1', '0_1'])
    alpha = trial.suggest_float('alpha', 1e-6, 1e-1, log=True)
    learning_rate_init = trial.suggest_float('learning_rate_init', 1e-5, 1e-1, log=True)
    batch_size = trial.suggest_categorical('batch_size', ['32', '64', '128', '256', 'auto'])
    batch_size = batch_size if batch_size == 'auto' else int(batch_size)
    n_iter_no_change = trial.suggest_int('n_iter_no_change', 10, 30)
    validation_fraction = trial.suggest_float('validation_fraction', 0.1, 0.2)

    xscaler = MinMaxScaler(feature_range=RANGE_MAP[norm_range])
    yscaler = MinMaxScaler(feature_range=RANGE_MAP[norm_range])
    Xn = xscaler.fit_transform(X_raw)
    yn = yscaler.fit_transform(y_raw)[:, 0]

    fold_scores = []
    for tr_idx, te_idx in FOLD_INDICES:
        seed_scores = []
        for seed in OPT_SEEDS:
            model = MLPRegressor(
                hidden_layer_sizes=arch, activation=activation, solver='adam',
                alpha=alpha, batch_size=batch_size,
                learning_rate_init=learning_rate_init,
                max_iter=MAX_ITER_CAP, shuffle=True, random_state=seed,
                early_stopping=True, n_iter_no_change=n_iter_no_change,
                validation_fraction=validation_fraction,
                beta_1=0.9, beta_2=0.999, epsilon=1e-08,
            )
            model.fit(Xn[tr_idx], yn[tr_idx])
            pred = model.predict(Xn[te_idx])
            rmse = np.sqrt(np.mean((pred - yn[te_idx]) ** 2))
            seed_scores.append(rmse)
        fold_scores.append(np.mean(seed_scores))

    return float(np.mean(fold_scores))


# ============================================================
# Run the study
# ============================================================

if __name__ == "__main__":

    study = optuna.create_study(direction='minimize', sampler=TPESampler(seed=0))
    study.optimize(objective, n_trials=N_TRIALS, show_progress_bar=True)

    print("\nBest trial:")
    print(f"  normalized RMSE = {study.best_value:.5f}")
    print(f"  params = {study.best_params}")

    # reconstruct architecture tuple for the saved JSON (was stored as a
    # string key for Optuna's categorical sampler)
    best_params = dict(study.best_params)
    best_params['hidden_layer_sizes'] = list(eval(best_params['hidden_layer_sizes']))

    with open(OUT_PARAMS_JSON, 'w') as f:
        json.dump(best_params, f, indent=2)
    print(f"\nSaved best params -> {OUT_PARAMS_JSON}")

    # all trials, for your own inspection / sanity-checking the search
    trials_df = study.trials_dataframe()
    trials_df.to_csv(OUT_TRIALS_CSV, sep='\t', index=False, float_format='%.6f')
    print(f"Saved all trials  -> {OUT_TRIALS_CSV}")

    # diagnostic plots
    try:
        from optuna.visualization.matplotlib import (
            plot_optimization_history, plot_param_importances)

        fig = plot_optimization_history(study)
        plt.tight_layout()
        plt.savefig(OUT_HISTORY_PNG, dpi=200, bbox_inches='tight')
        plt.close('all')

        fig = plot_param_importances(study)
        plt.tight_layout()
        plt.savefig(OUT_IMPORTANCE_PNG, dpi=200, bbox_inches='tight')
        plt.close('all')

        print(f"Saved plots       -> {OUT_HISTORY_PNG}, {OUT_IMPORTANCE_PNG}")
    except Exception as e:
        print(f"Plotting skipped ({e}); best_mlp_params.json and "
              f"optuna_trials.csv are still saved.")
