#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
optimize_NNMLP_gs.py

MLP-NN hyperparameter search with GridSearchCV inside a leave-one-storm-out
loop (outer: each storm is held out once; inner: storm-grouped 3-fold CV on
the remaining storms).

Usage:  python optimize_NNMLP_gs.py <scrmt> <rst>
    scrmt : 0 = RMSE, 1 = MAE, 2 = R2  -> metric used to pick the best model
            (all three metrics are always recorded)
    rst   : random seed of the MLP
Optional environment variable N_JOBS (default -1 = all cores) for GridSearchCV.

"""

import os
os.environ['PYTHONWARNINGS'] = 'ignore'   # also reaches the joblib workers
                                          # (warnings.filterwarnings only affects this process)
import sys
import warnings; warnings.filterwarnings("ignore")
import copy
import numpy as np
import pandas as pd
import netCDF4 as nc
from sklearn.base import BaseEstimator, RegressorMixin
from sklearn.neural_network import MLPRegressor
from sklearn.model_selection import GridSearchCV, GroupKFold


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
        raw[v] = arr  # shape (n_tlag, n_time)

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


def normalization(*args):
    ''' Normalization Run
    input X[time, variables], method, min(x), max(x)
    '''

    if len(args) == 3:
        method = 1
        x=np.copy(args[0]).astype('f')
        xp1=np.copy(args[1]).astype('f')
        xp2=np.copy(args[2]).astype('f')
    elif len(args) == 4:
        x=np.copy(args[0]).astype('f')
        method=int(np.copy(args[1]))
        xp1=np.copy(args[2]).astype('f')
        xp2=np.copy(args[3]).astype('f')

    elif len(args) < 3:
        sys.exit(' ERROR! More inputs are requested')

    x=np.atleast_2d(x)

    # normalization to create arrays between 0 and 1.
    if method == 1:
        x1 = ( x - xp1 ) / ( xp2 - xp1 )

    # same as 1 with log function applied
    elif method ==2:
        x = np.copy(np.log(x-xp1+1.))
        x1 = x / xp2

    # normalization to create arrays between -1 and 1.
    if method == 3:
        x1 = 2.*(( x - xp1 ) / ( xp2 - xp1 ))-1.

    # same as 3 with log function applied
    if method == 4:
        x = np.copy(np.log(x-xp1+1.))
        x1 = 2.*(x / xp2)-1.

    # subtract mean and divide by 3*standardDeviation
    elif method == 5:
        x1= (x-xp1) / (3.*xp2)

    # non-negative values
    elif method ==6:
        x1 = x - xp1

    return x1,xp1,xp2


def denormalization(*args):
    ''' DeNormalization
    input X[time, variables], method, min(x), max(x)
    '''
    x=np.copy(args[0]).astype('f')
    method=np.copy(args[1]).astype('i')
    xp1=np.copy(args[2]).astype('f')
    xp2=np.copy(args[3]).astype('f')
    if method==1:
        y=(x*(xp2-xp1)+xp1)
    elif method==2:
        y=x*xp2
        y=np.exp(y)-1.+xp1
    elif method==3:
        y=(((x+1.)/2.)*(xp2-xp1)+xp1)
    elif method==4:
        y=((x+1.)/2.)*xp2
        y=np.exp(y)-1.+xp1
    elif method==5:
        y=(x*(3*xp2)+xp1)
    elif method==6:
        y=(x+xp1)
    return y


class StormESMLP(BaseEstimator, RegressorMixin):
    """MLPRegressor whose early stopping uses whole held-out STORMS.

    X must carry the storm code in its LAST column (ignored by the network).
    Inside fit(): val_storm_frac of the storms are held out (random, fixed by
    random_state), the net is trained epoch by epoch (partial_fit) on the rest,
    and training stops when the RMSE on the held-out storms has not improved
    for `patience` epochs; the best weights are restored.
    """
    def __init__(self, hidden_layer_sizes=(20,), alpha=1e-3, learning_rate_init=1e-3,
                 batch_size=128, max_epochs=2000, patience=20,
                 val_storm_frac=0.15, random_state=0):
        self.hidden_layer_sizes = hidden_layer_sizes
        self.alpha = alpha
        self.learning_rate_init = learning_rate_init
        self.batch_size = batch_size
        self.max_epochs = max_epochs
        self.patience = patience
        self.val_storm_frac = val_storm_frac
        self.random_state = random_state

    def fit(self, X, y):
        g, Xf = X[:, -1], X[:, :-1]
        ug = np.unique(g)
        rng = np.random.RandomState(self.random_state)
        nval = max(1, int(round(self.val_storm_frac * len(ug))))
        isv = np.isin(g, rng.choice(ug, nval, replace=False))
        Xtr, ytr, Xva, yva = Xf[~isv], y[~isv], Xf[isv], y[isv]

        net = MLPRegressor(hidden_layer_sizes=self.hidden_layer_sizes, activation='tanh',
                           solver='adam', alpha=self.alpha,
                           learning_rate_init=self.learning_rate_init,
                           batch_size=self.batch_size, shuffle=True,
                           random_state=self.random_state)
        best, wait, best_w = np.inf, 0, None
        for ep in range(self.max_epochs):
            net.partial_fit(Xtr, ytr)                      # one epoch
            v = np.mean((net.predict(Xva) - yva) ** 2)
            if v < best - 1e-7:
                best, wait = v, 0
                best_w = (copy.deepcopy(net.coefs_), copy.deepcopy(net.intercepts_))
            else:
                wait += 1
                if wait >= self.patience:
                    break
        net.coefs_, net.intercepts_ = best_w
        self.net_, self.n_iter_, self.best_val_mse_ = net, ep + 1, best
        return self

    def predict(self, X):
        return self.net_.predict(X[:, :-1])


if __name__ == "__main__":

    # ---------------- metric and seed (command-line) ----------------
    scoring_metrics     = ['neg_root_mean_squared_error', 'neg_mean_absolute_error', 'r2']
    scoring_metrics_tag = ['nrmse', 'nmae', 'r2']
    scoring = dict(zip(scoring_metrics_tag, scoring_metrics))   # all three are recorded

    scrmt = int(sys.argv[1]) if len(sys.argv) > 1 else 0        # metric used to pick the best model
    rst   = int(sys.argv[2]) if len(sys.argv) > 2 else 6        # random seed
    refit_tag = scoring_metrics_tag[scrmt]
    N_JOBS = int(os.environ.get('N_JOBS', -1))

    # ---------------- config ----------------
    fname = "TC_ModelObs_CyclInfo.nc"
    feature_file = "FS_MLPNN_res.txt"
    target = 'res_hs'                       # switch to 'ohs' later (and use its own feature file)
    normet = 1                              # 1: arrays between 0 and 1 (best in the Optuna search)
    out_csv = f"GridSearch_best_params_{refit_tag}_rst{rst:03d}.csv"
    out_npz = f"GridSearch_targv_{refit_tag}_rst{rst:03d}.npz"

    # Grid narrowed with the Optuna results: tanh, batch 128, patience 20 fixed.
    # 5 x 3 x 3 = 45 combinations x 3 inner folds = 135 fits per held-out storm.
    # param_grid = {'hidden_layer_sizes': [(10,), (20,), (50,), (20, 20), (50, 50)],
    #    'alpha': [1e-3, 3e-3, 1e-2],
    #    'learning_rate_init': [3e-4, 1e-3, 3e-3]}

    param_grid = {
        'hidden_layer_sizes': [(2,),(10,),(20,),(50,),(70,),(100,),
            (5,5),(10,10),(20,20),(50,50),(70,70),(100,100)],
        'alpha': [0.001, 0.005, 0.01, 0.1],
        'learning_rate_init': [5e-5, 1e-4, 5e-4, 1e-3],
    }


    # ---------------- read data ----------------
    LVARS = [
        'clat', 'clon', 'vmx', 'rmw', 'vfm', 'cdir_sin', 'cdir_cos',
        'r34', 'r34_ne', 'r34_nw', 'r34_se', 'r34_sw',
        'r50', 'r50_ne', 'r50_nw', 'r50_se', 'r50_sw',
        'r64', 'r64_ne', 'r64_nw', 'r64_se', 'r64_sw']

    OMVARS = [
        'jlday', 'obs_ID', 'cid', 'mlat', 'mlon', 'csec', 'ccdist', 'ccangle',
        'ohs', 'ownd', 'mhs', 'mtp', 'mwnd', 'hhs', 'hwnd',
        'ghs_c', 'ghs_em', 'gwnd_c', 'gwnd_em']

    ds, df = read_tc_dataset(fname, OMVARS, LVARS)
    print(f"Loaded {df.shape[0]} observations x {df.shape[1]} columns")

    # Clean rmw IN df (cleaning the raw dict `ds` after df is built does not change df)
    for c in [c for c in df.columns if c.startswith('rmw_lag')]:
        n_bad = int((df[c] > 300000).sum())
        df.loc[df[c] > 300000, c] = np.nan
        print(f"  {c}: {n_bad} values > 300 km set to NaN")

    # residue and ccangle decomposition
    df['res_hs'] = df['mhs'] - df['ohs']
    df['ccangle_sin'] = np.sin(np.deg2rad(df['ccangle']))
    df['ccangle_cos'] = np.cos(np.deg2rad(df['ccangle']))

    # input variables
    selected_features = read_feature_list(feature_file)
    missing = [f for f in selected_features if f not in df.columns]
    if missing:
        raise ValueError(f"features not in df: {missing}")

    # --- Input / output arrays (2-D) ---
    mip = np.asarray(df[selected_features], dtype='float64')
    mop = np.asarray(df[[target]], dtype='float64')

    # remove rows with NaN/inf; every aligned array gets the same mask
    valid = np.isfinite(mip).all(axis=1) & np.isfinite(mop).all(axis=1)
    print(f"{valid.sum()}/{len(valid)} rows kept")
    mip, mop = mip[valid], mop[valid]
    cid   = df['cid'].values[valid]
    mhs   = df['mhs'].values[valid]       # kept for the final validation
    ohs   = df['ohs'].values[valid]
    otime = df['time'].values[valid]

    # --- Normalization ---
    npinpmin = np.nanmin(mip, axis=0)
    npinpmax = np.nanmax(mip, axis=0)
    const = npinpmax == npinpmin
    if const.any():
        raise ValueError("constant input columns: " +
                         str([selected_features[k] for k in np.where(const)[0]]))
    x1, mipnp1, mipnp2 = normalization(mip, normet, npinpmin, npinpmax)
    npoutmin = np.nanmin(mop, axis=0)
    npoutmax = np.nanmax(mop, axis=0)
    v1, mopnp1, mopnp2 = normalization(mop, normet, npoutmin, npoutmax)
    scale_m = float(np.ravel(mopnp2 - mopnp1)[0])   # normalized RMSE/MAE x scale_m = metres (normet = 1)

    # final target variable for validation (out-of-sample prediction, physical units)
    targv = np.full(len(cid), np.nan, dtype='f')

    # ---------------- leave-one-storm-out ----------------
    ucid = np.unique(cid)
    rows = []
    print(" "); print("-----")
    for i in range(len(ucid)):
        indtrain = np.where(cid != ucid[i])[0]
        indval = np.where(cid == ucid[i])[0]

        # --- Train / validation sets (y 1-D) ---
        # storm code appended as the LAST column (used only by StormESMLP's early stopping)
        gcode = np.unique(cid, return_inverse=True)[1].astype('float64')
        X_train = np.column_stack([x1[indtrain, :], gcode[indtrain]])
        X_test = np.column_stack([x1[indval, :], gcode[indval]])
        y_train = v1[indtrain, 0]
        y_test = v1[indval, 0]
        g_train = cid[indtrain]

        print(" CV cycle "+repr(i)+". Training records: "+repr(X_train.shape[0]))
        print(" CV cycle "+repr(i)+". Validation records: "+repr(X_test.shape[0])+"  "
              +repr(np.round((100*(len(indval)/len(indtrain))),2))+"%")

        # tanh, adam, batch 128 are fixed inside StormESMLP; early stopping on held-out STORMS
        model = StormESMLP(batch_size=128, max_epochs=2000, patience=20,
                           val_storm_frac=0.15, random_state=rst)

        # inner CV grouped by storm (no storm in both train and validation folds)
        grid_search = GridSearchCV(estimator=model, param_grid=param_grid,
                                   scoring=scoring, refit=refit_tag,
                                   cv=GroupKFold(n_splits=3), n_jobs=N_JOBS)
        grid_search.fit(X_train, y_train, groups=g_train)

        # best configuration + its inner-CV scores (all three metrics)
        b = grid_search.best_index_
        cvr = grid_search.cv_results_
        row = dict(storm=ucid[i], n_train=len(indtrain), n_val=len(indval),
                   **grid_search.best_params_,
                   cv_rmse_m=-cvr['mean_test_nrmse'][b] * scale_m,
                   cv_mae_m=-cvr['mean_test_nmae'][b] * scale_m,
                   cv_r2=cvr['mean_test_r2'][b])
        rows.append(row)
        pd.DataFrame(rows).to_csv(out_csv, sep='\t', index=False, float_format='%.6f')  # saved every storm
        print("NN Best Parameters ("+scoring_metrics[scrmt]+") random init "+repr(rst)+" :",
              grid_search.best_params_)

        best_model = grid_search.best_estimator_

        # Predict target variable in the held-out storm (MLPRegressor.predict, 1-D)
        aux = best_model.predict(X_test)
        targv[indval] = denormalization(aux, normet, mopnp1, mopnp2)

        del aux, grid_search, best_model, model

    np.savez(out_npz, targv=targv, cid=cid, mhs=mhs, ohs=ohs, otime=otime)


