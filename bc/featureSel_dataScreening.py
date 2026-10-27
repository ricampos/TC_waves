#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Data visualization, data screening, and feature selection.
"""

import matplotlib
matplotlib.use('Agg')
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler
from sklearn.feature_selection import RFECV
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import GroupKFold
from sklearn.inspection import permutation_importance
from sklearn.ensemble import RandomForestRegressor
from sklearn.svm import SVR
from sklearn.neural_network import MLPRegressor
from xgboost import XGBRegressor
from sklearn.linear_model import ElasticNet
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, WhiteKernel
from lightgbm import LGBMRegressor
from catboost import CatBoostRegressor
import netCDF4 as nc
import warnings; warnings.filterwarnings("ignore")
fnetcdf="NETCDF4"

# Font size and style
sl=13
matplotlib.rcParams.update({'font.size': sl}); plt.rc('font', size=sl) 
matplotlib.rc('xtick', labelsize=sl); matplotlib.rc('ytick', labelsize=sl); matplotlib.rcParams.update({'font.size': sl})

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

def timeseries(dtime,wvar,wname):
    '''
    Time-series plots of model (lines) and observations (discrete points).
    '''
    fig1, ax = plt.subplots(figsize=(13, 4))
    # ax.plot(dtime, wvar, color='dimgray', marker='.', linestyle='', linewidth=2., zorder=2)
    ax.plot(dtime, wvar, marker='.', markersize=5, markerfacecolor='k',
        markeredgewidth=0.5, markeredgecolor='black', linestyle='', linewidth=0., zorder=3)

    ax.set_xlim(dtime.min()-pd.Timedelta(days=10), dtime.max()+pd.Timedelta(days=10))
    ax.xaxis.set_major_formatter(DateFormatter('%b/%y'))
    ax.fmt_xdata = DateFormatter('%b/%y')
    ax.set_xlabel("Time"); ax.set_ylabel(wname, fontsize=15)
    plt.tight_layout()
    plt.grid(c='grey', ls='--', alpha=0.3, zorder=1)
    plt.savefig("TimeSeries_"+wname+".png", dpi=200, facecolor='w', edgecolor='w', orientation='portrait',
        format='png', bbox_inches='tight', pad_inches=0.1)
    plt.close(fig1)


def pdf(wvar,wname):
    '''
    Probability Density Plot
    '''
    from scipy.ndimage.filters import gaussian_filter
    from scipy.stats import gaussian_kde

    # remove NaN
    ind=np.where(wvar>-999.)
    if np.size(ind)>0:
        wvar=wvar[ind[0]]
    else:
        raise ValueError(' No quality data available.')

    # plot
    px=np.linspace(wvar.min(),wvar.max()*0.99,100)

    # plt.close('all')
    fig1 = plt.figure(1,figsize=(7,6)); ax = fig1.add_subplot(111)
    dx = gaussian_kde(wvar)
    ax.fill_between(px, 0., gaussian_filter(dx(px), 1.), color='grey', alpha=0.7,zorder=1)
    ax.plot(px,gaussian_filter(dx(px), 1.),color='k',linestyle='-',linewidth=1.,alpha=0.8,zorder=2)
    plt.grid(c='grey', ls=':', alpha=0.5,zorder=1)
    plt.locator_params(axis='y', nbins=7) ; plt.locator_params(axis='x', nbins=7) 
    ax.set_ylabel("Probability Density")
    ax.set_xlabel(wname)
    plt.tight_layout(); plt.axis('tight')

    plt.savefig("PDF_"+wname+".png", dpi=200, facecolor='w', edgecolor='w', orientation='portrait',
        format='png', bbox_inches='tight', pad_inches=0.1)

    #ax.set_yscale( "log" ); plt.axis('tight');plt.ylim(ymin = 0.001)
    #plt.savefig("PDFl_"+wname+".png", dpi=200, facecolor='w', edgecolor='w', orientation='portrait',
    #    format='png', bbox_inches='tight', pad_inches=0.1)

    plt.close(fig1); del fig1, ax


def plot_corr_clustermap(ndf, out_path='corr_clustermap.png', method='average'):
    corr = ndf.corr(method='pearson')

    g = sns.clustermap(
        corr,
        method=method,          # linkage method for the dendrogram
        cmap='coolwarm',
        vmin=-1, vmax=1,
        center=0,
        figsize=(24, 24),
        xticklabels=True,
        yticklabels=True,
        annot=False,
        cbar_kws={'label': 'Pearson r'},
    )
    g.ax_heatmap.tick_params(labelsize=6)
    g.savefig(out_path, dpi=200, bbox_inches='tight')
    print(f"Saved {out_path}")
    return corr, g


def drop_lags(ndf, keep_lag='lag0', lag_suffixes=('_lag3', '_lag6', '_lag12', '_lag24')):
    """Keep only lag0 (and non-lag) columns, dropping the redundant
    higher-lag duplicates (e.g. clon_lag3/6/12/24 vs clon_lag0)."""
    keep_cols = [c for c in ndf.columns if not c.endswith(lag_suffixes)]
    return ndf[keep_cols]


def top_correlated_pairs(corr, n=30, exclude_self=True):
    """List the n strongest |correlation| pairs — useful for spotting
    near-duplicate columns (e.g. r34 vs r34_ne) worth dropping before RFE."""
    c = corr.copy()
    if exclude_self:
        np.fill_diagonal(c.values, np.nan)
    pairs = (
        c.stack()
        .rename('r')
        .reset_index()
        .rename(columns={'level_0': 'var1', 'level_1': 'var2'})
    )
    pairs['abs_r'] = pairs['r'].abs()
    pairs = pairs[pairs['var1'] < pairs['var2']]  # drop mirrored duplicates
    return pairs.sort_values('abs_r', ascending=False).head(n)


def plot_corr_target_first(nfdf, target_vars=('res_hs', 'ohs'),
                            out_path='corr_matrix_target_first.png',
                            cluster_rest=True):
    """
    Pearson correlation heatmap with target_vars pinned as the first
    rows/columns (highlighted in bold), and the remaining variables
    hierarchically clustered so correlated blocks still group together.
    """
    corr = nfdf.corr(method='pearson')
    target_vars = [v for v in target_vars if v in corr.columns]
    other_vars = [c for c in corr.columns if c not in target_vars]

    if cluster_rest and len(other_vars) > 2:
        from scipy.cluster.hierarchy import linkage, dendrogram
        from scipy.spatial.distance import squareform
        d = 1 - corr.loc[other_vars, other_vars].abs()
        np.fill_diagonal(d.values, 0)
        condensed = squareform(d.values, checks=False)
        Z = linkage(condensed, method='average')
        order = dendrogram(Z, no_plot=True)['leaves']
        other_vars = [other_vars[i] for i in order]

    ordered_cols = target_vars + other_vars
    corr_ordered = corr.loc[ordered_cols, ordered_cols]

    fig, ax = plt.subplots(figsize=(24, 24))
    sns.heatmap(
        corr_ordered, cmap='coolwarm', vmin=-1, vmax=1, center=0,
        xticklabels=True, yticklabels=True, annot=False,
        cbar_kws={'label': 'Pearson r'}, ax=ax,
    )
    ax.tick_params(labelsize=6)

    # bold + enlarge the target variable tick labels
    for label in list(ax.get_xticklabels()) + list(ax.get_yticklabels()):
        if label.get_text() in target_vars:
            label.set_fontweight('bold')
            label.set_fontsize(10)

    # visually separate target rows/cols with a line
    n_target = len(target_vars)
    ax.axhline(n_target, color='black', linewidth=1.5)
    ax.axvline(n_target, color='black', linewidth=1.5)

    fig.savefig(out_path, dpi=200, bbox_inches='tight')
    print(f"Saved {out_path}")
    return corr_ordered
    """Full list of unique variable pairs and their Pearson r, sorted by
    the signed correlation value (ascending: most negative -> most
    positive). No top-N cutoff."""
    c = corr.copy()
    np.fill_diagonal(c.values, np.nan)
    pairs = (
        c.stack()
        .rename('r')
        .reset_index()
        .rename(columns={'level_0': 'var1', 'level_1': 'var2'})
    )
    pairs = pairs[pairs['var1'] < pairs['var2']]  # drop mirrored duplicates
    return pairs.sort_values('r', ascending=ascending).reset_index(drop=True)

def all_correlated_pairs(corr, ascending=True):
    """Full list of unique variable pairs and their Pearson r, sorted by
    the signed correlation value (ascending: most negative -> most
    positive). No top-N cutoff."""
    c = corr.copy()
    np.fill_diagonal(c.values, np.nan)
    pairs = (
        c.stack()
        .rename('r')
        .reset_index()
        .rename(columns={'level_0': 'var1', 'level_1': 'var2'})
    )
    pairs = pairs[pairs['var1'] < pairs['var2']]
    return pairs.sort_values('r', ascending=ascending).reset_index(drop=True)

def save_correlated_pairs(pairs, out_path='correlated_pairs.txt'):
    with open(out_path, 'w') as f:
        f.write(pairs.to_string(index=False))
    print(f"Saved {len(pairs)} pairs to {out_path}")


def corr_with_target(corr, target, n=None):
    """Sorted correlation of every variable with a single target column
    (e.g. 'ohs_lag0' or a bias column) — more actionable than the full matrix."""
    s = corr[target].drop(index=target).sort_values(key=np.abs, ascending=False)
    return s if n is None else s.head(n)


def plot_corr_with_target(nfdf, target, out_prefix=None):
    """
    Horizontal bar plot of every variable's Pearson correlation with
    `target`, sorted by value (not |value|) so negative and positive
    correlations separate cleanly to either side of zero. Also saves the
    full ranked list to a text file (variable, r).
    """
    if out_prefix is None:
        out_prefix = target

    corr = nfdf.corr(method='pearson')
    s = corr[target].drop(index=target).sort_values(ascending=True)

    colors = ['firebrick' if v < 0 else 'steelblue' for v in s.values]

    fig, ax = plt.subplots(figsize=(8, max(6, 0.22 * len(s))))
    ax.barh(s.index, s.values, color=colors, zorder=3)
    ax.axvline(0, color='k', linewidth=0.8)
    ax.set_xlabel(f"Pearson r with {target}")
    ax.set_xlim(-1, 1)
    ax.grid(c='grey', ls='--', alpha=0.3, zorder=1)
    ax.tick_params(labelsize=8)
    plt.tight_layout()
    plt.savefig(f"Corr_{out_prefix}.png", dpi=200, facecolor='w', edgecolor='w',
                bbox_inches='tight', pad_inches=0.1)
    plt.close(fig)

    out_txt = f"corr_with_{out_prefix}.txt"
    with open(out_txt, 'w') as f:
        for var, r in s.sort_values(key=np.abs, ascending=False).items():
            f.write(f"{var}\t{r:.4f}\n")
    print(f"Saved Corr_{out_prefix}.png and {out_txt}")

    return s

def plot_pps_with_target(nfdf, target, out_prefix=None, sample=None, random_seed=1):
    """
    Predictive Power Score (nonlinear, asymmetric) of every variable
    predicting `target`, via ppscore's decision-tree-based estimator.
    Horizontal bar plot (0-1, no sign) + saved ranked text file —
    same structure as plot_corr_with_target, but nonlinear.

    sample: pass an int to subsample rows for speed (ppscore fits a
    decision tree per variable); None uses all rows.
    """
    import ppscore as pps

    if out_prefix is None:
        out_prefix = target

    res = pps.predictors(nfdf, y=target, sample=sample, random_seed=random_seed)
    res = res[res['x'] != target]
    s = res.set_index('x')['ppscore'].sort_values(ascending=True)

    fig, ax = plt.subplots(figsize=(8, max(6, 0.22 * len(s))))
    ax.barh(s.index, s.values, color='seagreen', zorder=3)
    ax.set_xlabel(f"PPS predicting {target}")
    ax.set_xlim(0, 1)
    ax.grid(c='grey', ls='--', alpha=0.3, zorder=1)
    ax.tick_params(labelsize=8)
    plt.tight_layout()
    plt.savefig(f"PPS_{out_prefix}.png", dpi=200, facecolor='w', edgecolor='w',
                bbox_inches='tight', pad_inches=0.1)
    plt.close(fig)

    out_txt = f"pps_with_{out_prefix}.txt"
    with open(out_txt, 'w') as f:
        for var, p in s.sort_values(ascending=False).items():
            f.write(f"{var}\t{p:.4f}\n")
    print(f"Saved PPS_{out_prefix}.png and {out_txt}")

    return s


def plot_pps_matrix(nfdf, out_path='pps_matrix.png', sample=None, random_seed=1):
    """
    Full pairwise PPS matrix (all x predicting all y). Unlike Pearson,
    this is ASYMMETRIC (x predicting y != y predicting x), so no
    clustering/symmetric reordering is applied here — rows/cols kept in
    nfdf's original column order for a like-for-like comparison.

    Expensive: fits one decision tree per (x, y) pair -> n_vars^2 fits.
    Use `sample` to subsample rows if this is too slow for your n_vars.
    """
    import ppscore as pps

    matrix_df = pps.matrix(nfdf, sample=sample, random_seed=random_seed)
    pivot = matrix_df.pivot(columns='x', index='y', values='ppscore')
    pivot = pivot.loc[nfdf.columns, nfdf.columns]

    fig, ax = plt.subplots(figsize=(24, 24))
    sns.heatmap(
        pivot, cmap='viridis', vmin=0, vmax=1,
        xticklabels=True, yticklabels=True, annot=False,
        cbar_kws={'label': 'PPS (x predicts y)'}, ax=ax,
    )
    ax.set_xlabel('x (predictor)'); ax.set_ylabel('y (target)')
    ax.tick_params(labelsize=6)
    fig.savefig(out_path, dpi=200, bbox_inches='tight')
    print(f"Saved {out_path}")
    return pivot


def run_rfecv_linear(X, y, groups, ttag, step=1, cv_splits=3):
    valid = X.notna().all(axis=1) & y.notna()
    Xc, yc = X.loc[valid], y.loc[valid]
    gc = np.asarray(groups)[valid.values]
    print(f"[{ttag}] {valid.sum()}/{len(valid)} rows valid (no NaN) for RFECV")

    model = LinearRegression()
    cv = GroupKFold(n_splits=cv_splits)
    rfe = RFECV(estimator=model, step=step, cv=cv)
    rfe.fit(Xc.values, yc.values, groups=gc)

    # Selected: refit on the final subset, rank by |coefficient| descending
    sel_mask = rfe.support_
    final_model = LinearRegression().fit(Xc.values[:, sel_mask], yc.values)
    selected_df = pd.DataFrame({
        'feature': X.columns[sel_mask],
        'abs_coef': np.abs(final_model.coef_),
    }).sort_values('abs_coef', ascending=False).reset_index(drop=True)

    # Excluded: rank 2 = last eliminated (most important, "almost included")
    # ... higher rank = eliminated earlier = less important. Ascending sort
    # puts the most important excluded feature on top
    excluded_df = pd.DataFrame({
        'feature': X.columns[~sel_mask],
        'elimination_rank': rfe.ranking_[~sel_mask],
    }).sort_values('elimination_rank', ascending=True).reset_index(drop=True)

    print(f"[{ttag}] Selected ({len(selected_df)}):", list(selected_df['feature']))
    print(f"[{ttag}] Excluded ({len(excluded_df)}):", list(excluded_df['feature']))

    selected_df.to_csv(f"RFECV_LinearRegression_Selected_{ttag}.txt", sep='\t', index=False, float_format='%.6f')
    excluded_df.to_csv(f"RFECV_LinearRegression_Excluded_{ttag}.txt", sep='\t', index=False, float_format='%.6f')

    return rfe, selected_df, excluded_df

MODEL_FACTORIES = {
    'RandomForest': lambda seed: RandomForestRegressor(
        n_estimators=200, max_depth=None, n_jobs=-1, random_state=seed),
    'SVM': lambda seed: SVR(kernel='rbf', C=1.0, epsilon=0.1),
    'MLP': lambda seed: MLPRegressor(
        hidden_layer_sizes=(100, 100), activation='tanh', solver='adam',
        learning_rate='adaptive', learning_rate_init=1e-3, max_iter=2000,
        random_state=seed, early_stopping=True, n_iter_no_change=15),
    'XGBoost': lambda seed: XGBRegressor(
        n_estimators=300, max_depth=5, learning_rate=0.05, random_state=seed),
    'LightGBM': lambda seed: LGBMRegressor(
        n_estimators=300, max_depth=5, learning_rate=0.05,
        random_state=seed, verbose=-1),
    'CatBoost': lambda seed: CatBoostRegressor(
        iterations=300, depth=5, learning_rate=0.05,
        random_state=seed, verbose=False),
    'ElasticNet': lambda seed: ElasticNet(
        alpha=0.01, l1_ratio=0.5, random_state=seed, max_iter=5000),
    'GaussianProcess': lambda seed: GaussianProcessRegressor(
        kernel=RBF(length_scale=1.0) + WhiteKernel(noise_level=1.0),
        random_state=seed, normalize_y=True),
}

# Rows to subsample for GP specifically -- exact GP is O(n^3), infeasible at full size
GP_MAX_SAMPLES = 3000

def rank_features_permutation(X, y, groups, model_name, ttag,
                               seeds=(1, 4, 7, 10, 13), n_repeats=5, cv_splits=3):
    valid = X.notna().all(axis=1) & y.notna()
    Xc, yc = X.loc[valid].values, y.loc[valid].values
    gc = np.asarray(groups)[valid.values]

    gkf = GroupKFold(n_splits=cv_splits)
    factory = MODEL_FACTORIES[model_name]

    scores = np.zeros((len(seeds), X.shape[1]))
    for i, seed in enumerate(seeds):
        fold_scores = np.zeros((cv_splits, X.shape[1]))
        for k, (tr, te) in enumerate(gkf.split(Xc, yc, groups=gc)):
            tr_idx = tr
            if model_name == 'GaussianProcess' and len(tr) > GP_MAX_SAMPLES:
                rng = np.random.RandomState(seed)
                tr_idx = rng.choice(tr, size=GP_MAX_SAMPLES, replace=False)

            model = factory(seed)
            model.fit(Xc[tr_idx], yc[tr_idx])
            pi = permutation_importance(model, Xc[te], yc[te], n_repeats=n_repeats,
                                         random_state=seed, n_jobs=-1)  # n_jobs=1: avoid joblib disk temp files
            fold_scores[k] = pi.importances_mean
        scores[i] = fold_scores.mean(axis=0)
        print(f"[{ttag}] {model_name}: seed {seed} done")

    mean_imp = scores.mean(axis=0)
    std_imp = scores.std(axis=0)
    lower_bound = mean_imp - std_imp

    ranked = pd.DataFrame({'feature': X.columns, 'importance': mean_imp, 'importance_std': std_imp})

    selected_df = (ranked[lower_bound > 0]
                   .sort_values('importance', ascending=False).reset_index(drop=True))
    excluded_df = (ranked[lower_bound <= 0]
                   .sort_values('importance', ascending=False).reset_index(drop=True))

    selected_df.to_csv(f"FeatureSelected_{model_name}_{ttag}.txt", sep='\t', index=False, float_format='%.6f')
    excluded_df.to_csv(f"FeatureExcluded_{model_name}_{ttag}.txt", sep='\t', index=False, float_format='%.6f')

    print(f"[{ttag}] {model_name}: {len(selected_df)} selected, {len(excluded_df)} excluded")

    return ranked, selected_df, excluded_df

if __name__ == "__main__":

    # Read Obs and Model
    fname="TC_ModelObs_CyclInfo.nc"

    # Variables stored on ('tlag', 'time')
    LVARS = [
        'clat', 'clon', 'vmx', 'rmw', 'vfm', 'cdir_sin', 'cdir_cos',
        'r34', 'r34_ne', 'r34_nw', 'r34_se', 'r34_sw',
        'r50', 'r50_ne', 'r50_nw', 'r50_se', 'r50_sw',
        'r64', 'r64_ne', 'r64_nw', 'r64_se', 'r64_sw']
     
    # Variables stored on ('time',)
    OMVARS = [
        'jlday', 'obs_ID', 'cid', 'mlat', 'mlon', 'csec', 'ccdist', 'ccangle',
        'ohs', 'ownd', 'mhs', 'mtp', 'mwnd', 'hhs', 'hwnd',
        'ghs_c', 'ghs_em', 'gwnd_c', 'gwnd_em']
 
    ds, df = read_tc_dataset(fname,OMVARS,LVARS)
    print(f"Loaded {df.shape[0]} observations x {df.shape[1]} columns")
    print(df.columns.tolist())

    # Clean
    ds['rmw'][ds['rmw']>300000]=np.nan

    # Time series and PDF plots
    for v in list(ds.keys())[8::]:
        if len(ds[v].shape)==1:
            timeseries(ds['time_dt'],ds[v],v)
            pdf(ds[v],v)
        else:
            for i in range(0,ds[v].shape[0]):
                timeseries(ds['time_dt'],ds[v][i,:],str(v)+"_Lag"+str(ds['tlag'][i]))
                pdf(ds[v][i,:],str(v)+"_Lag"+str(ds['tlag'][i]))

    # Add residue to the dataframe
    df['res_hs'] = df['mhs'] - df['ohs']
    # print all the variables of df
    with open('df_variables.txt', 'w') as f:
        for col in df.columns:
            f.write(f"{col}\n")

    print(' '); print(df.columns.tolist()); print(' ')

    # Keep variables used for modeling.
    exclude = ['obs_ID', 'time', 'cid']
    cols = [c for c in df.columns if c not in exclude and not c.startswith(('g', 'h'))]
    fdf = df[cols]
    print(' '); print(fdf.columns.tolist()); print(' ')

    # Normalization
    scaler = MinMaxScaler(feature_range=(0, 1))
    nfdf = pd.DataFrame(scaler.fit_transform(fdf[cols]), columns=cols, index=fdf.index)
    # plot again
    # for v in list(nfdf.keys()):
    #    timeseries(ds['time_dt'],nfdf[v],v)
    #    pdf(nfdf[v],v)

    # Cross-correlation
    # (Pearson)
    corr_ordered = plot_corr_target_first(nfdf, target_vars=('res_hs', 'ohs'))
    pairs = all_correlated_pairs(corr_ordered, ascending=True)
    save_correlated_pairs(pairs, out_path='correlated_pairs.txt')
    # Saved 7381 pairs to correlated_pairs.txt
    plot_corr_with_target(nfdf, 'res_hs')
    plot_corr_with_target(nfdf, 'ohs')
    # len(list(nfdf.keys())) = 122
    # Predictive Power Score (nonlinear) — targeted, same as Pearson bar plots
    plot_pps_with_target(nfdf, 'res_hs', sample=10000)
    plot_pps_with_target(nfdf, 'ohs', sample=10000)

    # Organize arrays
    wtrgtobs = nfdf[['ohs', 'ownd', 'res_hs']] # isolate the obs and target variables
    winput = nfdf.drop(columns=['ohs', 'ownd', 'res_hs']) # input variables only

    # === Recursive feature elimination (RFECV) ===
    # storm id per row, aligned to winput/wtrgtobs index (nfdf's index == df's index)
    groups = df.loc[nfdf.index, 'cid'].values # cross-validation through cyclone id (each storm)
    # 40441/45848 rows valid (no NaN) for RFECV
    # Linear model
    rfe_hs,  sel_hs,  exc_hs  = run_rfecv_linear(winput, wtrgtobs['res_hs'], groups, ttag='res_hs')
    rfe_ohs, sel_ohs, exc_ohs = run_rfecv_linear(winput, wtrgtobs['ohs'],   groups, ttag='ohs')
    # Expanding to non-linear models
    for model_name in ['RandomForest', 'SVM', 'MLP', 'XGBoost', 'LightGBM', 'CatBoost', 'ElasticNet', 'GaussianProcess']:
        rank_features_permutation(winput, wtrgtobs['res_hs'], groups, model_name, ttag='res_hs')
        rank_features_permutation(winput, wtrgtobs['ohs'],   groups, model_name, ttag='ohs')


