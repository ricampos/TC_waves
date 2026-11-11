#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
optimize_NNMLP.py
MLP-NN training algorithm using scikit-learn with multiple input hyperparameters
Used for optimization after multiple runs.
Validation is separated.
"""

import sys
import matplotlib
# matplotlib.use('Agg')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.neural_network import MLPRegressor
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


if __name__ == "__main__":

    ftag="MLPNN_res_results"
    fpath="/work/noaa/marine/ricardo.campos/work/analysis/TC_waves/4postproc/AI_ML"
    # Read Obs and Model
    fname="/work/noaa/marine/ricardo.campos/work/analysis/TC_waves/4postproc/AI_ML/TC_ModelObs_CyclInfo.nc"

    activation = ['tanh', 'relu']
    hidden_layer_sizes = [ (2,),(10,),(20,),(50,),(80,),(120,),
            (10,10),(20,20),(50,50),(80,80),(120,120)]

    alpha = [0.001, 0.005, 0.01, 0.1] # L2 regularization strength
    learning_rate_init = [1e-5, 5e-5, 1e-4, 5e-4, 1e-3]

    # Input Arguments
    actv = 0 # activation function
    rst = int(sys.argv[1]) # random initialization
    hl = int(sys.argv[2]) # hidden layer choice
    alp = int(sys.argv[3]) # alpha index
    lri = int(sys.argv[4]) # alpha index

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

    # Add residue to the dataframe
    df['res_hs'] = df['ohs'] - df['mhs'] 
    # print all the variables of df
    with open('df_variables.txt', 'w') as f:
        for col in df.columns:
            f.write(f"{col}\n")

    print(' '); print(df.columns.tolist()); print(' ')

    # decompose ccangle in sin cos
    df['ccangle_sin'] = np.sin(np.deg2rad(df['ccangle']))
    df['ccangle_cos'] = np.cos(np.deg2rad(df['ccangle']))

    # read input variables
    selected_features = read_feature_list("FS_MLPNN_res.txt")

    # --- Input array ---
    mip = np.array(df[selected_features])
    # --- Output Space ---
    mop = np.array([df['res_hs']]).T
    # remove NaN
    valid = np.isfinite(mip).all(axis=1) & np.isfinite(mop).all(axis=1)
    ohs = df['ohs'][valid]; mhs = df['mhs'][valid]; res_hs = df['res_hs'][valid]
    mip, mop = mip[valid], mop[valid]
    cid = df['cid'][valid]

    # --- Normalization ---
    normet = 1 # arrays between 0 and 1.
    # -- Input --
    # npinpmin = np.array(np.min(mip,axis=0)-np.std(mip,axis=0))
    npinpmin = np.array(np.nanmin(mip,axis=0))
    # npinpmax = np.array(np.max(mip,axis=0)+np.std(mip,axis=0))
    npinpmax = np.array(np.nanmax(mip,axis=0))
    x1,mipnp1,mipnp2 = normalization(mip,normet,npinpmin,npinpmax)
    # -- Output --
    npoutmin = np.array(np.nanmin(mop,axis=0))
    # npoutmax = np.array(np.max(mop,axis=0)+np.std(mop,axis=0))
    npoutmax = np.array(np.nanmax(mop,axis=0))
    v1,mopnp1,mopnp2 = normalization(mop,normet,npoutmin,npoutmax)

    # final target variable for validation
    targv = np.zeros(len(cid),'f')*np.nan

    # Cross validation through cyclone ID. Trains with all the cyclones and select 1 for test.
    for i in range(0,len(np.unique(cid))):
        indtrain = np.where(cid != np.unique(cid)[i])[0]
        indval = np.where(cid == np.unique(cid)[i])[0]

        # --- Train/Test and Independent Validation set ---
        X_train = np.array(x1[indtrain,:])
        X_test = np.array(x1[indval,:])
        y_train = np.array(v1[indtrain,0])
        y_test = np.array(v1[indval,0])
        print(" ")
        print(" CV cycle "+repr(i)+". Training records: "+repr(X_train.shape[0])+". Validation records: "+repr(X_test.shape[0])+"  "+repr(np.round((100*(len(indval)/len(indtrain))),2))+"%")

        nn = MLPRegressor(hidden_layer_sizes=hidden_layer_sizes[hl], activation=activation[actv], \
                solver='adam', alpha=alpha[alp], batch_size='auto', learning_rate='adaptive', \
                learning_rate_init=learning_rate_init[lri], max_iter=8000, shuffle=True, \
                random_state=rst, tol=1e-7, verbose=False, warm_start=False, early_stopping=True)

        nn.fit(X_train, y_train)
        # Predict targert variable in the validation set
        aux = nn.predict(X_test)
        targv[indval] = denormalization(aux,normet,mopnp1,mopnp2)
        del aux, nn, X_train, X_test, y_train, y_test, indtrain, indval
        print(" Ok CV cycle "+repr(i)); print(" "); print("---------")

    # --- Save results ---
    ftag=ftag+"_"+str(hl).zfill(2)+"_"+str(actv).zfill(2)+"_"+str(alp).zfill(2)+"_"+str(lri).zfill(2)+"_"+str(rst).zfill(2)
    np.savetxt(fpath + "/" + ftag + ".txt",
               np.column_stack([ohs, mhs, res_hs, targv]),
               fmt='%.4f', delimiter='\t',
               header='ohs\tmhs\tres_hs\ttargv', comments='')

    print(" "); print(" Results saved. Done") 

