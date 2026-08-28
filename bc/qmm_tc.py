#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Bias correction of wave forecast for TCs
Cross Validation through the events
PModel, HAFS, and GEFS
"""

import numpy as np
import pandas as pd
import os
import netCDF4 as nc
import xarray as xr
import time
from calendar import timegm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy import stats
import mvalstats
import pvalstats
from pvalstats import *
import sys
import warnings; warnings.filterwarnings("ignore")

import matplotlib.colors as colors
sl=14 # plot style configuration
matplotlib.rcParams.update({'font.size': sl}); plt.rc('font', size=sl) 
matplotlib.rc('xtick', labelsize=sl); matplotlib.rc('ytick', labelsize=sl); matplotlib.rcParams.update({'font.size': sl})


# Quantile Mapping bias-correction
def qm_train(model=None,obs=None,prob=np.nan):
    """ 
    Univariate linear regression calibration using the Quantile Mapping Method.
    Based on:
     https://doi.org/10.3390/rs14194918
     https://doi.org/10.1080/01621459.2014.929522
     https://doi.org/10.1016/j.jhydrol.2010.10.024

    Fit module
    Input:
    - model data
    - observations ("truth")
    - probability array (optional) to define the percentiles 
    Output:
    - slope
    - intercept
    """

    if (np.size(model)>2)==False or (np.size(obs)>2)==False:
        raise ValueError("Model and Obs arrays must be informed.")

    if np.size(model) != np.size(obs):
        raise ValueError("Model and Obs arrays must have the same sizes.")

    if size(prob)==1:
        prob=np.array(np.arange(0.5,99.5+0.1,0.5),'f')  # Probability array

    slope,intercept = np.polyfit(np.nanpercentile(model,prob),np.nanpercentile(obs,prob),1)
    # model_cal=np.array((model*slope)+intercept)
    print(" QMM linear regression. Slope: "+repr(slope)+"  Intercept: "+repr(intercept))

    return float(slope),float(intercept)


def qmcal(model=None,slope=1.,intercept=0.,pprint='yes'):
    """ 
    Univariate linear regression calibration using the Quantile Mapping Method
    Calibration module based on previously trained qm_train
    Input:
    - model data
    - slope
    - intercept
    Output:
    - calibrated model data
    """

    if (np.size(model)>2)==False:
        raise ValueError("Model array must be informed.")

    model_cal=np.array((model*slope)+intercept)
    if pprint=='yes':
        print(" QMM linear regression. Slope: "+repr(slope)+"  Intercept: "+repr(intercept))

    return np.array(model_cal).astype('float')


def scdensityplot(obs,model,vmin,vmax,maxdiff,lobs,lmodel,nlabel):
    """ 
    Scatter Plot (density)
    """
    aux=np.linspace(0,vmax,100)

    fig1 = plt.figure(1, figsize=(7,6)); ax = fig1.add_subplot(111)
    ax.plot(aux, aux, 'k--', linewidth=1., alpha=0.9, zorder=4)
    hb = ax.hexbin(obs, model, gridsize=100, cmap='jet', mincnt=1,
                   bins='log', zorder=3)

    ax.locator_params(axis='y', nbins=7)
    ax.locator_params(axis='x', nbins=7)
    ax.set_xlabel(lobs); ax.set_ylabel(lmodel)
    ax.grid(c='grey', ls=':', alpha=0.5, zorder=1)
    ax.set_xlim(xmin=0.9, xmax=aux.max())
    ax.set_ylim(ymin=0.9, ymax=aux.max())

    merr = mvalstats.metrics(model,obs,vmin=vmin,vmax=vmax,maxdiff=maxdiff, pctlerr='no')
    a, b, cc, p_value, std_err = stats.linregress(obs,model)
    stats_text = (f"Bias = {merr[0]:.2f}\n"
                  f"RMSE = {merr[1]:.2f}\n"
                  f"SI = {merr[5]:.2f}\n"
                  f"CC = {merr[7]:.2f}\n"
                  f"Slope = {a:.2f}\n"
                  f"Intercept = {b:.2f}")

    ax.text(0.97, 0.03, stats_text,
            transform=ax.transAxes,
            fontsize=15,
            va='bottom', ha='right',
            bbox=dict(boxstyle='round', facecolor='white',
                      edgecolor='grey', alpha=0.85),
            zorder=5)

    plt.tight_layout()
    plt.savefig("ScatterPlot_"+nlabel+".png", dpi=200, facecolor='w', edgecolor='w',orientation='portrait', format='png',transparent=False, bbox_inches='tight', pad_inches=0.1)
    plt.close('all')


if __name__ == "__main__":

    # Read Obs and Model data
    dirdata="/home/ricardo/work/noaa/analysis/TC_Waves/4postproc/data/"
    # Headers
    hd = 'mean, variance, skewness, kurtosis, min, max, percentile80, percentile90, percentile95, percentile99, percentile99.9'
    merrname=["bias","RMSE","NBias","NRMSE","SCrmse","SI","HH","CC","N","bias_p95","RMSE_p95","N_p95"]

    # HAFS
    dh = pd.read_csv(dirdata+"Data_HAFS_fromObs.txt", sep='\t')
    hhs = np.array(dh['HTSGW_surface'])
    hwnd = np.sqrt( np.array(dh['UGRD_surface'])**2 + np.array(dh['VGRD_surface'])**2 )
    del dh

    # GEFS
    for i in range(0,31):
        dg = pd.read_csv(dirdata+"Data_GEFS_fromObs_m"+str(i).zfill(2)+".txt", sep='\t')
        if i==0:
            ghs = np.array([dg['HTSGW_surface']])
            gwnd = np.array([np.sqrt( np.array(dg['UGRD_surface'])**2 + np.array(dg['VGRD_surface'])**2 )])
        else:
            ghs = np.append(ghs, np.array([dg['HTSGW_surface']]), axis=0)
            gwnd = np.append(gwnd, np.array([np.sqrt( np.array(dg['UGRD_surface'])**2 + np.array(dg['VGRD_surface'])**2 )]), axis=0)

    # PModel
    df = pd.read_csv(dirdata+"Data_Obs_PModel_Default.txt", sep='\t')
    ot = df.time.values
    ot_dt = pd.to_datetime(ot.astype(str), format='%Y%m%d%H%M') # datetime
    ohs = np.array(df['obs_hs']); ownd = np.array(df['obs_wnd'])
    mhs = np.array(df['pm_hs']); mwnd = np.sqrt( np.array(df['pm_uwnd'])**2 + np.array(df['pm_vwnd'])**2 )
    oid = np.array(df['id']).astype('str')
    cmap = np.array(df['cmap']); csec = np.array(df['csec'])
    ccdist = np.array(df['ccdist']); ccangle = np.array(df['ccangle'])
    cid = np.array(df['cid'])

    ind = np.where( (np.mean(ghs,axis=0)>1.0) & (mhs>1.0) & (ohs>1.0) & (hhs>1.0) & (cmap==5) &
        (np.char.startswith(oid, 'DWSD')!=True) &
        (np.char.startswith(oid, 'MICROSWIFT')!=True) &
        (np.char.startswith(oid, 'WSRA')!=True) &
        (np.char.startswith(oid, 'SAILDRONE')!=True) )[0]

    ot = ot[ind]; ot_dt = ot_dt[ind]
    ohs = ohs[ind]; ownd = ownd[ind]
    mhs = mhs[ind]; mwnd = mwnd[ind]
    oid = oid[ind]; cid = cid[ind]
    csec = csec[ind]; ccdist = ccdist[ind]; ccangle = ccangle[ind]
    hhs = hhs[ind]; hwnd = hwnd[ind]
    ghs = ghs[:,ind]; gwnd = gwnd[:,ind]
    del ind

    fcid=[]
    for i in np.unique(cid):
        ind=np.where(cid==i)[0]
        if len(ind)>=10:
            fcid=np.append(fcid,int(i))

        del ind

    fcid = np.array(fcid).astype('int') # Total 84 cyclones with valid obs
    ind = np.where(np.isin(cid, fcid)==True)[0]
    ot = ot[ind]; ot_dt = ot_dt[ind]
    ohs = ohs[ind]; ownd = ownd[ind]
    mhs = mhs[ind]; mwnd = mwnd[ind]
    oid = oid[ind]; cid = cid[ind]
    csec = csec[ind]; ccdist = ccdist[ind]; ccangle = ccangle[ind]
    hhs = hhs[ind]; hwnd = hwnd[ind]
    ghs = ghs[:,ind]; gwnd = gwnd[:,ind]
    del ind

    # QMM Bias Correction and Cross Validation (84 cycles) throughout the events

    # --- PModel ---
    mhs_bc = np.copy(mhs)
    fsl = np.zeros(len(mhs))*np.nan
    fitc = np.zeros(len(mhs))*np.nan
    c=0
    for i in fcid:
        indtrain = np.where(cid!=int(i))[0]
        indval = np.where(cid==int(i))[0]

        sl,itc = qm_train(model=mhs[indtrain],obs=ohs[indtrain])
        mhs_bc[indval] = qmcal(model=mhs[indval],slope=sl,intercept=itc,pprint='no')
        fsl[c] = float(sl); fitc[c] = float(itc)
        c=c+1

    # Plots
    mop=ModelObsPlot(model=[mhs,mhs_bc],obs=ohs, mlabels=['PModel','PModel_BC'], axisnames=['PModel','Obs'],vaxisname="Obs and PModel",ftag="QMM_PModel_")
    mop.qqplot(); mop.taylordiagram(); mop.scatterplot()
    mop=ModelObsPlot(model=[mhs,mhs_bc],obs=ohs, mlabels=['PModel','PModel_BC'], axisnames=['PModel','Obs'],vaxisname="Obs (shaded) and PModel",ftag="QMM_PModel_")
    mop.pdf()

    scdensityplot(ohs,mhs,1,25.,20.,"Obs","PModel","PModel")
    scdensityplot(ohs,mhs_bc,1,25.,20.,"Obs","PModel","PModel_BC")

    # Table stat metrics
    merr_hs = mvalstats.metrics(mhs,ohs,vmin=1.,vmax=25.,maxdiff=20., pctlerr='yes')
    dst = pd.DataFrame(np.array([merr_hs]).round(4), columns=merrname)
    dst.to_csv("Merr_Hs_PModel.txt", index=False, header=True, sep='\t'); del dst, merr_hs

    merr_hs = mvalstats.metrics(mhs_bc,ohs,vmin=1.,vmax=25.,maxdiff=20., pctlerr='yes')
    dst = pd.DataFrame(np.array([merr_hs]).round(4), columns=merrname)
    dst.to_csv("Merr_Hs_PModel_BC.txt", index=False, header=True, sep='\t'); del dst, merr_hs


    # --- HAFS ---
    hhs_bc = np.copy(hhs)
    fsl = np.zeros(len(hhs))*np.nan
    fitc = np.zeros(len(hhs))*np.nan
    c=0
    for i in fcid:
        indtrain = np.where(cid!=int(i))[0]
        indval = np.where(cid==int(i))[0]

        sl,itc = qm_train(model=hhs[indtrain],obs=ohs[indtrain])
        hhs_bc[indval] = qmcal(model=hhs[indval],slope=sl,intercept=itc,pprint='no')
        fsl[c] = float(sl); fitc[c] = float(itc)
        c=c+1

    # Plots
    mop=ModelObsPlot(model=[hhs,hhs_bc],obs=ohs, mlabels=['HAFS','HAFS_BC'], axisnames=['HAFS','Obs'],vaxisname="Obs and HAFS",ftag="QMM_HAFS_")
    mop.qqplot(); mop.taylordiagram(); mop.scatterplot()
    mop=ModelObsPlot(model=[hhs,hhs_bc],obs=ohs, mlabels=['HAFS','HAFS_BC'], axisnames=['HAFS','Obs'],vaxisname="Obs (shaded) and HAFS",ftag="QMM_HAFS_")
    mop.pdf()

    scdensityplot(ohs,hhs,1,25.,20.,"Obs","HAFS","HAFS")
    scdensityplot(ohs,hhs_bc,1,25.,20.,"Obs","HAFS","HAFS_BC")

    # Table stat metrics
    merr_hs = mvalstats.metrics(hhs,ohs,vmin=1.,vmax=25.,maxdiff=20., pctlerr='yes')
    dst = pd.DataFrame(np.array([merr_hs]).round(4), columns=merrname)
    dst.to_csv("Merr_Hs_HAFS.txt", index=False, header=True, sep='\t'); del dst, merr_hs

    merr_hs = mvalstats.metrics(hhs_bc,ohs,vmin=1.,vmax=25.,maxdiff=20., pctlerr='yes')
    dst = pd.DataFrame(np.array([merr_hs]).round(4), columns=merrname)
    dst.to_csv("Merr_Hs_HAFS_BC.txt", index=False, header=True, sep='\t'); del dst, merr_hs


    # --- GEFSc ---
    ghs= ghs[0,:]
    ghs_bc = np.copy(ghs)
    fsl = np.zeros(len(ghs))*np.nan
    fitc = np.zeros(len(ghs))*np.nan
    c=0
    for i in fcid:
        indtrain = np.where(cid!=int(i))[0]
        indval = np.where(cid==int(i))[0]

        sl,itc = qm_train(model=ghs[indtrain],obs=ohs[indtrain])
        ghs_bc[indval] = qmcal(model=ghs[indval],slope=sl,intercept=itc,pprint='no')
        fsl[c] = float(sl); fitc[c] = float(itc)
        c=c+1

    # Plots
    mop=ModelObsPlot(model=[ghs,ghs_bc],obs=ohs, mlabels=['GEFSc','GEFSc_BC'], axisnames=['GEFSc','Obs'],vaxisname="Obs and GEFSc",ftag="QMM_GEFSc_")
    mop.qqplot(); mop.taylordiagram(); mop.scatterplot()
    mop=ModelObsPlot(model=[ghs,ghs_bc],obs=ohs, mlabels=['GEFSc','GEFSc_BC'], axisnames=['GEFSc','Obs'],vaxisname="Obs (shaded) and GEFSc",ftag="QMM_GEFSc_")
    mop.pdf()

    scdensityplot(ohs,ghs,1,25.,20.,"Obs","GEFSc","GEFSc")
    scdensityplot(ohs,ghs_bc,1,25.,20.,"Obs","GEFSc","GEFSc_BC")

    # Table stat metrics
    merr_hs = mvalstats.metrics(ghs,ohs,vmin=1.,vmax=25.,maxdiff=20., pctlerr='yes')
    dst = pd.DataFrame(np.array([merr_hs]).round(4), columns=merrname)
    dst.to_csv("Merr_Hs_GEFSc.txt", index=False, header=True, sep='\t'); del dst, merr_hs

    merr_hs = mvalstats.metrics(ghs_bc,ohs,vmin=1.,vmax=25.,maxdiff=20., pctlerr='yes')
    dst = pd.DataFrame(np.array([merr_hs]).round(4), columns=merrname)
    dst.to_csv("Merr_Hs_GEFSc_BC.txt", index=False, header=True, sep='\t'); del dst, merr_hs


