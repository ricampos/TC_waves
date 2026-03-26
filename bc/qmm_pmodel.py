#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

Bias correction of PModel outputs

"""

import numpy as np
import pandas as pd
import os
import netCDF4 as nc
import xarray as xr
import time
from calendar import timegm
import mvalstats
import pvalstats
from pvalstats import *
import sys
import warnings; warnings.filterwarnings("ignore")

import matplotlib.colors as colors
sl=13 # plot style configuration
matplotlib.rcParams.update({'font.size': sl}); plt.rc('font', size=sl) 
matplotlib.rc('xtick', labelsize=sl); matplotlib.rc('ytick', labelsize=sl); matplotlib.rcParams.update({'font.size': sl})


# Quantile Mapping bias-correction
def qm_train(model=None,obs=None,prob=np.nan):
    """ 
    Univariate linear regression calibration using the Quantile Mapping Method.
    Base on:
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


if __name__ == "__main__":

    pmrun="Default"

    # Read Obs
    df = pd.read_csv('Data_Obs_PModel_Default.txt', sep='\t')
    ot = df.time.values
    ohs = np.array(df['obs_hs']); mhs = np.array(df['pm_hs'])
    oid = np.array(df['id']).astype('str')
    cmap = np.array(df['cmap'])

    # Subsample certain obs
    # idx = np.where(np.char.startswith(oid, 'NDBC') | np.char.startswith(oid, 'CDIP') | np.char.startswith(oid, 'SPOT')
    # | np.char.startswith(oid, 'DWSD')  )[0]
    # ohs=ohs[idx]; mhs=mhs[idx]; cmap=cmap[idx]; ot=ot[idx]

    # idx = np.where( (np.char.startswith(oid, 'DWSD')!=True) | (np.char.startswith(oid, 'MICROSWIFT')!=True) | (np.char.startswith(oid, 'SAILDRONE')!=True))
    # oid=oid[idx]; ohs=ohs[idx]; mhs=mhs[idx]; cmap=cmap[idx]; ot=ot[idx]

    # Headers
    hd = 'mean, variance, skewness, kurtosis, min, max, percentile80, percentile90, percentile95, percentile99, percentile99.9'
    merrname=["bias","RMSE","NBias","NRMSE","SCrmse","SI","HH","CC","N","bias_p95","RMSE_p95","N_p95"]

    # loop through the cathegories
    acmap=np.array([2,3,4,5]); ncmap=np.array(["Disturbance","ExtrCycl","Subtrop","TropCycl"]).astype('str')
    

    for i in range(0,len(acmap)):

        # ---- Hs ----
        ind = np.where( (ohs>1.) & (mhs>1.) & (cmap==acmap[i]) )
        if np.size(ind)>0:
            ind=ind[0]

        cohs = ohs[ind]; cmhs = mhs[ind]; cot = ot[ind]
        cmhs_bc = np.copy(cmhs)*np.nan
        yrs = cot // 100_000_000; years = np.unique(yrs)

        fsl = np.zeros(len(years),'f')*np.nan; fitc = np.zeros(len(years),'f')*np.nan
        for j in range(0,len(years)):
            indval = np.where(yrs == years[j]); indtrain = np.where(yrs != years[j])
            sl,itc = qm_train(model=cmhs[indtrain],obs=cohs[indtrain])
            cmhs_bc[indval] = qmcal(model=cmhs[indval],slope=sl,intercept=itc,pprint='no')
            fsl[j] = float(sl); fitc[j] = float(itc)
            del indval, indtrain, sl, itc

        # Initial summary stats
        fname = "SlopeIntercepts_"+ncmap[i]+".txt"
        np.savetxt(fname, [fsl, fitc], delimiter='\t', fmt='%.8f')
        del fname

        # Initial summary stats
        fname = "SummaryStats_HsObs_"+ncmap[i]+".txt"
        result=np.round(np.array([mvalstats.smrstat(cohs,0.1,30.)]),3)
        ifile = open(fname,'w'); ifile.write(hd+' \n')
        np.savetxt(ifile,result.astype('str'),fmt="%s",delimiter='\t') 
        ifile.close(); del ifile, fname, result

        fname = "SummaryStats_HsModel_"+ncmap[i]+".txt"
        result=np.round(np.array([mvalstats.smrstat(cmhs,0.1,30.)]),3)
        ifile = open(fname,'w'); ifile.write(hd+' \n')
        np.savetxt(ifile,result.astype('str'),fmt="%s",delimiter='\t') 
        ifile.close(); del ifile, fname, result

        fname = "SummaryStats_HsModelBC_"+ncmap[i]+".txt"
        result=np.round(np.array([mvalstats.smrstat(cmhs_bc,0.1,30.)]),3)
        ifile = open(fname,'w'); ifile.write(hd+' \n')
        np.savetxt(ifile,result.astype('str'),fmt="%s",delimiter='\t') 
        ifile.close(); del ifile, fname, result

        # Table stat metrics
        merr_hs = mvalstats.metrics(cmhs,cohs,vmin=1.,vmax=25.,maxdiff=20., pctlerr='yes')
        dst = pd.DataFrame(np.array([merr_hs]).round(4), columns=merrname)
        dst.to_csv("Merr_Hs_"+ncmap[i]+".txt", index=False, header=True, sep='\t'); del dst, merr_hs

        merr_hs = mvalstats.metrics(cmhs_bc,cohs,vmin=1.,vmax=25.,maxdiff=20., pctlerr='yes')
        dst = pd.DataFrame(np.array([merr_hs]).round(4), columns=merrname)
        dst.to_csv("MerrBC_Hs_"+ncmap[i]+".txt", index=False, header=True, sep='\t'); del dst, merr_hs

        # Plots
        mop=ModelObsPlot(model=[cmhs,cmhs_bc],obs=cohs, mlabels=['PModel','PModel_BC'], axisnames=['PModel','Obs'],vaxisname="Obs and PModel",ftag="Eval_"+ncmap[i]+"_")
        mop.qqplot(); mop.taylordiagram(); mop.scatterplot()
        mop=ModelObsPlot(model=[cmhs,cmhs_bc],obs=cohs, mlabels=['PModel','PModel_BC'], axisnames=['PModel','Obs'],vaxisname="Obs (shaded) and PModel",ftag="Eval_"+ncmap[i]+"_")
        mop.pdf()

        del ind, cohs, cmhs, cmhs_bc, cot, yrs, years
        print(" "); print(" ---- "); print(" ")

