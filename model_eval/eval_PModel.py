#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

Evaluation of parametric wave model PModel

"""

import numpy as np
# from matplotlib.mlab import *
import pandas as pd
# from pylab import *
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
# netcdf format
fnetcdf="NETCDF4"


import matplotlib.colors as colors
sl=13 # plot style configuration
matplotlib.rcParams.update({'font.size': sl}); plt.rc('font', size=sl) 
matplotlib.rc('xtick', labelsize=sl); matplotlib.rc('ytick', labelsize=sl); matplotlib.rcParams.update({'font.size': sl})


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

    idx = np.where( (np.char.startswith(oid, 'DWSD')!=True) | (np.char.startswith(oid, 'MICROSWIFT')!=True) | (np.char.startswith(oid, 'SAILDRONE')!=True))
    oid=oid[idx]; ohs=ohs[idx]; mhs=mhs[idx]; cmap=cmap[idx]; ot=ot[idx]

    # Headers
    hd = 'mean, variance, skewness, kurtosis, min, max, percentile80, percentile90, percentile95, percentile99, percentile99.9'
    merrname=["bias","RMSE","NBias","NRMSE","SCrmse","SI","HH","CC","N","bias_p95","RMSE_p95","N_p95"]

    acmap=np.array([2,3,4,5]); ncmap=np.array(["Disturbance","ExtrCycl","Subtrop","TropCycl"]).astype('str')

    for i in range(0,len(acmap)):

        # ---- Hs ----
        ind = np.where( (ohs>1.) & (mhs>1.) & (cmap==acmap[i]) )
        if np.size(ind)>0:
            ind=ind[0]

            # Initial summary stats
            fname = "SummaryStats_HsObs_"+ncmap[i]+".txt"
            result=np.round(np.array([mvalstats.smrstat(ohs,0.1,30.)]),3)
            ifile = open(fname,'w'); ifile.write(hd+' \n')
            np.savetxt(ifile,result.astype('str'),fmt="%s",delimiter='\t') 
            ifile.close(); del ifile, fname, result

            # Table stat metrics
            merr_hs = mvalstats.metrics(mhs[ind],ohs[ind],vmin=1.,vmax=25.,maxdiff=20., pctlerr='yes')
            dst = pd.DataFrame(np.array([merr_hs]).round(4), columns=merrname)
            dst.to_csv("Merr_Hs_"+ncmap[i]+".txt", index=False, header=True, sep='\t'); del dst, merr_hs

            # Plots
            mop=ModelObsPlot(model=mhs[ind],obs=ohs[ind], axisnames=['PModel','Obs'],vaxisname="Obs and PModel",ftag="Eval_"+ncmap[i]+"_")
            mop.qqplot(); mop.taylordiagram(); mop.scatterplot()
            mop=ModelObsPlot(model=mhs[ind],obs=ohs[ind], axisnames=['PModel','Obs'],vaxisname="Obs (shaded) and PModel",ftag="Eval_"+ncmap[i]+"_")
            mop.pdf()

            # time-series
            # fig1, ax = plt.subplots(figsize=(9, 4))
            # ax.plot(pd.to_datetime(ot[ind], format='%Y%m%d%H%M'), ohs[ind], color='dimgray', marker='.', linestyle='', linewidth=2., label='Obs', zorder=2)
            # ax.plot(pd.to_datetime(ot[ind], format='%Y%m%d%H%M'), mhs[ind], color='firebrick', marker='.', linestyle='', linewidth=2., label='PModel', zorder=2)
            ## ax.set_xlim(self.dtime[0], self.dtime[-1])
            ## ax.xaxis.set_major_formatter(DateFormatter('%b%d'))
            ## ax.fmt_xdata = DateFormatter('%b%d')
            # ax.set_xlabel("Time"); ax.set_ylabel("Hs (m)")
            # plt.grid(c='grey', ls='--', alpha=0.3, zorder=1)
            # plt.savefig('TimeSeries.png', dpi=200, facecolor='w', edgecolor='w', orientation='portrait',
            #         format='png', bbox_inches='tight', pad_inches=0.1)
            # plt.close(fig1)

            del ind


