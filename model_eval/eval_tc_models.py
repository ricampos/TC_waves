#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

Evaluation of TC models

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

    # Read Obs and PMODEL
    df = pd.read_csv('Data_Obs_PModel_Default.txt', sep='\t')
    df = df.replace(-999., np.nan)
    ot = df.time.values; oid = np.array(df['id']).astype('str'); cmap = np.array(df['cmap'])
    ohs = np.array(df['obs_hs']); otp = np.array(df['obs_tp']); ownd = np.array(df['obs_wnd'])
    # PModel
    pmhs = np.array(df['pm_hs']); pmtp = np.array(df['pm_tp'])
    # HAFS
    dh = pd.read_csv('Data_HAFS_fromObs.txt', sep='\t')
    dh = dh.replace(-999., np.nan)
    hahs = np.array(dh['HTSGW_surface']); hatp = np.array(dh['PERPW_surface'])
    hawnd = np.sqrt( np.array(dh['UGRD_surface'])**2 + np.array(dh['VGRD_surface'])**2 )
    # GEFS
    dg = pd.read_csv('Data_GEFS_fromObs_m00.txt', sep='\t')
    dg = dg.replace(-999., np.nan)
    dg = pd.read_csv('Data_GEFS_fromObs_m00.txt', sep='\t')
    ghs = np.array(dg['HTSGW_surface']); gtp = np.array(dg['PERPW_surface'])
    gwnd = np.sqrt( np.array(dg['UGRD_surface'])**2 + np.array(dg['VGRD_surface'])**2 )

    # -------------------------------------
    # Headers
    hd = 'mean, variance, skewness, kurtosis, min, max, percentile80, percentile90, percentile95, percentile99, percentile99.9'
    merrname=["bias","RMSE","NBias","NRMSE","SCrmse","SI","HH","CC","N","bias_p95","RMSE_p95","N_p95"]

    acmap=np.array([2,3,4,5]); ncmap=np.array(["Disturbance","ExtrCycl","Subtrop","TropCycl"]).astype('str')

    # ---- Hs ----
    for i in range(0,len(acmap)):
        ind = np.where( (ohs>0.1) & (pmhs>0.1) & (hahs>0.1) & (ghs>0.1) & (cmap==acmap[i]) )
        if np.size(ind)>0:
            ind=ind[0]
            print( ncmap[i]+" "+repr(len(ind)) )
            # ----- Initial summary stats ------
            fname = "SummaryStats_Hs_Obs_"+ncmap[i]+".txt"
            result=np.round(np.array([mvalstats.smrstat(ohs,0.1,30.)]),3)
            ifile = open(fname,'w'); ifile.write(hd+' \n')
            np.savetxt(ifile,result.astype('str'),fmt="%s",delimiter='\t') 
            ifile.close(); del ifile, fname, result
            # 
            fname = "SummaryStats_Hs_PMODEL_"+ncmap[i]+".txt"
            result=np.round(np.array([mvalstats.smrstat(pmhs,0.1,30.)]),3)
            ifile = open(fname,'w'); ifile.write(hd+' \n')
            np.savetxt(ifile,result.astype('str'),fmt="%s",delimiter='\t') 
            ifile.close(); del ifile, fname, result
            # 
            fname = "SummaryStats_Hs_HAFS_"+ncmap[i]+".txt"
            result=np.round(np.array([mvalstats.smrstat(hahs,0.1,30.)]),3)
            ifile = open(fname,'w'); ifile.write(hd+' \n')
            np.savetxt(ifile,result.astype('str'),fmt="%s",delimiter='\t') 
            ifile.close(); del ifile, fname, result
            # 
            fname = "SummaryStats_Hs_GEFS_"+ncmap[i]+".txt"
            result=np.round(np.array([mvalstats.smrstat(ghs,0.1,30.)]),3)
            ifile = open(fname,'w'); ifile.write(hd+' \n')
            np.savetxt(ifile,result.astype('str'),fmt="%s",delimiter='\t') 
            ifile.close(); del ifile, fname, result

            # ---- Table stat metrics ----
            merr_hs = mvalstats.metrics(pmhs[ind],ohs[ind],vmin=1.,vmax=25.,maxdiff=20., pctlerr='yes')
            dst = pd.DataFrame(np.array([merr_hs]).round(4), columns=merrname)
            dst.to_csv("Merr_Hs_PMODEL_"+ncmap[i]+".txt", index=False, header=True, sep='\t'); del dst, merr_hs
            #
            merr_hs = mvalstats.metrics(hahs[ind],ohs[ind],vmin=1.,vmax=25.,maxdiff=20., pctlerr='yes')
            dst = pd.DataFrame(np.array([merr_hs]).round(4), columns=merrname)
            dst.to_csv("Merr_Hs_HAFS_"+ncmap[i]+".txt", index=False, header=True, sep='\t'); del dst, merr_hs
            #
            merr_hs = mvalstats.metrics(ghs[ind],ohs[ind],vmin=1.,vmax=25.,maxdiff=20., pctlerr='yes')
            dst = pd.DataFrame(np.array([merr_hs]).round(4), columns=merrname)
            dst.to_csv("Merr_Hs_GEFS_"+ncmap[i]+".txt", index=False, header=True, sep='\t'); del dst, merr_hs

            # ---- Plots ----
            mop=ModelObsPlot(model=[pmhs[ind],hahs[ind],ghs[ind]],obs=ohs[ind], mlabels=['PModel','HAFS','GEFS'], axisnames=['Model','Obs'],vaxisname="Obs and Model",ftag="Eval_"+ncmap[i]+"_")
            mop.qqplot(); mop.taylordiagram(); mop.scatterplot()
            mop=ModelObsPlot(model=[pmhs[ind],hahs[ind],ghs[ind]],obs=ohs[ind], mlabels=['PModel','HAFS','GEFS'], axisnames=['Model','Obs'],vaxisname="Hs (m)",ftag="Eval_"+ncmap[i]+"_")
            mop.pdf()

            # time-series
            fig1, ax = plt.subplots(figsize=(9, 4))
            ax.plot(ohs[ind], color='dimgray', marker='.', linestyle='', linewidth=2., zorder=2)
            ax.plot(ohs[ind], color='k', linewidth=2., label='Obs', zorder=2)
            ax.plot(pmhs[ind], color='blue', linewidth=2., label='PModel', zorder=2)
            ax.plot(hahs[ind], color='red', linewidth=2., label='HAFS', zorder=2)
            ax.plot(ghs[ind], color='green', linewidth=2., label='GEFS', zorder=2)
            ax.legend()
            ax.set_xlabel("Time"); ax.set_ylabel("Hs (m)")
            plt.grid(c='grey', ls='--', alpha=0.3, zorder=1)
            plt.savefig('TimeSeries_Hs_1.png', dpi=200, facecolor='w', edgecolor='w', orientation='portrait',format='png', bbox_inches='tight', pad_inches=0.1)
            # plt.close(fig1)

            del ind


    # ---- Tp ----
    for i in range(0,len(acmap)):
        ind = np.where( (otp>0.1) & (pmtp>0.1) & (hatp>0.1) & (gtp>0.1) & (cmap==acmap[i]) )
        if np.size(ind)>0:
            ind=ind[0]
            print( ncmap[i]+" "+repr(len(ind)) )
            # ----- Initial summary stats ------
            fname = "SummaryStats_Tp_Obs_"+ncmap[i]+".txt"
            result=np.round(np.array([mvalstats.smrstat(otp,0.1,30.)]),3)
            ifile = open(fname,'w'); ifile.write(hd+' \n')
            np.savetxt(ifile,result.astype('str'),fmt="%s",delimiter='\t') 
            ifile.close(); del ifile, fname, result
            # 
            fname = "SummaryStats_Tp_PMODEL_"+ncmap[i]+".txt"
            result=np.round(np.array([mvalstats.smrstat(pmtp,0.1,30.)]),3)
            ifile = open(fname,'w'); ifile.write(hd+' \n')
            np.savetxt(ifile,result.astype('str'),fmt="%s",delimiter='\t') 
            ifile.close(); del ifile, fname, result
            # 
            fname = "SummaryStats_Tp_HAFS_"+ncmap[i]+".txt"
            result=np.round(np.array([mvalstats.smrstat(hatp,0.1,30.)]),3)
            ifile = open(fname,'w'); ifile.write(hd+' \n')
            np.savetxt(ifile,result.astype('str'),fmt="%s",delimiter='\t') 
            ifile.close(); del ifile, fname, result
            # 
            fname = "SummaryStats_Tp_GEFS_"+ncmap[i]+".txt"
            result=np.round(np.array([mvalstats.smrstat(gtp,0.1,30.)]),3)
            ifile = open(fname,'w'); ifile.write(hd+' \n')
            np.savetxt(ifile,result.astype('str'),fmt="%s",delimiter='\t') 
            ifile.close(); del ifile, fname, result

            # ---- Table stat metrics ----
            merr_tp = mvalstats.metrics(pmtp[ind],otp[ind],vmin=1.,vmax=25.,maxdiff=20., pctlerr='yes')
            dst = pd.DataFrame(np.array([merr_tp]).round(4), columns=merrname)
            dst.to_csv("Merr_Tp_PMODEL_"+ncmap[i]+".txt", index=False, header=True, sep='\t'); del dst, merr_tp
            #
            merr_tp = mvalstats.metrics(hatp[ind],otp[ind],vmin=1.,vmax=25.,maxdiff=20., pctlerr='yes')
            dst = pd.DataFrame(np.array([merr_tp]).round(4), columns=merrname)
            dst.to_csv("Merr_Tp_HAFS_"+ncmap[i]+".txt", index=False, header=True, sep='\t'); del dst, merr_tp
            #
            merr_tp = mvalstats.metrics(gtp[ind],otp[ind],vmin=1.,vmax=25.,maxdiff=20., pctlerr='yes')
            dst = pd.DataFrame(np.array([merr_tp]).round(4), columns=merrname)
            dst.to_csv("Merr_Tp_GEFS_"+ncmap[i]+".txt", index=False, header=True, sep='\t'); del dst, merr_tp

            # ---- Plots ----
            mop=ModelObsPlot(model=[pmtp[ind],hatp[ind],gtp[ind]],obs=otp[ind], mlabels=['PModel','HAFS','GEFS'], axisnames=['Model','Obs'],vaxisname="Obs and Model",ftag="Eval_"+ncmap[i]+"_")
            mop.qqplot(); mop.taylordiagram(); mop.scatterplot()
            mop=ModelObsPlot(model=[pmtp[ind],hatp[ind],gtp[ind]],obs=otp[ind], mlabels=['PModel','HAFS','GEFS'], axisnames=['Model','Obs'],vaxisname="Tp (s)",ftag="Eval_"+ncmap[i]+"_")
            mop.pdf()

            # time-series
            fig1, ax = plt.subplots(figsize=(9, 4))
            ax.plot(otp[ind], color='dimgray', marker='.', linestyle='', linewidth=2., zorder=2)
            ax.plot(otp[ind], color='k', linewidth=2., label='Obs', zorder=2)
            ax.plot(pmtp[ind], color='blue', linewidth=2., label='PModel', zorder=2)
            ax.plot(hatp[ind], color='red', linewidth=2., label='HAFS', zorder=2)
            ax.plot(gtp[ind], color='green', linewidth=2., label='GEFS', zorder=2)
            ax.legend()
            ax.set_xlabel("Time"); ax.set_ylabel("Tp (s)")
            plt.grid(c='grey', ls='--', alpha=0.3, zorder=1)
            plt.savefig('TimeSeries_Tp_1.png', dpi=200, facecolor='w', edgecolor='w', orientation='portrait',format='png', bbox_inches='tight', pad_inches=0.1)
            # plt.close(fig1)

            del ind


    # ---- WND ----
    hawnd[hawnd>150]=np.nan; gwnd[gwnd>150]=np.nan; ownd[ownd>150]=np.nan
    for i in range(0,len(acmap)):
        ind = np.where( (ownd>0.1) & (hawnd>0.1) & (gwnd>0.1) & (cmap==acmap[i]) )
        if np.size(ind)>0:
            ind=ind[0]
            print( ncmap[i]+" "+repr(len(ind)) )
            # ----- Initial summary stats ------
            fname = "SummaryStats_WND_Obs_"+ncmap[i]+".txt"
            result=np.round(np.array([mvalstats.smrstat(ownd,0.1,30.)]),3)
            ifile = open(fname,'w'); ifile.write(hd+' \n')
            np.savetxt(ifile,result.astype('str'),fmt="%s",delimiter='\t') 
            ifile.close(); del ifile, fname, result
            # 
            fname = "SummaryStats_WND_HAFS_"+ncmap[i]+".txt"
            result=np.round(np.array([mvalstats.smrstat(hawnd,0.1,30.)]),3)
            ifile = open(fname,'w'); ifile.write(hd+' \n')
            np.savetxt(ifile,result.astype('str'),fmt="%s",delimiter='\t') 
            ifile.close(); del ifile, fname, result
            # 
            fname = "SummaryStats_WND_GEFS_"+ncmap[i]+".txt"
            result=np.round(np.array([mvalstats.smrstat(gwnd,0.1,30.)]),3)
            ifile = open(fname,'w'); ifile.write(hd+' \n')
            np.savetxt(ifile,result.astype('str'),fmt="%s",delimiter='\t') 
            ifile.close(); del ifile, fname, result

            # ---- Table stat metrics ----
            merr_wnd = mvalstats.metrics(hawnd[ind],ownd[ind],vmin=1.,vmax=25.,maxdiff=20., pctlerr='yes')
            dst = pd.DataFrame(np.array([merr_wnd]).round(4), columns=merrname)
            dst.to_csv("Merr_WND_HAFS_"+ncmap[i]+".txt", index=False, header=True, sep='\t'); del dst, merr_wnd
            #
            merr_wnd = mvalstats.metrics(gwnd[ind],ownd[ind],vmin=1.,vmax=25.,maxdiff=20., pctlerr='yes')
            dst = pd.DataFrame(np.array([merr_wnd]).round(4), columns=merrname)
            dst.to_csv("Merr_WND_GEFS_"+ncmap[i]+".txt", index=False, header=True, sep='\t'); del dst, merr_wnd

            # ---- Plots ----
            mop=ModelObsPlot(model=[hawnd[ind],gwnd[ind]],obs=ownd[ind], mlabels=['HAFS','GEFS'], axisnames=['Model','Obs'],vaxisname="Obs and Model",ftag="Eval_"+ncmap[i]+"_")
            mop.qqplot(); mop.taylordiagram(); mop.scatterplot()
            mop=ModelObsPlot(model=[hawnd[ind],gwnd[ind]],obs=ownd[ind], mlabels=['HAFS','GEFS'], axisnames=['Model','Obs'],vaxisname="U10 (m/s)",ftag="Eval_"+ncmap[i]+"_")
            mop.pdf()

            # time-series
            fig1, ax = plt.subplots(figsize=(9, 4))
            ax.plot(ownd[ind], color='dimgray', marker='.', linestyle='', linewidth=2., zorder=2)
            ax.plot(ownd[ind], color='k', linewidth=2., label='Obs', zorder=2)
            ax.plot(hawnd[ind], color='blue', linewidth=2., label='HAFS', zorder=2)
            ax.plot(gwnd[ind], color='red', linewidth=2., label='GEFS', zorder=2)
            ax.legend()
            ax.set_xlabel("Time"); ax.set_ylabel("U10 (m/s)")
            plt.grid(c='grey', ls='--', alpha=0.3, zorder=1)
            plt.savefig('TimeSeries_WND_1.png', dpi=200, facecolor='w', edgecolor='w', orientation='portrait',format='png', bbox_inches='tight', pad_inches=0.1)
            # plt.close(fig1)

            del ind


