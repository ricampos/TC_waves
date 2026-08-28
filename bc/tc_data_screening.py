#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
 Initial Data Screening
"""

import numpy as np
import pandas as pd
import xarray as xr
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

    # Initial Plot
    mop=ModelObsPlot(model=[mhs,hhs,ghs[0,:],np.nanmean(ghs,axis=0)],obs=ohs, mlabels=['PModel','HAFS','GEFS_c','GEFS_EM'], axisnames=['Model','Obs'],vaxisname="Obs and Models",ftag="Eval_Models_")
    mop.qqplot(); mop.taylordiagram(); mop.scatterplot()

    mop=ModelObsPlot(model=[mhs,hhs,ghs[0,:],np.nanmean(ghs,axis=0)],obs=ohs, mlabels=['PModel','HAFS','GEFS_c','GEFS_EM'], axisnames=['Model','Obs'],vaxisname="Obs (shaded) and Models",ftag="Eval_Models_")
    mop.pdf()

    # Density plots
    plt.close('all')
    scdensityplot(ohs,mhs,1,25.,20.,"Obs","PModel","PModel")
    scdensityplot(ohs,hhs,1,25.,20.,"Obs","HAFS","HAFS")
    scdensityplot(ohs,ghs[0,:],1,25.,20.,"Obs","GEFS_c","GEFSc")
    scdensityplot(ohs,np.nanmean(ghs,axis=0),1,25.,20.,"Obs","GEFS_EM","GEFSEM")

    # Pmodel against HAFS
    scdensityplot(hhs,mhs,1,25.,20.,"HAFS","PModel","PModelHAFS")
    scdensityplot(hhs,ghs[0,:],1,25.,20.,"HAFS","GEFS_c","PModelGEFSc")
    scdensityplot(hhs,np.nanmean(ghs,axis=0),1,25.,20.,"HAFS","GEFS_EM","PModelGEFSEM")

    # Sectors
    for i in np.unique(csec):
        ind=np.where(csec==i)[0]
        scdensityplot(ohs[ind],mhs[ind],1,25.,20.,"Obs","PModel","CSEC_"+repr(i)+"_PModel")
        scdensityplot(ohs[ind],hhs[ind],1,25.,20.,"Obs","HAFS","CSEC_"+repr(i)+"_HAFS")
        scdensityplot(ohs[ind],ghs[0,:][ind],1,25.,20.,"Obs","GEFS_c","CSEC_"+repr(i)+"_GEFSc")
        scdensityplot(ohs[ind],np.nanmean(ghs,axis=0)[ind],1,25.,20.,"Obs","GEFS_EM","CSEC_"+repr(i)+"_GEFSEM")
        del ind

    # RESIDUE
    res=mhs-ohs

    fig1 = plt.figure(1, figsize=(7,6)); ax = fig1.add_subplot(111)
    ax.axhline(y=0, color='black', linestyle='--')
    hb = ax.hexbin(ccdist, res, gridsize=100, cmap='jet', mincnt=1,
                   bins='log', zorder=3)

    ax.locator_params(axis='y', nbins=7)
    ax.locator_params(axis='x', nbins=7)
    ax.set_xlabel("Distance to the Center"); ax.set_ylabel("Residual (Model-Obs)")
    ax.grid(c='grey', ls=':', alpha=0.5, zorder=1)
    plt.tight_layout()
    plt.savefig("Residue_CDIST.png", dpi=200, facecolor='w', edgecolor='w',orientation='portrait', format='png',transparent=False, bbox_inches='tight', pad_inches=0.1)
    plt.close('all')


    fig1 = plt.figure(1, figsize=(7,6)); ax = fig1.add_subplot(111)
    ax.axhline(y=0, color='black', linestyle='--')
    hb = ax.hexbin(ccangle, res, gridsize=100, cmap='jet', mincnt=1,
                   bins='log', zorder=3)

    ax.locator_params(axis='y', nbins=7)
    ax.locator_params(axis='x', nbins=7)
    ax.set_xlabel("R (angle)"); ax.set_ylabel("Residual (Model-Obs)")
    ax.grid(c='grey', ls=':', alpha=0.5, zorder=1)
    plt.tight_layout()
    plt.savefig("Residue_CANGLE.png", dpi=200, facecolor='w', edgecolor='w',orientation='portrait', format='png',transparent=False, bbox_inches='tight', pad_inches=0.1)
    plt.close('all')


    fig1 = plt.figure(1, figsize=(7,6)); ax = fig1.add_subplot(111)
    ax.axhline(y=0, color='black', linestyle='--')
    hb = ax.hexbin(csec, res, gridsize=100, cmap='jet', mincnt=1,
                   bins='log', zorder=3)

    ax.locator_params(axis='y', nbins=7)
    ax.locator_params(axis='x', nbins=7)
    ax.set_xlabel("Sector"); ax.set_ylabel("Residual (Model-Obs)")
    ax.grid(c='grey', ls=':', alpha=0.5, zorder=1)
    plt.tight_layout()
    plt.savefig("Residue_CSEC.png", dpi=200, facecolor='w', edgecolor='w',orientation='portrait', format='png',transparent=False, bbox_inches='tight', pad_inches=0.1)
    plt.close('all')

    # Individual TC
    with open("/home/ricardo/work/noaa/analysis/TC_Waves/1preproc/cyclonemap/cyclone_names.txt", "r") as f:
        cnames = [line.strip() for line in f]

    fcid=[]
    for i in np.unique(cid):
        ind=np.where(cid==i)[0]
        if len(ind)>=10:
            scdensityplot(ohs[ind],mhs[ind],1,25.,20.,"Obs","PModel","TC_"+cnames[i]+"_PModel")
            scdensityplot(ohs[ind],hhs[ind],1,25.,20.,"Obs","HAFS","TC_"+cnames[i]+"_HAFS")
            scdensityplot(ohs[ind],ghs[0,:][ind],1,25.,20.,"Obs","GEFS_c","TC_"+cnames[i]+"_GEFSc")
            scdensityplot(ohs[ind],np.nanmean(ghs,axis=0)[ind],1,25.,20.,"Obs","GEFS_EM","TC_"+cnames[i]+"_GEFSEM")
            fcid=np.append(fcid,int(i))

        del ind

    fcid = np.array(fcid).astype('int') # Total 84 cyclones with valid obs

