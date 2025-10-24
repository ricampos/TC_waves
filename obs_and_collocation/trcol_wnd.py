#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from matplotlib.mlab import *
import pandas as pd
from pylab import *
import os
import time
import timeit
from calendar import timegm
import sys
# from ww3tools https://github.com/NOAA-EMC/WW3-tools
import mvalstats
import pvalstats
from pvalstats import *
# ---------------------------------------------------
import warnings; warnings.filterwarnings("ignore")
# netcdf format
fnetcdf="NETCDF4"

# Font size and style
sl=13
matplotlib.rcParams.update({'font.size': sl}); plt.rc('font', size=sl) 
matplotlib.rc('xtick', labelsize=sl); matplotlib.rc('ytick', labelsize=sl); matplotlib.rcParams.update({'font.size': sl})


# Auxiliar functions for triple collocation

def nregr_B(X,Y,Z,e_x,e_y,e_z):
    """
    Neutral regression (Marsden, 1999), from Janssen et al (2007) 
    https://doi.org/10.1175/JTECH2069.1
    Section 2, equation (6)
    """

    gamma = np.nanmean(e_x**2) / np.nanmean(e_y**2)
    a = gamma * np.nanmean(X*Y)
    b = np.nanmean(X**2) - (gamma * np.nanmean(Y**2))
    c = - np.nanmean(X*Y)
    By = (-b + np.sqrt(b**2 - 4*a*c)) / (2.0 * a)
    del gamma,a,b,c

    gamma = np.nanmean(e_x**2) / np.nanmean(e_z**2)
    a = gamma * np.nanmean(X*Z)
    b = np.nanmean(X**2) - (gamma * np.nanmean(Z**2))
    c = - np.nanmean(X*Z)
    Bz = (-b + np.sqrt(b**2 - 4*a*c)) / (2.0 * a)
    del gamma,a,b,c

    gamma = np.nanmean(e_y**2) / np.nanmean(e_x**2)
    a = gamma * np.nanmean(Y * X)
    b = np.nanmean(Y**2) - (gamma * np.nanmean(X**2))
    c = - np.nanmean(Y * X)
    Bx = (-b + np.sqrt(b**2 - 4*a*c)) / (2.0 * a)
    del gamma,a,b,c

    return Bx,By,Bz


def roerr(X,Y,Z,Bx,By,Bz):
    """
    Residual errors (random observational errors)
    From Janssen et al (2007) Section 2
    https://doi.org/10.1175/JTECH2069.1
    And Houghton et al (2021) Section 2
    https://doi.org/10.1175/JTECH-D-20-0187.1
    """

    Xl = X/Bx
    Yl = Y/By
    Zl = Z/Bz

    e_xl = np.nanmean( (Xl - Yl)*(Xl - Zl) )
    e_yl = np.nanmean( (Yl - Xl)*(Yl - Zl) )
    e_zl = np.nanmean( (Zl - Xl)*(Zl - Yl) )

    e_x = e_xl * Bx
    e_y = e_yl * By
    e_z = e_zl * Bz

    return e_x, e_y, e_z

# visualization
def qqplot(X,T,qcolor,auxmax,xlname,varname,fname):
    '''
    Quantile-Quantile plot.
    '''
    plt.close('all')

    qobs = np.sort(T)
    qm = np.sort(X)

    a=0. ; b=auxmax
    aux=np.linspace(a-0.2*a,b+0.2*a,100)
    # plot
    fig1 = plt.figure(1,figsize=(5,4.5)); ax = fig1.add_subplot(111)
    ax.plot(aux,aux,'k', linewidth=2.,alpha=0.4,zorder=2); del a,b
    ax.plot(qobs,qm, color=qcolor, marker='.', linestyle=' ',linewidth=1.,alpha=0.8,zorder=3)

    plt.grid(c='grey', ls=':', alpha=0.5,zorder=1)
    for i in np.array([50,80,90,95,99]):
        plt.axvline(x= np.nanpercentile(X,int(i)),ls='--',color='grey',linewidth=1.,alpha=0.9,zorder=1)
        plt.text(np.nanpercentile(X,int(i)),(aux.max()-aux.min())/15 + aux.min(),str(int(i))+'th',color='dimgrey',fontsize=sl-7,zorder=4)
        plt.text(np.nanpercentile(X,int(i)),(aux.max()-aux.min())/1.05 + aux.min(),str(int(i))+'th',color='dimgrey',fontsize=sl-7)

    plt.ylim(ymax = auxmax, ymin = 0.7)
    plt.xlim(xmax = auxmax, xmin = 0.7)
    plt.locator_params(axis='y', nbins=7) ; plt.locator_params(axis='x', nbins=7) 

    ax.set_xlabel(xlname); ax.set_ylabel(varname)
    plt.tight_layout()
    plt.savefig(fname, dpi=200, facecolor='w', edgecolor='w',orientation='portrait', format='png',transparent=False, bbox_inches='tight', pad_inches=0.1)
    plt.close(fig1); del fig1, ax


def scatterplot(X,T,qcolor,auxmax,xlname,varname,fname):
    '''
    Scatter plot.
    '''
    plt.close('all')

    qobs = np.array(T)
    qm = np.array(X)

    a=0. ; b=auxmax
    aux=np.linspace(a-0.2*a,b+0.2*a,100)
    # plot
    fig1 = plt.figure(1,figsize=(5,4.5)); ax = fig1.add_subplot(111)
    ax.plot(aux,aux,'k', linewidth=1.,alpha=0.9,zorder=1); del a,b
    ax.plot(aux,aux,'k', linewidth=0.5,alpha=0.4,zorder=3)
    ax.scatter(qobs,qm, color=qcolor, marker='.', linewidth=1.,alpha=0.8,zorder=2)

    plt.grid(c='grey', ls=':', alpha=0.5,zorder=1)
    for i in np.array([50,80,90,95,99]):
        plt.axvline(x= np.nanpercentile(X,int(i)),ls='--',color='grey',linewidth=1.,alpha=0.9,zorder=1)
        plt.text(np.nanpercentile(X,int(i)),(aux.max()-aux.min())/15 + aux.min(),str(int(i))+'th',color='dimgrey',fontsize=sl-7,zorder=4)
        plt.text(np.nanpercentile(X,int(i)),(aux.max()-aux.min())/1.05 + aux.min(),str(int(i))+'th',color='dimgrey',fontsize=sl-7)

    plt.ylim(ymax = auxmax, ymin = 0.7)
    plt.xlim(xmax = auxmax, xmin = 0.7)
    plt.locator_params(axis='y', nbins=7) ; plt.locator_params(axis='x', nbins=7) 

    ax.set_xlabel(xlname); ax.set_ylabel(varname)
    plt.tight_layout()
    plt.savefig(fname+".png", dpi=200, facecolor='w', edgecolor='w',orientation='portrait', format='png',transparent=False, bbox_inches='tight', pad_inches=0.1)
    plt.close(fig1); del fig1, ax

    # Density plot
    fig2 = plt.figure(1,figsize=(5,4.5)); ax = fig2.add_subplot(111)
    ax.plot(aux,aux,'k', linewidth=1.,alpha=0.9,zorder=1)
    ax.plot(aux,aux,'k', linewidth=0.5,alpha=0.4,zorder=3)
    xy = np.vstack([qobs, qm])
    z = gaussian_kde(xy)(xy)
    ax.scatter(qobs, qm, c=z, s=5, cmap=plt.cm.jet, zorder=2)

    plt.grid(c='grey', ls=':', alpha=0.5,zorder=1)
    for i in np.array([50,80,90,95,99]):
        plt.axvline(x= np.nanpercentile(X,int(i)),ls='--',color='grey',linewidth=1.,alpha=0.9,zorder=1)
        plt.text(np.nanpercentile(X,int(i)),(aux.max()-aux.min())/15 + aux.min(),str(int(i))+'th',color='dimgrey',fontsize=sl-7,zorder=4)
        plt.text(np.nanpercentile(X,int(i)),(aux.max()-aux.min())/1.05 + aux.min(),str(int(i))+'th',color='dimgrey',fontsize=sl-7)

    plt.ylim(ymax = auxmax, ymin = 0.7)
    plt.xlim(xmax = auxmax, xmin = 0.7)
    plt.locator_params(axis='y', nbins=7) ; plt.locator_params(axis='x', nbins=7) 

    ax.set_xlabel(xlname); ax.set_ylabel(varname)
    plt.tight_layout()
    plt.savefig(fname+"_denst.png", dpi=200, facecolor='w', edgecolor='w',orientation='portrait', format='png',transparent=False, bbox_inches='tight', pad_inches=0.1)
    plt.close(fig2); del fig2, ax


if __name__ == "__main__":

    # Inputs -----
    # Initial and final date 
    datemin='2022010100'
    datemax='2025010100'

    MTYPE="NDBC"

    # Minimum value for the analysis - to select (or not) more severe events
    vmin = 10.
    # Cyclone only
    cycl='no'
    # -----

    adatemin= np.double(timegm( time.strptime(datemin, '%Y%m%d%H')))
    adatemax= np.double(timegm( time.strptime(datemax, '%Y%m%d%H')))

    # ------- Buoy ---------
    dfx = pd.read_csv("/home/ricardo/work/noaa/analysis/TC_Waves/2collocation/Data_"+MTYPE+".txt", sep='\t')
    btime = np.array( (pd.to_datetime(dfx['time'][:], format='%Y%m%d%H%M') - pd.Timestamp('1970-01-01')) // pd.Timedelta('1s') ).astype('double')
    blat = np.array(dfx['lat'][:]); blon = np.array(dfx['lon'][:])
    cmap = dfx['cmap'].values[:]
    bwnd = dfx['wnd'].values[:]
    bwnd[bwnd<0.]=np.nan
    del dfx
    # -----------

    # -------- Satellite ---------
    satms = np.array(['CFOSAT','HY2B','SARAL','SENTINEL3B','CRYOSAT2','JASON3','SENTINEL3A','SENTINEL6A']).astype('str')
    for i in range(0,len(satms)):
        dfy = pd.read_csv("/home/ricardo/work/noaa/analysis/TC_Waves/2collocation/data_proc/"+MTYPE+"/alt/Data_REF_"+satms[i]+".txt", sep='\t') 
        if i==0:
            swnd = np.zeros((len(dfy)),'f')*np.nan

        aux = np.array(dfy['wnd_avr'].values[:]); aux[aux<0]=np.nan
        swnd = np.nanmean(([swnd],[aux]),axis=0)[0,:]; del aux

        del dfy

    # ----------------                  

    # --------- Model -----------
    dfz = pd.read_csv("/home/ricardo/work/noaa/analysis/TC_Waves/2collocation/data_proc/"+MTYPE+"/gdas/ColData_GDAS.txt", sep='\t') 
    mwnd = dfz['wnd'].values[:]
    mwnd[mwnd<0.]=np.nan
    del dfz
    # ---------------------------

    # Quality Control Altimeter (exclude outliers that contaminate the statistics)
    indqq = np.where( (swnd>mwnd*1.35+3.) | (swnd<mwnd*0.68-3.) )
    if size(indqq)>0.:
        swnd[indqq[0]] = np.nan

    indqq = np.where( (swnd>bwnd*1.35+3.) | (swnd<bwnd*0.68-3.) )
    if size(indqq)>0.:
        swnd[indqq[0]] = np.nan

    # Max diff for Quality Control
    diflim = 7.

    emean=np.nanmean(([bwnd],[swnd],[mwnd]),axis=0)[0,:]
    if cycl=='no':
        ind = np.where( (np.abs(bwnd-emean)<diflim) & (np.abs(swnd-emean)<diflim) & (bwnd>vmin) & (swnd>vmin) & (mwnd>vmin) & (btime>=adatemin) & (btime<=adatemax) )[0]
    if cycl=='yes':
        ind = np.where( (cmap>1) & (np.abs(bwnd-emean)<diflim) & (np.abs(swnd-emean)<diflim) & (bwnd>vmin) & (swnd>vmin) & (mwnd>vmin) & (btime>=adatemin) & (btime<=adatemax) )[0]

    nsz = int(size(ind))
    print(" Number of matchups for the triple collocation: "+repr(nsz))


    # ========== Triple Collocation =======================

    X = bwnd[ind]
    Y = swnd[ind]
    Z = mwnd[ind]

    # first guess of calibration constants
    Bx=1; By=1; Bz=1
    e_x,e_y,e_z = roerr(X,Y,Z,Bx,By,Bz)

    ferr=1.; c=0
    while ferr>1e-6 and c<1000:

        Bx,By,Bz = nregr_B(X,Y,Z,e_x,e_y,e_z)
        # Residual errors estimation (random observational errors)
        ne_x,ne_y,ne_z = roerr(X,Y,Z,Bx,By,Bz)

        # difference per step
        ferr = np.abs(ne_x-e_x)/np.abs(ne_x) + np.abs(ne_y-e_y)/np.abs(ne_y) + np.abs(ne_z-e_z)/np.abs(ne_z)
        # update
        e_x = ne_x; e_y = ne_y; e_z = ne_z

        # print(repr(ferr))
        c=c+1


    # final residual errors (random observational errors)
    rex = np.sqrt(e_x); rey = np.sqrt(e_y); rez = np.sqrt(e_z)

    # Final Approximation of the truth
    Tx = (X/Bx); Ty = (Y/By); Tz = (Z/Bz)
    txyz = (Tx*(1/e_x) + Ty*(1/e_y) + Tz*(1/e_z) ) / ((1/e_x)+(1/e_y)+(1/e_z))

    # Save results
    tpresults = np.array([[nsz, e_x, e_y, e_z, rex, rey, rez, Bx, By, Bz]])
    header = "Size e2_x e2_y e2_z rex rey rez Bx By Bz"
    np.savetxt("tpresults_Wnd_"+str(int(vmin)).zfill(2)+".txt", tpresults, header=header, fmt="%.6f")

    # Save dataset used for the analysis
    df = pd.DataFrame({
        'time': pd.to_datetime(btime[ind], unit='s').strftime('%Y%m%d%H%M'),
        'lat': np.round(blat[ind],5),
        'lon': np.round(blon[ind],5),
        'wnd_buoy': X,
        'wnd_sat': Y,
        'wnd_gdas': Z,
        'wnd_RefTruth': txyz
    })

    df.to_csv("Dataset_TrpCol_"+str(int(vmin)).zfill(2)+".txt", sep='\t', index=False, header=True)


    # Apply traditional validation to compare and assess each estimate.

    # initial data plot
    bdtime = pd.to_datetime(btime[ind], unit='s')
    fig1, ax = plt.subplots(figsize=(11, 5))
    # ax.fill_between(bdtime, 0., txyz, color='silver', alpha=0.5, zorder=1)
    # ax.plot(bdtime, txyz, color='silver', marker='.', linestyle='', linewidth=1.,zorder=1)
    ax.plot(bdtime, X, color="#1f77b4", marker='.', linestyle='', linewidth=2., label='Buoy', zorder=3)
    ax.plot(bdtime, Y, color="darkgreen", marker='.', linestyle='', linewidth=2., label='Altimeter', zorder=3)
    ax.plot(bdtime, Z, color="firebrick", marker='.', linestyle='', linewidth=2., label='Model', zorder=3)
    ax.set_xlabel("Time")
    ax.set_ylabel("WND (m/s)")
    ax.grid()
    ax.set_xlim(bdtime.min(), bdtime.max())
    ax.xaxis.set_major_formatter(DateFormatter('%b%Y'))
    ax.fmt_xdata = DateFormatter('%b%Y')
    ax.legend(fontsize=sl - 3)
    plt.grid(c='grey', ls='--', alpha=0.3, zorder=1)
    plt.savefig("TimeSeries_TrpCol_"+str(int(vmin)).zfill(2)+".png", dpi=200, facecolor='w', edgecolor='w', orientation='portrait',
        format='png', bbox_inches='tight', pad_inches=0.1)
    plt.close(fig1)

    # summary stats
    st_X = mvalstats.smrstat(X)
    st_Y = mvalstats.smrstat(Y)
    st_Z = mvalstats.smrstat(Z)

    ifile = open("Table_SummaryStats_TrpCol_"+str(int(vmin)).zfill(2)+".txt",'w')
    ifile.write("# Summary Stats, total amount of data: "+repr(nsz)+" \n")
    ifile.write('# Buoy, Altimeter, Model \n')
    ifile.write('# mean, variance, skewness, kurtosis, min, max, percentile80, percentile90, percentile95, percentile99, percentile99.9 \n')
    np.savetxt(ifile,np.atleast_2d(st_X) ,fmt="%12.4f",delimiter='	') 
    np.savetxt(ifile,np.atleast_2d(st_Y) ,fmt="%12.4f",delimiter='	') 
    np.savetxt(ifile,np.atleast_2d(st_Z) ,fmt="%12.4f",delimiter='	') 
    ifile.close(); del ifile

    # error
    merr_X = np.array(mvalstats.metrics(X,txyz, pctlerr='yes'))
    merr_Y = np.array(mvalstats.metrics(Y,txyz, pctlerr='yes'))
    merr_Z = np.array(mvalstats.metrics(Z,txyz, pctlerr='yes'))
    ifile = open("Table_ErrMetrics_TrpCol_"+str(int(vmin)).zfill(2)+".txt",'w')
    ifile.write("# Error metrics, total amount of data: "+repr(nsz)+" \n")
    ifile.write('# Buoy, Altimeter, Model \n')
    ifile.write('# bias, RMSE, NBias, NRMSE, SCrmse, SI, HH, CC, N, Bias95p, RMSE95p, N95p \n')
    np.savetxt(ifile,np.atleast_2d(merr_X) ,fmt="%12.4f",delimiter='	') 
    np.savetxt(ifile,np.atleast_2d(merr_Y) ,fmt="%12.4f",delimiter='	') 
    np.savetxt(ifile,np.atleast_2d(merr_Z) ,fmt="%12.4f",delimiter='	') 
    ifile.close(); del ifile

    # QQ-plots, Scatter Plots, Taylor Diagram

    auxmax = 30.
    qqplot(X,txyz,"#1f77b4",auxmax,"Reference truth (T)","Buoy","QQplot_X_TrpCol_"+str(int(vmin)).zfill(2)+".png")
    qqplot(Y,txyz,"darkgreen",auxmax,"Reference truth (T)","Altimeter","QQplot_Y_TrpCol_"+str(int(vmin)).zfill(2)+".png")
    qqplot(Z,txyz,"firebrick",auxmax,"Reference truth (T)","Model","QQplot_Z_TrpCol_"+str(int(vmin)).zfill(2)+".png")
    qqplot(Y,X,"k",auxmax,"Buoy","Altimeter","QQplot_XY_TrpCol_"+str(int(vmin)).zfill(2)+".png")

    scatterplot(X,txyz,"#1f77b4",auxmax,"Reference truth (T)","Buoy","Scatterplot_X_TrpCol_"+str(int(vmin)).zfill(2))
    scatterplot(Y,txyz,"darkgreen",auxmax,"Reference truth (T)","Altimeter","Scatterplot_Y_TrpCol_"+str(int(vmin)).zfill(2))
    scatterplot(Z,txyz,"firebrick",auxmax,"Reference truth (T)","Model","Scatterplot_Z_TrpCol_"+str(int(vmin)).zfill(2))
    scatterplot(Y,X,"dimgrey",auxmax,"Buoy","Altimeter","Scatterplot_XY_TrpCol_"+str(int(vmin)).zfill(2))

    mop=ModelObsPlot(model=[X,Z,Y],obs=txyz, mlabels=['Buoy','Model','Altimeter'],ftag="TrpCol_"+str(int(vmin)).zfill(2)+"_")
    mop.taylordiagram()

    print("🌊 Triple Collocation completed :) 🌊")

