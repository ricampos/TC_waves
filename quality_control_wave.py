#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
quality_control_wave.py

VERSION AND LAST UPDATE:
 v1.0  03/19/2025

PURPOSE:
 Quality control of wave observations. Checks include:
 - Data Range
 - Redundancy and Duplicates
 - Rate-of-Change Tests
 - Exclude records based on water depth and distance to the coast
 - Use model data (GDAS) to evaluate and exclude discrepant observations

 https://github.com/ioos/ioos_qc/blob/master/ioos_qc/qartod.py#L46-L77

USAGE:
 functions
   data_range
   duplicates
   rate_of_change
   landcoast_exclude
   model_compare
   wsra (dedicated to P3 airbone radar WSRA)

 They read a dictionary, from wread.py containing lat,lon,time, and the wave parameters.

OUTPUT:
 Apart from wsra that returns the dictionary, all the other functions return a numpy array of Hs.

DEPENDENCIES:
 See setup.py and the imports below.

AUTHOR and DATE:
 03/19/2025: Ricardo M. Campos, first version.

PERSON OF CONTACT:
 Ricardo M Campos: ricardo.campos@noaa.gov

"""

from pylab import *
import pandas as pd
import xarray as xr
import netCDF4 as nc
import numpy as np
import yaml
from haversine import haversine, Unit
import warnings; warnings.filterwarnings("ignore")


# Data Range Checks
def data_range(wdic,var='hs',vmin=None,vmax=None):
    """

    """

    if vmin==None:
        if var=='hs' or var=='Hs' or var=='swh':
            vmin=0.1
        elif var=='tp' or var=='Tp':
            vmin=1.
        elif var=='u10' or var=='U10' or var=='WSP' or var=='wsp':
            vmin=0.
        elif var=='dir' or var=='dp' or var=='Dp' or var=='dm' or var=='Dm':
            vmin=-180.;

    if vmax==None:
        if var=='hs' or var=='Hs' or var=='swh':
            vmax=20.
        elif var=='tp' or var=='Tp':
            vmax=30.
        elif var=='u10' or var=='U10' or var=='WSP' or var=='wsp':
            vmax=80.
        elif var=='dir' or var=='dp' or var=='Dp' or var=='dm' or var=='Dm':
            vmax=360.;

    fmvar = np.array(wdic[var])
    ind = np.where((fmvar<vmin)|(fmvar>vmax))
    if np.size(ind)>0:
        for i in range(0,len(ind[0])):
            print(" QC data_range - Cleaned data: index "+repr(ind[0][i])+" Hs "+repr(fmvar[ind[0][i]]))
            fmvar[ind[0][i]]=np.nan

    return fmvar
    del fmvar


# Redundancy and Duplicates
def duplicates(wdic):
    """

    """

    fhs = np.array(wdic['hs'][:])
    for i in range(0,len(wdic['time'])):
        qf=0

        if np.isnan(wdic['hs'][i])==False:

            c=0; stp=0; j=np.nan
            while np.isnan(j)==True and stp==0:
                c=c+1
                if i-c>=0:
                    if np.isnan(wdic['hs'][i-c])==False:
                        j=c; stp=1
                else:
                    stp=1

            if i-j>=0:
                if wdic['hs'][i]==wdic['hs'][i-j]:
                    qf=1

            if qf>0:
                print(" QC duplicates - Cleaned data: index "+repr(i)+" Hs "+repr(wdic['hs'][i]))
                fhs[i]=np.nan

        del qf

    return fhs
    del fhs


# Rate-of-Change Tests
def rate_of_change(wdic):
    """

    """

    fhs = np.array(wdic['hs'][:])
    cfhs=fhs*np.nan; fqf=fhs*np.nan

    for i in range(0,len(wdic['time'])):
        qf=0

        if np.isnan(wdic['hs'][i])==False:

            c=0; stp=0; j=np.nan
            while np.isnan(j)==True and stp==0:
                c=c+1
                if i-c>=0:
                    if np.isnan(wdic['hs'][i-c])==False:
                        j=c; stp=1
                else:
                    stp=1

            c=0; stp=0; k=np.nan
            while np.isnan(k)==True and stp==0:
                c=c+1
                if i+c<len(wdic['time']):
                    if np.isnan(wdic['hs'][i+c])==False:
                        k=c; stp=1
                else:
                    stp=1

            if i-j>=0:

                dist1 = float(haversine((wdic['latitude'][i], wdic['longitude'][i]), (wdic['latitude'][i-j], wdic['longitude'][i-j]), unit=Unit.KILOMETERS))
                dtime1 = float(wdic['time'][i]-wdic['time'][i-j])/3600

                dmhs1=wdic['hs'][i]-wdic['hs'][i-j]

                if dtime1<=1.1 and dist1<=15:
                    if dmhs1>4:
                        qf=1

                elif dtime1<=3.1 and dist1<=30:
                    if dmhs1>5.5:
                        qf=1

                elif dtime1<=6.1 and dist1<=60:
                    if dmhs1>7:
                        qf=1

            # Spiking Test
            if qf==0 and i+k<len(wdic['time']):

                dist2 = float(haversine((wdic['latitude'][i], wdic['longitude'][i]), (wdic['latitude'][i+k], wdic['longitude'][i+k]), unit=Unit.KILOMETERS))
                dtime2 = float(wdic['time'][i+k]-wdic['time'][i])/3600

                dmhs2=wdic['hs'][i]-wdic['hs'][i+k]

                if dtime2<=1.1 and dist2<=15:
                    if dmhs2>4:
                        qf=2

                elif dtime2<=3.1 and dist2<=30:
                    if dmhs2>5.5:
                        qf=2

                elif dtime2<=6.1 and dist2<=60:
                    if dmhs2>7:
                        qf=2

            if qf==0 and i-j>=0 and i+k<len(wdic['time']):

                # up and down
                if dtime2<=1.1 and dist2<=15:
                    if dmhs1+dmhs2> 4+int(wdic['hs'][i]/7):
                        qf=3

                elif dtime2<=3.1 and dist2<=30:
                    if dmhs1+dmhs2> 5+int(wdic['hs'][i]/7):
                        qf=3

                elif dtime2<=6.1 and dist2<=60:
                    if dmhs1+dmhs2> 7+int(wdic['hs'][i]/7):
                        qf=3

            if qf>0:
                print(" QC rate_of_change - Cleaned data: index "+repr(i)+" Hs "+repr(wdic['hs'][i]))
                cfhs[i]=wdic['hs'][i]
                fhs[i]=np.nan
                fqf[i]=qf

        del qf

    # plot(wdic['hs'],'bs'); plot(fhs,'gx'); plot(cfhs,'r.')
    # plot(wdic['date'],wdic['hs'],'bs'); plot(wdic['date'],fhs,'gx'); plot(wdic['date'],cfhs,'r.'); plot(wdic['date'],fqf,'y+')

    return fhs
    del fhs, cfhs, fqf


# Exclude records based on water depth and distance to the coast 
def landcoast_exclude(wdic,gpath=None,mdepth=None,mdfc=None,mdist=None):
    """

    """

    if gpath==None:
       gpath='/home/ricardo/work/noaa/analysis/TC_Waves/2collocation/gridInfo.nc'

    if mdepth==None:
        mdepth=20.

    if mdfc==None:
        mdfc=5.

    if mdist==None:
        mdist=50.

    f=nc.Dataset(gpath)
    glat=f.variables['latitude'][:]; glon=f.variables['longitude'][:]; glon[glon>180]=glon[glon>180]-360.
    dfc=f.variables['distcoast'][:,:]; depth=f.variables['depth'][:,:]
    f.close(); del f

    fhs = np.array(wdic['hs'][:])
    for i in range(0,len(wdic['time'])):
        if np.isnan(wdic['hs'][i])==False:
            indlat=np.min(np.where( np.abs(wdic['latitude'][i]-glat)==np.nanmin(np.abs(wdic['latitude'][i]-glat))))
            indlon=np.min(np.where( np.abs(wdic['longitude'][i]-glon)==np.nanmin(np.abs(wdic['longitude'][i]-glon))))

            dist = float(haversine((wdic['latitude'][i], wdic['longitude'][i]), (glat[indlat], glon[indlon]), unit=Unit.KILOMETERS))
            if dist<mdist:
                if depth[indlat,indlon]<mdepth or dfc[indlat,indlon]<mdfc:
                    print(" QC landcoast_exclude - Cleaned data: index "+repr(i)+" Hs "+repr(wdic['hs'][i])+", wdepth "+repr(depth[indlat,indlon])+"  dist_to_coast "+repr(dfc[indlat,indlon])) 
                    fhs[i]=np.nan

    return fhs
    del fhs


# Use model data (GDAS) to evaluate and exclude observations
def model_compare(wdic,gpath=None,mdist=None):
    """

    """

    if gpath==None:
       gpath="/home/ricardo/work/noaa/analysis/TC_Waves/data/GDAS"

    if mdist==None:
        mdist=50.

    # Read GDAS
    ds = xr.open_mfdataset(gpath+"/gdaswave.*.global.0p16.nc", combine='by_coords')
    gtime = np.array( (pd.to_datetime(ds['time'][:], format='%Y%m%d%H') - pd.Timestamp('1970-01-01')) // pd.Timedelta('1s') ).astype('double')
    glat = np.array(ds['latitude'][:]).astype('float'); glon = np.array(ds['longitude'][:]).astype('float'); glon[glon>180]=glon[glon>180]-360.
    ghs = ds['HTSGW_surface']
    # ---------------

    fhs = np.array(wdic['hs'][:])
    for i in range(0,len(wdic['time'])):
        if np.isnan(wdic['hs'][i])==False:
            indt = np.where( abs(gtime[:]-wdic['time'][i]) < 3600. )
            if np.size(indt)>0:
                indlat=np.min(np.where( np.abs(wdic['latitude'][i]-glat)==np.nanmin(np.abs(wdic['latitude'][i]-glat))))
                indlon=np.min(np.where( np.abs(wdic['longitude'][i]-glon)==np.nanmin(np.abs(wdic['longitude'][i]-glon))))
                dist = float(haversine((wdic['latitude'][i], wdic['longitude'][i]), (glat[indlat], glon[indlon]), unit=Unit.KILOMETERS))
                if dist<mdist:
                    fghs=float(ghs[np.nanmin(indt),indlat,indlon])
                    if wdic['hs'][i]>fghs*1.3+1. or wdic['hs'][i]<fghs*0.75-0.8:
                        print(" QC model_compare - Cleaned data: index "+repr(i)+" Hs "+repr(wdic['hs'][i])+", model value "+repr(fghs)) 
                        fhs[i]=np.nan

                    del fghs

    return fhs
    del fhs


# P-3 airbone radar
def wsra(wdic):
    """

    """

    # Range
    ind = np.where( (wdic['wind_spd']<0) | (wdic['wind_spd']>80) |
        (wdic['wind_dir']<-180) | (wdic['wind_dir']>360) |
        (wdic['hs']<0.3) | (wdic['hs']>20) |
        (wdic['dhs']<0.3) | (wdic['dhs']>20) |
        (wdic['dp']<-180) | (wdic['dp']>360) |
        (np.nanmean(wdic['rainfall_rate'],axis=1)<0) | (np.nanmean(wdic['rainfall_rate'],axis=1)>200) |
        (wdic['rainfall_rate_median']<0) | (wdic['rainfall_rate_median']>50) |
        (wdic['pralt']<1000) | (wdic['pralt']>4000) |
        (wdic['pseed']<80) | (wdic['pseed']>250) |
        (wdic['wcroll'][:,0]<-2.5) | (wdic['wcroll'][:,0]>2.5) |
        (wdic['wcroll'][:,1]<-2.5) | (wdic['wcroll'][:,1]>2.5) |
        (wdic['wcroll'][:,2]<-2.5) | (wdic['wcroll'][:,2]>2.5) |
        (wdic['wcroll'][:,3]<-2.5) | (wdic['wcroll'][:,3]>2.5) |
        (wdic['wcroll'][:,4]<-2.5) | (wdic['wcroll'][:,4]>2.5) )

    if np.size(ind)>0:
        for i in range(0,len(ind[0])):
            print(" QC WSRA data_range - Cleaned data: index "+repr(ind[0][i])+" Hs "+repr(wdic['hs'][ind[0][i]]))
            wdic['hs'][ind[0][i]]=np.nan
            wdic['dhs'][ind[0][i]]=np.nan
            wdic['dp'][ind[0][i]]=np.nan
            wdic['wind_spd'][ind[0][i]]=np.nan
            wdic['wind_dir'][ind[0][i]]=np.nan

    return wdic
    del wdic

