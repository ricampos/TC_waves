#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
 Extract reanalysis data to pair with observations
"""

import numpy as np
import pandas as pd
import netCDF4 as nc

if __name__ == "__main__":

    pmrun="Default"

    # Read Obs
    df = pd.read_csv('Data_Cyclones.txt', sep='\t')
    ot = np.array(df['time'].values[:]).astype('str'); ot = pd.to_datetime(ot, format='%Y%m%d%H%M'); ot = np.array(ot.view('int64') // 1_000_000_000).astype('double')
    flat = np.array(df['lat'].values[:]); flon = np.array(df['lon'].values[:])
    gidlat = np.array(df['gidlat'].values[:]); gidlon = np.array(df['gidlon'].values[:])
    # final arrays
    ft=np.array(df['time'].values[:]).astype('str'); fbt = np.array(df['obs_time'].values[:]).astype('str')
    fid=df['id'].values[:]; fcid=df['cid'].values[:]
    fcmap=df['cmap'].values[:]; fcsec=df['csec'].values[:]; 
    fccdist=df['ccdist'].values[:]; fccangle=df['ccangle'].values[:]
    # flat=np.zeros(len(ot),'float')-999; flon=np.zeros(len(ot),'float')-999
    # obs arrays
    ohs=df['hs'].values[:]; otm=df['tm'].values[:]; otp=df['tp'].values[:]; ownd=df['wnd'].values[:]

    # PModel arrays
    mhs=np.zeros(len(ot),'float')-999.; mtp=np.zeros(len(ot),'float')-999.; muwnd=np.zeros(len(ot),'float')-999.; mvwnd=np.zeros(len(ot),'float')-999.
    # Read PModel and allocate data
    for y in [2022,2023,2024,2025]:
        for m in [6,7,8,9,10,11]:
            f = nc.Dataset("Pmodel_reanalysis_"+str(y)+str(m).zfill(2)+"_Default.nc")
            # latm = f.variables['lat'][:]; lonm = f.variables['lon'][:]
            at = f.variables['time'][:]
            ahs = f.variables['hs'][:,:,:]; atp = f.variables['tp'][:,:,:]
            auwnd = f.variables['uwnd'][:,:,:]; avwnd = f.variables['vwnd'][:,:,:]
            for i in range(0,len(ot)):
                if np.min(np.abs( ot[i] - at )) < 1800.:
                    # Index array of PModel data
                    indt = np.min(np.where( np.abs( ot[i] - at ) == np.min(np.abs( ot[i] - at )) )[0])
                    # Allocate model data
                    mhs[i] = float(ahs[indt,gidlat[i],gidlon[i]])
                    mtp[i] = float(atp[indt,gidlat[i],gidlon[i]])
                    muwnd[i] = float(auwnd[indt,gidlat[i],gidlon[i]])
                    mvwnd[i] = float(avwnd[indt,gidlat[i],gidlon[i]])
                    # flat[i] = float(latm[gidlat[i]]); flon[i] = float(lonm[gidlon[i]])
                    del indt

            f.close(); del f
            del at,ahs,atp,auwnd,avwnd # latm,lonm
            print(" Ok "+repr(y)+str(m).zfill(2))

    mhs[np.isnan(mhs)==True]=-999.; mtp[np.isnan(mtp)==True]=-999.
    muwnd[mhs<0.]=-999.; mvwnd[mhs<0.]=-999.;
    muwnd[np.isnan(muwnd)==True]=-999.; mvwnd[np.isnan(mvwnd)==True]=-999.

    # Save results
    df = pd.DataFrame({
        'time': ft,
        'obs_time': fbt,
        'lat': np.round(flat,4),
        'lon': np.round(flon,4),
        'id': fid,
        'cmap': fcmap,
        'csec': fcsec,
        'ccdist': fccdist,
        'ccangle': fccangle,
        'cid': fcid,
        'obs_hs': np.round(ohs,3),
        'obs_tm': np.round(otm,3),
        'obs_tp': np.round(otp,3),
        'obs_wnd': np.round(ownd,3),
        'pm_hs': np.round(mhs,3),
        'pm_tp': np.round(mtp,3),
        'pm_uwnd': np.round(muwnd,3),
        'pm_vwnd': np.round(mvwnd,3),
    })
    df.to_csv("Data_Obs_PModel_"+pmrun+".txt", sep='\t', index=False, header=True)

