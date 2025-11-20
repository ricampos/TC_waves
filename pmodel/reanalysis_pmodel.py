#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
 Script to run a reanalysis of tropical cyclones using PModel, for the Atlantic and Pacific domains.
 It is based on ibtracs tracks. Results are produced for a high resolution 5km grid (gridInfo.nc).
 See TCpwaves.py and pmodel_config.yaml for the parametric modeling
 This code is run for one month. Input arguments are: year month string_tag
 Output parameters are significant wave height (Hs), peak period (Tp), and wind (U and V components)
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
import time
import timeit
import sys
import pandas as pd
from matplotlib.dates import DateFormatter
import TCpwaves
from TCpwaves import *
fnetcdf="NETCDF4_CLASSIC"

if __name__ == "__main__":

    wdname=str('pmodel_config.yaml')

    # Input parameters, year and month. And string to identify the configuration/calibration experiment (yaml)
    cyear=int(sys.argv[1])
    cmonth=int(sys.argv[2])
    ctag=str(sys.argv[3])

    start = timeit.default_timer()

    # Land/sea mask, water depth and distance to the coast.
    f=nc.Dataset('gridInfo.nc')
    latp=f.variables['latitude'][:]; mres=np.diff(latp).mean()
    lonp=f.variables['longitude'][:]; # lonp[lonp>180]=lonp[lonp>180]-360.
    mask=f.variables['mask'][:]
    # Grid where the parametric cyclone will be generated
    FHs = np.zeros(mask.shape,'f'); FHs[np.isnan(mask)==True]=np.nan
    FTp = np.zeros(mask.shape,'f'); FTp[np.isnan(mask)==True]=np.nan
    FUwnd = np.zeros(mask.shape,'f'); FUwnd[np.isnan(mask)==True]=np.nan
    FVwnd = np.zeros(mask.shape,'f'); FVwnd[np.isnan(mask)==True]=np.nan
    W={'lat':latp,'lon':lonp,'Hs':FHs,'Tp':FTp,'Uwnd':FUwnd,'Vwnd':FVwnd}

    # Read ibtracks cyclone data
    wibtr = read_ibtracs('ibtracs.last3years.list.v04r01.csv')

    # Cyclones inside the domain of interest
    ailon = np.copy(wibtr['ilon']); ailon[ailon<0]=ailon[ailon<0]+360.
    ind = np.where( (wibtr['ilat']>(latp.min()+1)) & (wibtr['ilat']<(latp.max()-1)) & (ailon>(lonp.min()+1)) & (ailon<(lonp.max()-1)))[0]
    for k in wibtr:
        wibtr[k] = wibtr[k][ind]

    ctime = (wibtr['itime'].values - np.datetime64('1970-01-01T00:00:00')) / np.timedelta64(1, 's')

    # Plot cyclone tracks
    # ibtcyclones = np.unique(wibtr['icyid'])
    # for i in range(0,len(ibtcyclones)):
    #    indt = np.where( wibtr['icyid'] == ibtcyclones[i] )[0]
    #    # plot cyclone track
    #    ctrack(wibtr['ilat'][indt],wibtr['ilon'][indt],wibtr['itime'][indt],wibtr['iVmax'][indt],flabel=wibtr['iname'][indt][0],ftitle=wibtr['iname'][indt][0])

    # --------------------
    # Time: 3 years. Jun,Jul,Aug,Sep,Oct,Nov.  ~240 steps per month
    tstart = np.datetime64(f'{cyear}-{cmonth:02d}-01T00:00:00'); tend = np.datetime64(f'{cyear}-{cmonth+1:02d}-01T00:00:00')
    time_array = np.arange(tstart, tend, np.timedelta64(3, 'h'))
    ftime = (time_array - np.datetime64('1970-01-01T00:00:00')) / np.timedelta64(1, 's')
    ftime = ftime.astype('float64'); del tstart, tend, time_array

    # final arrays
    fhs=np.zeros((len(ftime),mask.shape[0],mask.shape[1]),'f')*np.nan
    ftp=np.zeros((len(ftime),mask.shape[0],mask.shape[1]),'f')*np.nan
    fuwnd=np.zeros((len(ftime),mask.shape[0],mask.shape[1]),'f')*np.nan
    fvwnd=np.zeros((len(ftime),mask.shape[0],mask.shape[1]),'f')*np.nan

    for t in range(0,len(ftime)):
        indt = np.where( ctime == ftime[t] )[0]
        if np.any(indt):
            # in case of more than one cyclone at the same time, loop through cyclones
            ahs=np.zeros((len(indt),mask.shape[0],mask.shape[1]),'f')
            atp=np.zeros((len(indt),mask.shape[0],mask.shape[1]),'f')
            auwnd=np.zeros((len(indt),mask.shape[0],mask.shape[1]),'f')
            avwnd=np.zeros((len(indt),mask.shape[0],mask.shape[1]),'f')

            for i in range(0,len(indt)):

                # heading directions
                indc = np.where(wibtr['icyid'] == wibtr['icyid'][indt][i])[0]
                rangle = cbearing(wibtr['ilat'][indc],wibtr['ilon'][indc])
                rangle = rangle[ np.nanmin(np.where( wibtr['itime'][indt][i]==wibtr['itime'][indc] )) ]
                del indc

                tcw = TCWaves(Vmax=wibtr['iVmax'][indt][i],Vfm=wibtr['iVfm'][indt][i],Rmax=wibtr['iRmax'][indt][i],R34=wibtr['iR34'][indt][i],
                    Lat=wibtr['ilat'][indt][i],Lon=wibtr['ilon'][indt][i],wdname=wdname)

                pwm = tcw.PWModel()
                pwm = tcw.xy_to_lonlat(pwm)
                pwm = tcw.rotate(pwm,rangle)
                pwm,WW = tcw.blend(pwm,W)

                if np.any(np.where(WW['Hs']>0.1)):
                    ahs[i,:,:] = np.array(WW['Hs'])
                    atp[i,:,:] = np.array(WW['Tp'])
                    auwnd[i,:,:] = np.array(WW['Uwnd'])
                    avwnd[i,:,:] = np.array(WW['Vwnd'])

                del pwm,WW,tcw,rangle,iVmax

            if np.any(np.where(np.nanmax(ahs[:,:,:],axis=0)>0.1)):
                fhs[t,:,:] = np.array(np.nanmax(ahs[:,:,:],axis=0))
                ftp[t,:,:] = np.array(np.nanmax(atp[:,:,:],axis=0))
                fuwnd[t,:,:] = np.array(np.nanmean(auwnd[:,:,:],axis=0))
                fvwnd[t,:,:] = np.array(np.nanmean(avwnd[:,:,:],axis=0))

            del ahs,atp,auwnd,avwnd

        print(repr(t))

    fhs[fhs<0.1]=np.nan; ftp[ftp<0.1]=np.nan
    fuwnd[fuwnd<-999.]=np.nan; fvwnd[fvwnd<-999.]=np.nan
    fhs[:,np.isnan(mask)]=np.nan; ftp[:,np.isnan(mask)]=np.nan
    fuwnd[:,np.isnan(mask)]=np.nan; fvwnd[:,np.isnan(mask)]=np.nan

    ttag=str(time.gmtime(ftime.min())[0])+str(time.gmtime(ftime.min())[1]).zfill(2)

    # Save netcdf
    ncfile = nc.Dataset("Pmodel_reanalysis_"+ttag+"_"+ctag+".nc", "w", format=fnetcdf) 
    ncfile.history='Pmodel reanalysis using cyclone information from IBTracks. The .yaml has the parametric model parameters.'
    # create  dimensions.
    ncfile.createDimension('time', ftime.shape[0])
    ncfile.createDimension('lat', latp.shape[0])
    ncfile.createDimension('lon', lonp.shape[0])
    # Model results
    vhs = ncfile.createVariable('hs',np.dtype('float32'),('time','lat','lon'))
    vtp = ncfile.createVariable('tp',np.dtype('float32'),('time','lat','lon'))
    vuwnd = ncfile.createVariable('uwnd',np.dtype('float32'),('time','lat','lon'))
    vvwnd = ncfile.createVariable('vwnd',np.dtype('float32'),('time','lat','lon'))
    #
    vt = ncfile.createVariable('time',np.dtype('float64'),('time'))
    vlat = ncfile.createVariable('lat',np.dtype('float32'),('lat',))
    vlon = ncfile.createVariable('lon',np.dtype('float32'),('lon'))
    # Units
    vlat.units = 'degrees_north' ; vlon.units = 'degrees_east'
    vt.units = "seconds since 1970-01-01 00:00:00.0 0:00"
    # Allocate Data
    vt[:] = ftime[:]; vlat[:] = latp[:]; vlon[:] = lonp[:]
    vhs[:,:,:] = fhs[:,:,:]; vtp[:,:,:] = ftp[:,:,:]; vuwnd[:,:,:] = fuwnd[:,:,:]; vvwnd[:,:,:] = fvwnd[:,:,:]
    ncfile.close()
    print('netcdf ok ')

    stop = timeit.default_timer()
    print('Concluded in '+repr(int(round(stop - start,0)))+' seconds')

