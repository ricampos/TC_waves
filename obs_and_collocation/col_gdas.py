#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
# from matplotlib.mlab import *
import pandas as pd
# from pylab import *
import os
import netCDF4 as nc
import xarray as xr
import time
import timeit
from calendar import timegm
import sys
import warnings; warnings.filterwarnings("ignore")
# netcdf format
fnetcdf="NETCDF4"

if __name__ == "__main__":

    # Inputs, number of split segments. Segment selected.
    nseg=int(sys.argv[1])
    seg=int(sys.argv[2])
    maxti=1800.
    # ------------------------

    start = timeit.default_timer()
    print(" Processing segment "+str(seg)+" of total "+str(nseg))

    # Read buoys Reference
    df = pd.read_csv('Data_REF.txt', sep='\t')
    lines_to_read = len(df) // nseg
    if seg<nseg:
        df = df.iloc[(seg-1)*lines_to_read:(seg)*lines_to_read, :]
    else:
        df = df.iloc[(seg-1)*lines_to_read::, :]
    # bhs = np.array(df['hs'][:]); bhs[bhs<0.1]=np.nan
    blat = np.array(df['lat'][:]); blon = np.array(df['lon'][:])
    btime = np.array( (pd.to_datetime(df['time'][:], format='%Y%m%d%H%M') - pd.Timestamp('1970-01-01')) // pd.Timedelta('1s') ).astype('double')
    del df
    # ---------------

    # Read GDAS
    ds = xr.open_mfdataset("/home/ricardo/work/noaa/analysis/TC_Waves/data/GDAS/gdaswave.*.global.0p16.nc", combine='by_coords')
    gtime = np.array( (pd.to_datetime(ds['time'][:], format='%Y%m%d%H') - pd.Timestamp('1970-01-01')) // pd.Timedelta('1s') ).astype('double')
    glat = np.array(ds['latitude'][:]).astype('float'); glon = np.array(ds['longitude'][:]).astype('float'); glon[glon>180]=glon[glon>180]-360.
    ghs = ds['HTSGW_surface']; gwnd = ds['WIND_surface']
    # ---------------

    fhs=np.zeros((btime.shape[0]),'f')*np.nan
    fwnd=np.zeros((btime.shape[0]),'f')*np.nan
    for t in range(0,btime.shape[0]):

        indt = np.where( abs(gtime[:]-btime[t]) < maxti )
        if np.size(indt)>0:

            indlat = np.where( np.abs(glat-blat[t]) == np.nanmin(np.abs(glat-blat[t])) )[0][0]
            indlon = np.where( np.abs(glon-blon[t]) == np.nanmin(np.abs(glon-blon[t])) )[0][0]

            ahs = float(ghs[indt[0][0],indlat,indlon].values)
            if ahs>0:
                fhs[t]=ahs

            awnd = float(gwnd[indt[0][0],indlat,indlon].values)
            if awnd>0:
                fwnd[t]=awnd

        print(repr(t)+" of "+repr(btime.shape[0]))


    fhs[np.isnan(fhs)==True]=-999.999
    fwnd[np.isnan(fwnd)==True]=-999.999

    # Save final data
    df = pd.DataFrame({
        'time': pd.to_datetime(btime, unit='s').strftime('%Y%m%d%H'),
        'lat': np.round(blat,5),
        'lon': np.round(blon,5),
        'hs': np.round(fhs,4),
        'wnd': np.round(fwnd,4)
    })

    fname="ColData_GDAS_"+str(seg).zfill(2)+"_"+str(nseg).zfill(2)+".txt"
    if seg==1:
        df.to_csv(fname, sep='\t', index=False, header=True)
    else:
        df.to_csv(fname, sep='\t', index=False, header=False)

    stop = timeit.default_timer()
    print('Concluded in '+repr(int(round(stop - start,0)))+' seconds')

