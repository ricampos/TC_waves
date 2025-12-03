#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import numpy as np
import pandas as pd
import netCDF4 as nc
import timeit
import calendar
import warnings; warnings.filterwarnings("ignore")

if __name__ == "__main__":

    membr=int(sys.argv[1])

    gpath = "/work/noaa/marine/ricardo.campos/data/archiveOPruns/GEFSv12Waves_AWS/netcdf/tc_waves"

    # start time
    start = timeit.default_timer()

    # Read Obs
    df = pd.read_csv('Data_Obs_PModel_Default.txt', sep='\t')
    ot = np.array(df['time']).astype('str')
    ot = pd.to_datetime(ot, format='%Y%m%d%H%M'); ot = ot.astype("int64") // 1_000_000_000
    lat = np.array(df['lat']); lon = np.array(df['lon']); lon[lon<0]=lon[lon<0]+360.

    # GEFSv12 sample
    f = nc.Dataset(gpath+"/gefs.wave.2023113012.23.global.0p25.nc")
    glat = f.variables['latitude'][:]; glon = f.variables['longitude'][:]
    f.close(); del f

    gvarn = np.array(['UGRD_surface','VGRD_surface','HTSGW_surface','PERPW_surface','IMWF_surface','MWSPER_surface','DIRPW_surface','SWELL_1insequence','SWELL_2insequence','SWELL_3insequence','SWPER_1insequence','SWPER_2insequence','SWPER_3insequence',
        'SWDIR_1insequence','SWDIR_2insequence','SWDIR_3insequence','WVHGT_surface','WVPER_surface','WWSDIR_surface'])

    gres = np.zeros((len(ot),3+len(gvarn)),dtype=object)-999.

    for y in [2022,2023,2024]:
        for m in [6,7,8,9,10,11]:
            for d in range(1,calendar.monthrange(y, m)[1]+1):
                for c in [0,12]:

                    f = nc.Dataset(gpath+"/gefs.wave."+repr(y)+str(m).zfill(2)+str(d).zfill(2)+str(c).zfill(2)+"."+str(membr).zfill(2)+".global.0p25.nc")
                    t = f.variables['time'][0:4]

                    for i in range(0,len(t)):
                        ind = np.where(  np.abs(t[i]-ot) <= 5400. )
                        if np.size(ind)>0:
                            ind=ind[0]
                            for j in range(0,len(ind)):
                                gres[ind[j],0] = int(pd.to_datetime( ot[ind[j]] , unit="s").strftime("%Y%m%d%H"))

                                indlat = np.min(np.where( np.abs(glat-lat[ind[j]]) == np.min(np.abs(glat-lat[ind[j]])) )[0])
                                indlon = np.min(np.where( np.abs(glon-lon[ind[j]]) == np.min(np.abs(glon-lon[ind[j]])) )[0])
                                gres[ind[j],1] = float(glat[indlat])
                                gres[ind[j],2] = float(glon[indlon])

                                for k in range(0,len(gvarn)):
                                    gaux = float(f.variables[gvarn[k]][i,indlat,indlon])
                                    if f.variables[gvarn[k]][i,indlat,indlon] > -999.:
                                        gres[ind[j],3+k] = gaux
                                    else:
                                        gres[ind[j],3+k] = float(-999.)

                                    del gaux

                                del indlat,indlon

                        del ind

                    f.close(); del f, t

                    print("ok "+repr(y)+str(m).zfill(2)+str(d).zfill(2)+str(c).zfill(2)+"."+str(membr).zfill(2))


    print("Done")
    df = pd.DataFrame(gres, columns=np.append(['time','lat','lon'],gvarn))
    # df_rounded = df.round(4)
    for col in df.columns:
        if col != "time":
            df[col] = pd.to_numeric(df[col], errors="coerce")
            df[col] = df[col].round(5)

    df.to_csv("Data_GEFS_fromObs_m"+str(membr).zfill(2)+".txt", sep="\t", index=False, header=True)

    stop = timeit.default_timer()
    print('Concluded in '+repr(int(round(stop - start,0)))+' seconds')

