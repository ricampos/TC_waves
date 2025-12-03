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

    # start time
    start = timeit.default_timer()

    hpath = "/work/noaa/marine/ricardo.campos/data/archiveOPruns/HAFS_AWS/netcdf"

    hlist = np.loadtxt(hpath + "/list.txt", dtype=str)

    # Read Obs
    df = pd.read_csv('Data_Obs_PModel_Default.txt', sep='\t')
    ot = np.array(df['time']).astype('str')
    ot = pd.to_datetime(ot, format='%Y%m%d%H%M'); ot = ot.astype("int64") // 1_000_000_000
    lat = np.array(df['lat']); lon = np.array(df['lon']); lon[lon<0]=lon[lon<0]+360.

    hvarn = np.array(['UGRD_surface','VGRD_surface','HTSGW_surface','PERPW_surface','MWSPER_surface','DIRPW_surface','SWELL_1insequence','SWELL_2insequence','SWELL_3insequence','SWPER_1insequence','SWPER_2insequence','SWPER_3insequence',
        'SWDIR_1insequence','SWDIR_2insequence','SWDIR_3insequence','WVHGT_surface','WVPER_surface','WWSDIR_surface'])

    hres = np.zeros((len(ot),3+len(hvarn)),dtype=object)-999.

    for l in range(0,len(hlist)):

        f = nc.Dataset(hpath+"/"+hlist[l])
        t = f.variables['time'][0:4]
        hlat = f.variables['latitude'][:]; hlon = f.variables['longitude'][:]; hlon[hlon<0]=hlon[hlon<0]+360.
        for i in range(0,len(t)):
            ind = np.where(  (np.abs(t[i]-ot)<= 5400.) & (lat>=hlat.min()) & (lat<=hlat.max()) & (lon>=hlon.min()) & (lon<=hlon.max()) )
            if np.size(ind)>0:
                ind=ind[0]
                for j in range(0,len(ind)):

                    hres[ind[j],0] = int(pd.to_datetime( ot[ind[j]] , unit="s").strftime("%Y%m%d%H"))

                    indlat = np.min(np.where( np.abs(hlat-lat[ind[j]]) == np.min(np.abs(hlat-lat[ind[j]])) )[0])
                    indlon = np.min(np.where( np.abs(hlon-lon[ind[j]]) == np.min(np.abs(hlon-lon[ind[j]])) )[0])
                    hres[ind[j],1] = float(hlat[indlat])
                    hres[ind[j],2] = float(hlon[indlon])

                    for k in range(0,len(hvarn)):
                        haux = float(f.variables[hvarn[k]][i,indlat,indlon])
                        if f.variables[hvarn[k]][i,indlat,indlon] > -999.:
                            hres[ind[j],3+k] = haux
                        else:
                            hres[ind[j],3+k] = float(-999.)

                        del haux

                    del indlat,indlon

                del ind

        f.close(); del f, t

        print("ok "+hlist[l])


    print("Done")
    df = pd.DataFrame(hres, columns=np.append(['time','lat','lon'],hvarn))
    # df_rounded = df.round(4)
    for col in df.columns:
        if col != "time":
            df[col] = pd.to_numeric(df[col], errors="coerce")
            df[col] = df[col].round(5)

    df.to_csv("Data_HAFS_fromObs.txt", sep="\t", index=False, header=True)

    stop = timeit.default_timer()
    print('Concluded in '+repr(int(round(stop - start,0)))+' seconds')

