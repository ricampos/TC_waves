#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import netCDF4 as nc

if __name__ == "__main__":

    pmrun="Default"

    # Read Obs
    df = pd.read_csv('Data_TC.txt', sep='\t')
    ot = np.array(df['time']).astype('str')
    ot = pd.to_datetime(ot, format='%Y%m%d%H%M'); ot = np.array(ot.view('int64') // 1_000_000_000).astype('double')
    gidlat = np.array(df['gidlat']); gidlon = np.array(df['gidlon'])

    # Read PModel
    ft=[]; flat=[]; flon=[]; fid=[]; fcmap=[]; fcsec=[]; fcid=[]; ohs=[]; otm=[]; otp=[]; ownd=[]
    mhs=[]; mtp=[]; muwnd=[]; mvwnd=[]
    for y in [2022,2023,2024]:
        for m in [6,7,8,9,10,11]:
            f = nc.Dataset("Pmodel_reanalysis_"+str(y)+str(m).zfill(2)+"_Default.nc")
            at = f.variables['time'][:]
            ahs = f.variables['hs'][:,:,:]
            atp = f.variables['tp'][:,:,:]
            auwnd = f.variables['uwnd'][:,:,:]
            avwnd = f.variables['vwnd'][:,:,:]
            # hs = ds.hs.load()
            idt=[];ilat=[];ilon=[]
            for i in range(0,len(ot)):
                if np.min(np.abs( ot[i] - at )) < 1800.:
                    # Index array of PModel data
                    indt = np.min(np.where( np.abs( ot[i] - at ) == np.min(np.abs( ot[i] - at )) )[0])
                    idt=np.append(idt,indt)
                    ilat=np.append(ilat,gidlat[i])
                    ilon=np.append(ilon,gidlon[i])
                    del indt

                    # Allocate and keep observation data
                    ft=np.append(ft,ot[i]); fid=np.append(fid,df['id'][i])
                    flat=np.append(flat,df['lat'][i]); flon=np.append(flon,df['lon'][i])
                    fcmap=np.append(fcmap,df['cmap'][i]); fcsec=np.append(fcsec,df['csec'][i]); fcid=np.append(fcid,df['cid'][i])
                    ohs=np.append(ohs,df['hs'][i]); otm=np.append(otm,df['tm'][i]); otp=np.append(otp,df['tp'][i]); ownd=np.append(ownd,df['wnd'][i])

            # Apply index to retrieve PModel data, for observation points (time and location)
            idt=np.array(idt).astype('int')
            ilat=np.array(ilat).astype('int')
            ilon=np.array(ilon).astype('int')

            mhs = np.append(mhs,ahs[idt,ilat,ilon]); mtp = np.append(mtp,atp[idt,ilat,ilon])
            muwnd = np.append(muwnd,auwnd[idt,ilat,ilon]); mvwnd = np.append(mvwnd,avwnd[idt,ilat,ilon])

            f.close(); del f
            del idt,ilat,ilon,at,ahs,atp,auwnd,avwnd
            print(" Ok "+repr(y)+str(m).zfill(2))

    ind=np.where( (mhs>0.2) & (ohs>0.2) )[0]

    # Save results
    df = pd.DataFrame({
        'time': pd.to_datetime(ft[ind], unit='s').strftime('%Y%m%d%H%M'),
        'lat': np.round(flat[ind],4),
        'lon': np.round(flon[ind],4),
        'id': fid[ind],
        'cmap': np.array(fcmap[ind]).astype('int'),
        'csec': np.array(fcsec[ind]).astype('int'),
        'cid': np.array(fcid[ind]).astype('int'),
        'obs_hs': np.round(ohs[ind],3),
        'obs_tm': np.round(otm[ind],3),
        'obs_tp': np.round(otp[ind],3),
        'obs_wnd': np.round(ownd[ind],3),
        'pm_hs': np.round(mhs[ind],3),
        'pm_tp': np.round(mtp[ind],3),
        'pm_uwnd': np.round(muwnd[ind],3),
        'pm_vwnd': np.round(mvwnd[ind],3),
    })
    df.to_csv("Data_Obs_PModel_"+pmrun+".txt", sep='\t', index=False, header=True)

