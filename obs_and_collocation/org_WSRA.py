#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
org_WSRA.py

VERSION AND LAST UPDATE:
 v1.0  10/09/2025

PURPOSE:
 Data processing, quality control, and organizing of wave obs:
  WSRA

OUTPUT:
 Text file Data_*.txt. See the last line df.to_csv(

DEPENDENCIES:
 wread.py and quality_control_wave.py
 See the imports below.

AUTHOR and DATE:
 10/09/2025: Ricardo M. Campos, first version.

PERSON OF CONTACT:
 Ricardo M Campos: ricardo.campos@noaa.gov

"""

import math
from xarray import open_dataset
import netCDF4 as nc
import numpy as np
from pylab import *
import pandas as pd
import time
from wread import *
import quality_control_wave


if __name__ == "__main__":

    # Input settings
    # Data intervall of the final array (original time is also saved)
    wdt = 3600.
    # Buoy data type
    buoyd="WSRA"
    # data path
    dpath="/data/WSRA_L4"
    events = np.array(pd.read_csv(dpath+"/list_events_2024.txt", header=None).values).astype('str')[:,0]

    # GridMask
    f=nc.Dataset('gridInfo_TGPM.nc')
    latm=f.variables['latitude'][:]; lonm=f.variables['longitude'][:]; lonm[lonm>180]=lonm[lonm>180]-360.
    maskm=f.variables['mask'][:,:]; depthm=f.variables['depth'][:,:]; depthm=f.variables['depth'][:,:]
    f.close(); del f

    # Cyclone Info
    f = nc.Dataset('CycloneMap.nc')
    latc=f.variables['lat'][:]; lonc=f.variables['lon'][:]; lonc[lonc>180]=lonc[lonc>180]-360.
    ctime = np.array(f.variables['time'][:]).astype('double')
    cmap = np.array(f.variables['cmap'][:,:,:]).astype('float')
    csec = np.array(f.variables['csec'][:,:,:]).astype('float')
    cid = np.array(f.variables['cid'][:,:,:]).astype('float')
    f.close(); del f
    # cmap[cmap<0]=np.nan; cid[cid<0]=np.nan; csec[csec<0]=np.nan

    ftime=[]; frtime=[]; bid=[]
    lat=[]; lon=[]; gidlat=[]; gidlon=[] 
    hs=[]; tp=[]; tm=[]; wnd=[] 
    bcmap=[]; bcid=[]; bcsec=[]

    for i in range(0,len(events)):
        dname = dpath+"/"+events[i]
        fnames = np.array(pd.read_csv(dname+"/list.txt", header=None).values).astype('str')[:,0]
        for j in range(0,len(fnames)):
            # Read netcdf file
            fname=dpath+"/"+events[i]+"/"+fnames[j]

            try:
                wdic = tseriesnc_wsra(fname)
            except:
                print(" Cannot open "+fname)
            else:

                print(" Processing QC for "+fname)
                wdic['hs']=quality_control_wave.data_range(wdic,var='hs',vmin=0.3,vmax=20.)
                wdic['hs']=quality_control_wave.duplicates(wdic)
                wdic['hs']=quality_control_wave.rate_of_change(wdic)
                wdic['hs']=quality_control_wave.landcoast_exclude(wdic,gpath='gridInfo_TGPM.nc',mdepth=80,mdfc=5)
                wdic=quality_control_wave.wsra(wdic)
                wdic['hs']=quality_control_wave.model_compare(wdic,gpath='/data/GDAS',mdist=None)
                print(" OK - QC for "+fname+" "+fname)

                ind=np.where(wdic['hs']>0.1)
                if np.size(ind)>0:

                    ind=ind[0]

                    at=wdic['time'][ind]; adate=wdic['date'][ind]
                    ahs=wdic['hs'][ind];

                    if "wind_spd" in wdic:
                        awnd=wdic['wind_spd'][ind]
                    else:
                        awnd=np.zeros(len(ahs),'f')-999.999

                    if "tp" in wdic:
                        atp=wdic['tp'][ind]
                    else:
                        atp=np.zeros(len(ahs),'f')-999.999

                    if "tm" in wdic:
                        atm=wdic['tm'][ind]
                    else:
                        atm=np.zeros(len(ahs),'f')-999.999

                    alat=wdic['latitude'][ind]; alon=wdic['longitude'][ind]

                    initime = str(pd.to_datetime(adate.min()).year)+str(pd.to_datetime(adate.min()).month).zfill(2)+str(pd.to_datetime(adate.min()).day).zfill(2)+str(pd.to_datetime(adate.min()).hour).zfill(2)
                    fintime = str(pd.to_datetime(adate.max()).year)+str(pd.to_datetime(adate.max()).month).zfill(2)+str(pd.to_datetime(adate.max()).day).zfill(2)+str(pd.to_datetime(adate.max()).hour).zfill(2)
                    aftime = np.array(np.arange(float(timegm( time.strptime(initime, '%Y%m%d%H') )),float(timegm( time.strptime(fintime, '%Y%m%d%H') ))+1,wdt)).astype('double')

                    blat=[]; blon=[]; aindt=np.array([]).astype('int')
                    bhs=[]; bwnd=[]; btp=[]; btm=[]
                    btime=np.double([]); brt=np.double([])
                    abcmap=[]; abcid=[]; abcsec=[]
                    agidlat=[]; agidlon=[]; abid=[]
                    c=0
                    for t in range(0,len(at)):

                        # organize time and allocate data
                        indt=np.where(np.abs(at[t]-aftime)<=1800.)
                        if np.size(indt)>0:

                            bhs = np.append(bhs,ahs[t])
                            bwnd = np.append(bwnd,awnd[t])
                            btp = np.append(btp,atp[t])
                            btm = np.append(btm,atm[t])
                            btime = np.append(btime,double(np.min(at[t])))
                            brt = np.append(brt,double(np.min(aftime[indt])))
                            blat = np.append(blat,alat[t])
                            blon = np.append(blon,alon[t])
                            del indt

                            # Model position index
                            indlat = np.where( np.abs(latm-alat[t]) == np.nanmin(np.abs(latm-alat[t])) )[0][0]
                            indlon = np.where( np.abs(lonm-alon[t]) == np.nanmin(np.abs(lonm-alon[t])) )[0][0]
                            agidlat = np.append(agidlat,int(indlat)); agidlon = np.append(agidlon,int(indlon))
                            del indlat, indlon

                            # check cyclone presence
                            indc=np.where(np.abs(ctime-at[t])<=5400.)
                            if np.size(indc)>0:
                                indlat = np.where( np.abs(latc-alat[t]) == np.nanmin(np.abs(latc-alat[t])) )[0][0]
                                indlon = np.where( np.abs(lonc-alon[t]) == np.nanmin(np.abs(lonc-alon[t])) )[0][0]
                                if cmap[np.min(indc[0]),indlat,indlon]>0:
                                    abcmap = np.append(abcmap,int(cmap[np.min(indc[0]),indlat,indlon]))
                                    abcid = np.append(abcid,int(cid[np.min(indc[0]),indlat,indlon]))
                                    abcsec = np.append(abcsec,int(csec[np.min(indc[0]),indlat,indlon]))
                                else:
                                    abcmap = np.append(abcmap,0)
                                    abcid = np.append(abcid,0)
                                    abcsec = np.append(abcsec,0)

                                del indlat,indlon

                            else:
                                abcmap = np.append(abcmap,0)
                                abcid = np.append(abcid,0)
                                abcsec = np.append(abcsec,0)

                            del indc

                        print(" Ok data allocation "+repr(t)+" "+fname)

                    del ind

                    # Final arrays
                    if np.size(abcid)>0:
                        abot="WSRA"+events[i]+fname[-5:-3]
                        aux=np.full(btime.shape[0], abot, dtype=f'<U{len(abot)}')
                        abid = np.append(abid,aux); del aux, abot

                        btime=np.array(btime).astype('double'); brt=np.array(brt).astype('double')
                        agidlat=np.array(agidlat).astype('int'); agidlon=np.array(agidlon).astype('int')

                        ftime=np.append(ftime,btime); frtime=np.append(frtime,brt)
                        lat=np.append(lat,blat); lon=np.append(lon,blon); bid=np.append(bid,abid)
                        # glat=np.append(glat,bglat); glon=np.append(glon,bglon)
                        gidlat=np.append(gidlat,agidlat); gidlon=np.append(gidlon,agidlon)
                        hs=np.append(hs,bhs); tp=np.append(tp,btp); tm=np.append(tm,btm); wnd=np.append(wnd,bwnd)
                        bcmap = np.append(bcmap,abcmap)
                        abcid[abcid<0]=0.; abcid[np.isnan(abcid)==True]=0.; bcid = np.append(bcid,abcid)
                        bcsec = np.append(bcsec,abcsec)

                    print(fname+" done")

                    del bhs,btp,btm,bwnd,btime,brt,aftime,blat,blon,abcmap,abcid,abid #,bglat,bglon

            del wdic


    # Save results
    ind=np.where(hs>0.01)
    if np.size(ind)>0:

        hs[np.isnan(hs)==True]=-999.999; wnd[np.isnan(wnd)==True]=-999.999
        tp[np.isnan(tp)==True]=-999.999; tm[np.isnan(tm)==True]=-999.999

        hs=np.round(hs,4); tp=np.round(tp,4); tm=np.round(tm,4); wnd=np.round(wnd,4)
        lat=np.round(lat,5); lon=np.round(lon,5); # glat=np.round(glat,5); glon=np.round(glon,5)
        gidlat=np.array(gidlat).astype('int'); gidlon=np.array(gidlon).astype('int')
        bcmap=np.array(bcmap).astype('int'); bcid=np.array(bcid).astype('int'); bcsec=np.array(bcsec).astype('int')

        # Save wdics
        df = pd.DataFrame({
            'time': pd.to_datetime(frtime, unit='s').strftime('%Y%m%d%H%M'),
            'buoy_time': pd.to_datetime(ftime, unit='s').strftime('%Y%m%d%H%M'),
            'lat': lat,
            'lon': lon,
            # 'glat': glat,
            # 'glon': glon,
            'gidlat': gidlat,
            'gidlon': gidlon,
            'id': bid,
            'cmap': bcmap,
            'csec': bcsec,
            'cid': bcid,
            'hs': hs,
            'tm': tm,
            'tp': tp,
            'wnd': wnd
        })

        df.to_csv('Data_WSRA.txt', sep='\t', index=False, header=True)

