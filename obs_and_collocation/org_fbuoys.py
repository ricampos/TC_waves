#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
org_fbuoys.py

VERSION AND LAST UPDATE:
 v1.0  08/28/2025

PURPOSE:
 Data processing, quality control, and organizing of wave (or metocean) buoys (fixed):
  CDIP, NDBC

OUTPUT:
 Text file Data_*.txt. See the last line df.to_csv(

DEPENDENCIES:
 wread.py and quality_control_wave.py
 See the imports below.

AUTHOR and DATE:
 08/28/2025: Ricardo M. Campos, first version.

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
    buoyd="CDIP"

    # List buoy data
    if buoyd=="NDBC":
        dpath="/data/NDBC/wparam"
        bnames = np.array(pd.read_csv(dpath+"/list.txt", header=None).values).astype('str')
    else:
        dpath="/data/CDIP"
        bnames = np.array(pd.read_csv(dpath+"/CDIP_buoy_selection.txt", header=None).values).astype('str')

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
    lat=[]; lon=[]; gidlat=[]; gidlon=[]; # glat=[]; glon=[]; 
    hs=[]; tp=[]; tm=[]; wnd=[] 
    bcmap=[]; bcid=[]; bcsec=[]

    for i in range(0,len(bnames)):

        # Read netcdf file
        if buoyd=="NDBC":
            fname=dpath+"/"+bnames[i][0]
        else:
            fname=dpath+"/CDIP_buoy_"+bnames[i][0]+"_historic.nc" # CDIP

        try:
            if buoyd=="NDBC":
                wdic = tseriesnc_ndbc(fname)
            else:
                wdic = tseriesnc_cdip(fname)

        except:
            print(" Cannot open "+fname)
        else:

            if np.any(wdic['hs']>0.1):

                print(" Processing QC for "+bnames[i][0])
                wdic['hs']=quality_control_wave.data_range(wdic,var='hs',vmin=0.3,vmax=20.)
                wdic['hs']=quality_control_wave.duplicates(wdic)
                wdic['hs']=quality_control_wave.rate_of_change(wdic)
                wdic['hs']=quality_control_wave.landcoast_exclude(wdic,gpath='/2collocation/gridInfo_TGPM.nc',mdepth=80,mdfc=5)
                # wdic['hs']=quality_control_wave.model_compare(wdic,gpath=None,mdist=None)
                print(" OK - QC for "+bnames[i][0])

                ind=np.where(wdic['hs']>0.1)

                if np.size(ind)>0:

                    ind=ind[0]

                    at=wdic['time'][ind]; adate=wdic['date'][ind]
                    ahs=wdic['hs'][ind]

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

                    alat=wdic['latitude'][0]; alon=wdic['longitude'][0]

                    indlat = np.where( np.abs(latc-alat) == np.nanmin(np.abs(latc-alat)) )[0][0]
                    indlon = np.where( np.abs(lonc-alon) == np.nanmin(np.abs(lonc-alon)) )[0][0]

                    initime = str(pd.to_datetime(adate.min()).year)+str(pd.to_datetime(adate.min()).month).zfill(2)+str(pd.to_datetime(adate.min()).day).zfill(2)+str(pd.to_datetime(adate.min()).hour).zfill(2)
                    fintime = str(pd.to_datetime(adate.max()).year)+str(pd.to_datetime(adate.max()).month).zfill(2)+str(pd.to_datetime(adate.max()).day).zfill(2)+str(pd.to_datetime(adate.max()).hour).zfill(2)
                    aftime = np.array(np.arange(float(timegm( time.strptime(initime, '%Y%m%d%H') )),float(timegm( time.strptime(fintime, '%Y%m%d%H') ))+1,wdt)).astype('double')

                    bhs=[]; bwnd=[]; btp=[]; btm=[]
                    btime=np.double([]); brt=np.double([])
                    abcmap=[]; abcid=[]; abcsec=[]
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

                            del indt

                            # check cyclone presence
                            indc=np.where(np.abs(ctime-at[t])<=5400.)
                            if np.size(indc)>0:
                                if cmap[np.min(indc[0]),indlat,indlon]>0:
                                    abcmap = np.append(abcmap,int(cmap[np.min(indc[0]),indlat,indlon]))
                                    abcid = np.append(abcid,int(cid[np.min(indc[0]),indlat,indlon]))
                                    abcsec = np.append(abcsec,int(csec[np.min(indc[0]),indlat,indlon]))
                                else:
                                    abcmap = np.append(abcmap,0)
                                    abcid = np.append(abcid,0)
                                    abcsec = np.append(abcsec,0)
                            else:
                                abcmap = np.append(abcmap,0)
                                abcid = np.append(abcid,0)
                                abcsec = np.append(abcsec,0)

                            del indc


                        print(" Ok data allocation "+repr(t)+" "+bnames[i][0])

                    del ind, indlat, indlon
                    # insert the position in the model grid 
                    indlat = np.where( np.abs(latm-alat) == np.nanmin(np.abs(latm-alat)) )[0][0]
                    indlon = np.where( np.abs(lonm-alon) == np.nanmin(np.abs(lonm-alon)) )[0][0]

                    # Final arrays
                    btime=np.array(btime).astype('double'); brt=np.array(brt).astype('double')
                    blat=np.zeros((brt.shape[0]),'f')+alat; blon=np.zeros((brt.shape[0]),'f')+alon
                    # bglat=np.zeros((brt.shape[0]),'f')+latm[indlat]; bglon=np.zeros((brt.shape[0]),'f')+lonm[indlon]
                    agidlat=np.zeros((brt.shape[0]),'int')+int(indlat); agidlon=np.zeros((brt.shape[0]),'int')+int(indlon)
                    abid=np.array(np.zeros((brt.shape[0]),'i')).astype('str')
                    if buoyd=="NDBC":
                        abid[:]="NDBC"+fname.split('/')[-1].split('h')[0]
                    else:
                        abid[:]="CDIP"+fname.split('/')[-1].split('_')[2]

                    ftime=np.append(ftime,btime); frtime=np.append(frtime,brt)
                    lat=np.append(lat,blat); lon=np.append(lon,blon); bid=np.append(bid,abid)
                    # glat=np.append(glat,bglat); glon=np.append(glon,bglon)
                    gidlat=np.append(gidlat,agidlat); gidlon=np.append(gidlon,agidlon)
                    hs=np.append(hs,bhs); tp=np.append(tp,btp); tm=np.append(tm,btm); wnd=np.append(wnd,bwnd)
                    bcmap = np.append(bcmap,abcmap)
                    abcid[abcid<0]=0.; abcid[np.isnan(abcid)==True]=0.; bcid = np.append(bcid,abcid)
                    abcsec[abcsec<0]=0.; abcsec[np.isnan(abcsec)==True]=0.; bcsec = np.append(bcsec,abcsec)

                    print(bnames[i][0]+" done")

                    del indlat,indlon,bhs,btp,btm,bwnd,btime,brt,aftime,blat,blon,abcmap,abcid #,bglat,bglon

            del wdic

        del fname


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

        if buoyd=="NDBC":
            df.to_csv('Data_NDBC.txt', sep='\t', index=False, header=True)
        else:
            df.to_csv('Data_CDIP.txt', sep='\t', index=False, header=True)


