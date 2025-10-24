#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
 Tests with syntheric data.
 Mail goal is to check the spatial distribution and rotation of wave fields inside TCs.
 See TCpwaves.py
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
import pandas as pd
from matplotlib.dates import DateFormatter
import TCpwaves
from TCpwaves import *

if __name__ == "__main__":

    wdname=str('/home/ricardo/work/noaa/github/TC_waves/pmodel/pmodel_config.yaml')

    # Large fixed domain
    latp=np.arange(2.,67.1,0.1)
    # latp=np.arange(15.,50.,0.1)
    lonp=np.arange(270.,340.1,0.1)
    # lonp=np.arange(270.,355.1,0.1)
    FHs = np.zeros((len(latp),len(lonp)),'f')+1.
    FTp = np.zeros((len(latp),len(lonp)),'f')+1.
    FUwnd = np.zeros((len(latp),len(lonp)),'f')+1.
    FVwnd = np.zeros((len(latp),len(lonp)),'f')+1.
    W={'lat':latp,'lon':lonp,'Hs':FHs, 'Tp':FTp, 'Uwnd':FUwnd, 'Vwnd':FVwnd}

    # Just taking the fixed arrays from H Lee.
    dfibtr = pd.read_csv('ibtracs_Hlee.csv', header=[0,1])
    atime = pd.to_datetime(dfibtr['ISO_TIME'].values[:, 0])
    aRmax = np.array(dfibtr['USA_RMW'].values[:,0])*1852.
    aR34 = np.mean(np.array([dfibtr['USA_R34_NE'].values[:,0],dfibtr['USA_R34_SE'].values[:,0],
        dfibtr['USA_R34_SW'].values[:,0],dfibtr['USA_R34_NW'].values[:,0]]).astype('float'),axis=0)*1.852*1000.

    atime = atime[::4][0:-4]; aRmax = aRmax[::4][0:-4]; aR34 = aR34[::4][0:-4]

    # Synthetic Cyclone
    vfm_sim = np.array([2.,5.,10.])
    vmax_sim = np.array([70., 100., 150.])

    hmax = np.zeros((len(vfm_sim),len(vmax_sim),len(atime)),'f')*np.nan

    alon = np.full(len(atime), -35.0)
    # alat = np.full(len(atime), 33.0)

    for i in range(0,len(vfm_sim)):

        l_step = vfm_sim[i] * (1.852) * (1. / 111.12) * np.mean(np.diff(atime) / np.timedelta64(1, 'h'))
        alat = 30. -vfm_sim[i] + np.arange(len(atime)) * l_step
        # alon = 345. - np.arange(len(atime)) * l_step
        # alon = 285. + np.arange(len(atime)) * l_step

        aVfm = np.full(len(atime), vfm_sim[i])

        # heading directions
        rangle = cbearing(alat,alon)

        for j in range(0,len(vmax_sim)):
            aVmax = np.full(len(atime), vmax_sim[j])/1.94384
            for t in range(0,len(atime)):
                tcw = TCWaves(Vmax=aVmax[t],Vfm=aVfm[t],Rmax=aRmax[t],R34=aR34[t],Lat=alat[t],Lon=alon[t],wdname=wdname)
                pwm = tcw.PWModel()
                pwm = tcw.xy_to_lonlat(pwm)
                pwm = tcw.rotate(pwm,rangle[t])
                pwm,WW = tcw.blend(pwm,W)

                flabel="Northwards_"+atime[t].strftime('%Y%m%d%H')+"_"+repr(int(vfm_sim[i]))+"_"+repr(int(vmax_sim[j]))
                # flabel="Westwards_"+atime[t].strftime('%Y%m%d%H')+"_"+repr(int(vfm_sim[i]))+"_"+repr(int(vmax_sim[j]))
                # flabel="Eastwards_"+atime[t].strftime('%Y%m%d%H')+"_"+repr(int(vfm_sim[i]))+"_"+repr(int(vmax_sim[j]))
                ftitle=flabel
                tcw.hwplot(pwm,WW,wvar="Hs",flabel=flabel,ftitle=ftitle)

                hmax[i,j,t] = float(np.nanmax(pwm['Hs']))

                print(flabel)
                del pwm,WW,tcw,flabel,ftitle

