#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
 Tropical Cyclone Parametric Wave Fields.
 Example covering different combinations of Vmax, Vfm, Rmax, R34
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


    # Input parameters
    # Vmax=float(sys.argv[1])
    # Vfm=float(sys.argv[2])
    # Rmax=float(sys.argv[3])
    # R34=float(sys.argv[4])
    # Lat=float(sys.argv[5])
    # Lon=float(sys.argv[6])
    ## yaml config name
    # wdname=str(sys.argv[7])

    Vmax=float(50)
    Vfm=float(5)
    Rmax=float(30000)
    R34=float(300000)
    Lat=float(29.4)
    Lon=float(-77.3)
    rangle=float(0.)
    wdname=str('/home/ricardo/work/noaa/github/TC_waves/pmodel/pmodel_config.yaml')

    if (Vmax > 78 or Vfm > 15 or Rmax > 60e3 or R34 > 400e3 or Vmax < 17 or Rmax < 15e3 or R34 < 200e3):
        warnings.warn('Parameter(s) over the limits. Default(s) will be used instead.')

    # Land/sea mask, water depth and distance to the coast.
    f=nc.Dataset('gridInfo.nc')
    latp=f.variables['latitude'][:]; mres=np.diff(latp).mean()
    lonp=f.variables['longitude'][:]; # lonp[lonp>180]=lonp[lonp>180]-360.
    mask=f.variables['mask'][:]

    FHs = np.zeros(mask.shape,'f'); FHs[np.isnan(mask)==True]=np.nan
    FTp = np.zeros(mask.shape,'f'); FTp[np.isnan(mask)==True]=np.nan
    FUwnd = np.zeros(mask.shape,'f'); FUwnd[np.isnan(mask)==True]=np.nan
    FVwnd = np.zeros(mask.shape,'f'); FVwnd[np.isnan(mask)==True]=np.nan

    W={'lat':latp,'lon':lonp,'Hs':FHs, 'Tp':FTp, 'Uwnd':FUwnd, 'Vwnd':FVwnd}

    # Build object
    tcw = TCWaves(Vmax=Vmax,Vfm=Vfm,Rmax=Rmax,R34=R34,Lat=Lat,Lon=Lon,wdname=wdname)
    # Apply the parametric model. HS,TP,U,V,XX,YY
    pwm = tcw.PWModel()
    # Add lat lon, from the x and y arrays in km
    pwm = tcw.xy_to_lonlat(pwm)
    # Rotate
    pwm = tcw.rotate(pwm,rangle)
    # Blend cyclone into the large final domain, with 5km resolution
    pwm,W = tcw.blend(pwm,W)
    flabel="test";  ftitle="test"
    tcw.hwplot(pwm,W,flabel,ftitle)

    # Tests with diff combinations
    aVmax=np.array([17,30,40,50,65,78])
    aVfm=np.array([0,2.5,5.0,7.5,10.0,12.5,15.0])
    aRmax=np.array([15000,30000,60000])
    aR34=np.array([200000,300000,400000])

    hmax = np.zeros((len(aVmax),len(aVfm),len(aRmax),len(aR34)),'f')*np.nan

    for i in range(0,len(aVmax)):
        for j in range(0,len(aVfm)):
            for k in range(0,len(aRmax)):
                for l in range(0,len(aR34)):

                    Vmax=aVmax[i]; Vfm=aVfm[j]
                    Rmax=aRmax[k]; R34=aR34[l]

                    tcw = TCWaves(Vmax=Vmax,Vfm=Vfm,Rmax=Rmax,R34=R34,Lat=Lat,Lon=Lon,wdname=wdname)
                    pwm = tcw.PWModel()
                    pwm = tcw.xy_to_lonlat(pwm)
                    # pwm = tcw.rotate(pwm,rangle)
                    # pwm,WW = tcw.blend(pwm,W)

                    flabel="Vmax"+repr(Vmax)+"_Vfm10"+repr(int(Vfm*10))+"_Rmax"+repr(int(Rmax/1000))+"_R34."+repr(int(R34/1000))
                    tcw.hwplot(pwm,WW,flabel)

                    hmax[i,j,k,l] = float(np.nanmax(pwm['Hs']))

                    print(flabel)
                    del pwm,WW,tcw,Vmax,Vfm,Rmax,R34,flabel

    # Plot Hmax
    hmax.max() # 27.22

    # Visualize it
    plt.figure(figsize=(6,6))
    plt.imshow(hmax[:,2,:,:], origin='lower', cmap='viridis')
    plt.colorbar(label='Hmax')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.show()
    plt.savefig("Teste.png", dpi=200, facecolor='w', edgecolor='w',
        orientation='portrait', format='png',transparent=False, bbox_inches='tight', pad_inches=0.1)


