#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
 Tropical Cyclone Parametric Wave Fields.
 Example for Hurricane Lee
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

    # Example Hurricane Lee
    # https://www.ncei.noaa.gov/sites/default/files/2021-07/IBTrACS_v04_column_documentation.pdf
    dfibtr = pd.read_csv('ibtracs_Hlee.csv', header=[0,1])
    atime = pd.to_datetime(dfibtr['ISO_TIME'].values[:, 0])
    alat = np.array(np.array(dfibtr['LAT'].values[:,0]).astype('float'))
    alon = np.array(np.array(dfibtr['LON'].values[:,0]).astype('float'))
    aVmax = np.array(np.array(dfibtr['USA_WIND'].values[:,0]).astype('float')) # wind in knots
    aVfm = np.array(dfibtr['STORM_SPEED'].values[:,0])/1.94384
    aRmax = np.array(dfibtr['USA_RMW'].values[:,0])*1852.
    aR34 = np.mean(np.array([dfibtr['USA_R34_NE'].values[:,0],dfibtr['USA_R34_SE'].values[:,0],
        dfibtr['USA_R34_SW'].values[:,0],dfibtr['USA_R34_NW'].values[:,0]]).astype('float'),axis=0)*1.852*1000.

    # plot cyclone track
    ctrack(alat,alon,atime,aVmax,flabel="HLee",ftitle="Hurricane Lee")

    # PModel uses m/s
    aVmax=aVmax/1.94384 # wind in m/s

    # heading directions
    rangle = cbearing(alat,alon)

    hmax = np.zeros((len(atime)),'f')*np.nan
    for i in range(0,len(atime)):
        tcw = TCWaves(Vmax=aVmax[i],Vfm=aVfm[i],Rmax=aRmax[i],R34=aR34[i],Lat=alat[i],Lon=alon[i],wdname=wdname)
        pwm = tcw.PWModel()
        pwm = tcw.xy_to_lonlat(pwm)
        pwm = tcw.rotate(pwm,rangle[i])
        pwm,WW = tcw.blend(pwm,W)

        flabel="HLee_"+atime[i].strftime('%Y%m%d%H')
        ftitle=" "+atime[i].strftime('%Y%m%d')+" "+atime[i].strftime('%H')+"Z"
        tcw.hwplot(pwm,WW,flabel,ftitle)

        hmax[i] = float(np.nanmax(pwm['Hs']))

        print(flabel)
        del pwm,WW,tcw,flabel,ftitle


    fig1, ax = plt.subplots(figsize=(9, 4))
    ax.plot(atime, hmax, color='dimgray', linewidth=2., zorder=2)
    ax.plot(atime, hmax, color='k', marker='.', linestyle='', linewidth=2., zorder=2)
    ax.set_xlim(atime[0], atime[-1])
    ax.xaxis.set_major_formatter(DateFormatter('%b%d'))
    ax.fmt_xdata = DateFormatter('%b%d')
    ax.set_xlabel("Time (2023)")
    ax.set_ylabel("Hs (m)")
    plt.tight_layout()
    plt.grid(c='grey', ls='--', alpha=0.3, zorder=1)
    plt.savefig('TimeSeries_HsMax_Hlee.png', dpi=200, facecolor='w', edgecolor='w', orientation='portrait',
        format='png', bbox_inches='tight', pad_inches=0.1)


    fig, ax1 = plt.subplots(figsize=(13, 5))
    ax1.plot(atime, hmax, color='dimgray', linewidth=2., label='Hs', zorder=3)
    ax1.plot(atime, hmax, color='dimgray', marker='.', linewidth=2., zorder=3)
    ax1.set_ylabel("Hs (m)", color='dimgray')
    ax1.tick_params(axis='y', labelcolor='dimgray')
    ax1.set_xlim(atime[0], atime[-1])
    ax1.xaxis.set_major_formatter(DateFormatter('%b%d'))
    ax1.fmt_xdata = DateFormatter('%b%d')
    ax1.set_xlabel("Time (2023)")
    ax1.grid(c='grey', ls='--', alpha=0.3, zorder=1)

    ax2 = ax1.twinx()
    ax2.plot(atime, aVmax, color='royalblue', linewidth=2., zorder=2)
    ax2.set_ylabel("Wsp max (m/s)", color='royalblue')
    ax2.tick_params(axis='y', labelcolor='royalblue')

    ax3 = ax1.twinx()
    ax3.spines["left"].set_position(("axes", -0.09))
    ax3.plot(atime, aVfm, color='firebrick', linewidth=2., )
    ax3.set_ylabel("Vfm (m/s)", color='firebrick')
    ax3.tick_params(axis='y', labelcolor='firebrick')
    ax3.yaxis.set_label_position("left"); ax3.yaxis.set_ticks_position("left")

    ax4 = ax1.twinx()
    ax4.spines["right"].set_position(("axes", 1.1))
    ax4.plot(atime, aR34/1000., color='green', linewidth=2., zorder=1)
    ax4.set_ylabel("R34 (km)", color='green')
    ax4.tick_params(axis='y', labelcolor='green')

    ax5 = ax1.twinx()
    ax5.spines["right"].set_position(("axes", 1.22))
    ax5.plot(atime, aRmax/1000., color='orange', linewidth=2.)
    ax5.set_ylabel("Rmax (km)", color='orange')
    ax5.tick_params(axis='y', labelcolor='orange')

    plt.title("Hurricane Lee",fontsize=13)
    fig.tight_layout()
    plt.savefig('TimeSeries_Hlee.png', dpi=200, facecolor='w', edgecolor='w', orientation='portrait',
        format='png', bbox_inches='tight', pad_inches=0.1)

    # convert -delay 15 -loop 0 wspectrum_*.png wspectrum.gif

