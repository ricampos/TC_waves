#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib
matplotlib.use('Agg')
import math
from xarray import open_dataset
import netCDF4 as nc
import numpy as np
from pylab import *
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import numpy.ma as ma
import sys
import time
from matplotlib.dates import date2num
import matplotlib.ticker as ticker
import datetime
import pickle
from scipy import signal
from matplotlib import mlab
from scipy.ndimage.filters import gaussian_filter
from mpl_toolkits.basemap import Basemap, shiftgrid, cm
import mpl_toolkits.basemap
from matplotlib import ticker
colormap = cm.GMT_polar
palette = plt.cm.jet
palette.set_bad('aqua', 1.0)
from sklearn.utils import shuffle
import scipy.stats
import os
# matplotlib.use('Agg')

pgdas="/home/ricardo/work/noaa/analysis/TC_Waves/data/GDAS"

# Font size and style
sl=14
matplotlib.rcParams.update({'font.size': sl}); plt.rc('font', size=sl) 
matplotlib.rc('xtick', labelsize=sl); matplotlib.rc('ytick', labelsize=sl); matplotlib.rcParams.update({'font.size': sl})

# Read Cyclone Info
ds = xr.open_dataset('CycloneMap.nc')
ctime = ds.time.values; lat = ds['lat'].values; lon = ds['lon'].values
cmap = ds['cmap'].values; csec = ds['csec'].values; cid = ds['cid'].values
ds.close(); del ds

cnames = np.array(pd.read_csv('cyclone_names.txt', header=None).values).astype('str')
# csel = ['FIONA','IAN','IDALIA','LEE','FRANKLIN','HELENE','MILTON']
csel = csel = ['FIONA','IAN','IDALIA','FRANKLIN','HELENE','MILTON']

levels=np.array([-1,0,1,2])
slevels=np.array([0,1,2,3,4,5])
hlevels = np.linspace(0,11.2,113)
wlevels = np.linspace(0,40,41)

for t in range(0,len(ctime)):

    # select a few cyclones to plot
    goahead=0
    aux = np.unique(cid[t,:,:]); aux=aux[aux>=0]
    for i in range(0,len(aux)):
        if cnames[aux[i]] in csel:
            goahead=1

    if goahead==1:

        cptime=np.datetime_as_string(ctime[t], unit='h').replace('-', '').replace('T', '_')

        # read GDAS
        ds = xr.open_dataset(pgdas+"/gdaswave."+np.datetime_as_string(ctime[t], unit='D').replace('-', '')+".global.0p16.nc")
        mtime = ds.time.values
        indt=np.where(mtime==ctime[t])
        if np.size(indt)>0:
            indt=indt[0][0]
            hs = ds.HTSGW_surface.values[indt,:,:]; wnd=ds.WIND_surface.values[indt,:,:]
            # cyclone map, sectors, and mask
            fcmap = cmap[t,:,:]
            fcsec = csec[t,:,:]

            mask = np.zeros(hs.shape, dtype=bool)
            mask[fcmap>0]=True
            chs = ma.masked_where(~mask, hs)
            cwnd = ma.masked_where(~mask, wnd)

            [mnlon,mnlat]=np.meshgrid(lon,lat)

            # Plot Hs
            plt.figure(figsize=(10,10))
            map = Basemap(projection='ortho',lat_0=15,lon_0=-80,area_thresh=1.0,resolution='l')
            map.drawcoastlines(linewidth=0.25, zorder=3)
            map.drawcountries(linewidth=0.25, zorder=3)
            map.fillcontinents(color='grey', zorder=3)
            map.drawlsmask(ocean_color='w',land_color='grey', zorder=1)   
            xc,yc = map(mnlon,mnlat)  
            map.contourf(xc,yc,hs,hlevels,extend="max", vmin=0., cmap=plt.cm.jet,alpha=0.3, antialiased=True, zorder=2)
            map.contour(xc,yc,fcmap,levels,vmin=0.5, vmax=1.5, zorder=3,colors='k', linewidths=1)
            cs=map.contourf(xc,yc,chs,hlevels,extend="max", vmin=0., cmap=plt.cm.jet, zorder=3)
            map.fillcontinents(color='silver', zorder=3)
            map.drawcoastlines(linewidth=0.8)
            map.drawcountries(linewidth=0.5, linestyle='solid', color='k', antialiased=1, ax=None, zorder=4)
            map.drawmeridians(np.arange(0,360,30), zorder=5)
            map.drawparallels(np.arange(-90,90,30), zorder=5)
            plt.title(cptime[0:8]+" "+str(int(cptime[9:11])).zfill(2)+"Z", fontsize=18)
            cbar = plt.colorbar(cs, orientation='horizontal', pad=0.05, fraction=0.02, shrink=0.5, aspect=40)
            plt.tight_layout()
            savefig("CycloneMap_Hs_"+cptime+".png", dpi=150, facecolor='w', edgecolor='w',orientation='portrait', format='png',transparent=False, bbox_inches='tight', pad_inches=0.1)
            plt.close('all'); del map, cs

            # Plot Wnd
            plt.figure(figsize=(10,10))
            map = Basemap(projection='ortho',lat_0=15,lon_0=-80,area_thresh=1.0,resolution='l')
            map.drawcoastlines(linewidth=0.25, zorder=3)
            map.drawcountries(linewidth=0.25, zorder=3)
            map.fillcontinents(color='grey', zorder=3)
            map.drawlsmask(ocean_color='w',land_color='grey', zorder=1) 
            xc,yc = map(mnlon,mnlat)    
            map.contourf(xc,yc,wnd,wlevels,extend="max", vmin=0., cmap=plt.cm.jet,alpha=0.3, antialiased=True, zorder=2)
            map.contour(xc,yc,fcmap,levels,vmin=0.5, vmax=1.5, zorder=3,colors='k', linewidths=1)
            cs=map.contourf(xc,yc,cwnd,wlevels,extend="max", vmin=0., cmap=plt.cm.jet, zorder=3)
            map.fillcontinents(color='silver', zorder=3)
            map.drawcoastlines(linewidth=0.8)
            map.drawcountries(linewidth=0.5, linestyle='solid', color='k', antialiased=1, ax=None, zorder=4)
            map.drawmeridians(np.arange(0,360,30), zorder=5)
            map.drawparallels(np.arange(-90,90,30), zorder=5)
            plt.title(cptime[0:8]+" "+str(int(cptime[9:11])).zfill(2)+"Z", fontsize=18)
            cbar = plt.colorbar(cs, orientation='horizontal', pad=0.05, fraction=0.02, shrink=0.5, aspect=40)
            plt.tight_layout()
            savefig("CycloneMap_Wnd_"+cptime+".png", dpi=150, facecolor='w', edgecolor='w',orientation='portrait', format='png',transparent=False, bbox_inches='tight', pad_inches=0.1)
            plt.close('all')

            # Plot Sectors
            plt.figure(figsize=(10,10))
            map = Basemap(projection='ortho',lat_0=15,lon_0=-80,area_thresh=1.0,resolution='l')
            map.drawcoastlines(linewidth=0.25, zorder=3)
            map.drawcountries(linewidth=0.25, zorder=3)
            map.fillcontinents(color='grey', zorder=3)
            map.drawlsmask(ocean_color='w',land_color='grey', zorder=1)     
            map.contourf(xc,yc,fcsec,slevels,extend="max", vmin=0., cmap=plt.cm.jet, zorder=2)
            map.contour(xc,yc,fcmap,levels,vmin=0.5, vmax=1.5, zorder=3,colors='k', linewidths=1)
            map.fillcontinents(color='silver', zorder=3)
            map.drawcoastlines(linewidth=0.8)
            map.drawcountries(linewidth=0.5, linestyle='solid', color='k', antialiased=1, ax=None, zorder=4)
            map.drawmeridians(np.arange(0,360,30), zorder=5)
            map.drawparallels(np.arange(-90,90,30), zorder=5)
            plt.title(cptime[0:8]+" "+str(int(cptime[9:11])).zfill(2)+"Z", fontsize=18)
            plt.tight_layout()
            savefig("CycloneMap_Sectors_"+cptime+".png", dpi=150, facecolor='w', edgecolor='w',orientation='portrait', format='png',transparent=False, bbox_inches='tight', pad_inches=0.1)
            plt.close('all')

            print("Ok "+cptime)

        else:
            print("GDAS problem "+cptime)


# All Afected Areas during the 3 years
cmap[cmap<5]=0; cmap[cmap==5]=1 # Cyclone only.
mcmap = np.array(np.sum(cmap,axis=0)).astype('float')
mcmap[mcmap==0.]=np.nan
# Somatory of number of hours affected in 3 years.
flevels=np.array(np.linspace(0,510,52)).astype('int')
plt.figure(figsize=(10,10))
map = Basemap(projection='ortho',lat_0=15,lon_0=-80,area_thresh=1.0,resolution='l')
map.drawcoastlines(linewidth=0.25, zorder=3)
map.drawcountries(linewidth=0.25, zorder=3)
map.fillcontinents(color='grey', zorder=3)
map.drawlsmask(ocean_color='w',land_color='grey', zorder=1) 
[mnlon,mnlat]=np.meshgrid(lon,lat)
xc,yc = map(mnlon,mnlat)    
cs=map.contourf(xc,yc,mcmap,flevels, vmin=0., cmap=plt.cm.jet, zorder=2)
map.fillcontinents(color='silver', zorder=3)
map.drawcoastlines(linewidth=0.8)
map.drawcountries(linewidth=0.5, linestyle='solid', color='k', antialiased=1, ax=None, zorder=4)
map.drawmeridians(np.arange(0,360,30), zorder=5)
map.drawparallels(np.arange(-90,90,30), zorder=5)
cbar = plt.colorbar(cs, orientation='horizontal', pad=0.05, fraction=0.02, shrink=0.5, aspect=40)
plt.tight_layout()
savefig("CycloneMap_AreasAffected.png", dpi=150, facecolor='w', edgecolor='w',orientation='portrait', format='png',transparent=False, bbox_inches='tight', pad_inches=0.1)
plt.close('all')

