#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
 Tropical Cyclone Parametric Wave Fields.
 Example with 1000-member ensemble from the Monte Carlo Simulation.
 One forecast cycle.
 See TCpwaves.py
"""

import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
import pandas as pd
from matplotlib.dates import DateFormatter
from scipy.ndimage import gaussian_filter
import timeit
import TCpwaves
from TCpwaves import *

if __name__ == "__main__":

    skg=2

    # PModel configuration file
    wdname=str('/home/ricardo/work/noaa/github/TC_waves/pmodel/pmodel_config.yaml')

    # Land/sea mask, water depth and distance to the coast.
    f=nc.Dataset('gridInfo.nc')
    latp=f.variables['latitude'][:]; mres=np.diff(latp).mean()
    lonp=f.variables['longitude'][:]; # lonp[lonp>180]=lonp[lonp>180]-360.
    mask=f.variables['mask'][:]

    # Read mcrall_al092022_092400.feather file
    fdate = "2022092400"
    fpath="/home/ricardo/work/noaa/analysis/TC_Waves/modeling/PModel/data"
    tcs=read_mcrall(fdate,fpath)
    # Plot tracks
    ctrack(tcs['lat'],tcs['lon']-360.,tcs['dtime'][:],tcs['Vmax'],flabel="MCsim",ftitle="no")
    ctrack(tcs['lat'],tcs['lon']-360.,tcs['dtime'][:],tcs['Vmax'],flabel="MCsim2",ftitle="no",style=2)

    # Select only the area of interest affected by the TCs
    indlatmin=np.min(np.where( np.abs(latp-(np.nanmin(tcs['lat'])-3.))==np.min(np.abs(latp-(np.nanmin(tcs['lat'])-3.))) ))
    indlatmax=np.min(np.where( np.abs(latp-(np.nanmax(tcs['lat'])+3.))==np.min(np.abs(latp-(np.nanmax(tcs['lat'])+3.))) ))
    indlonmin=np.min(np.where( np.abs(lonp-(np.nanmin(tcs['lon'])-3.))==np.min(np.abs(lonp-(np.nanmin(tcs['lon'])-3.))) ))
    indlonmax=np.min(np.where( np.abs(lonp-(np.nanmax(tcs['lon'])+3.))==np.min(np.abs(lonp-(np.nanmax(tcs['lon'])+3.))) ))

    latp = latp[indlatmin:indlatmax+1][::skg]; lonp = lonp[indlonmin:indlonmax+1][::skg]
    mask = mask[indlatmin:indlatmax+1,:][:,indlonmin:indlonmax+1][::skg,:][:,::skg]
    mask[mask>=0]=1.
    # mask[mask!=1]=np.nan

    WFHS = np.zeros((len(tcs['idx']),len(tcs['time']),len(latp),len(lonp)),'f')*np.nan
    W={'lat':latp,'lon':lonp,'Hs':WFHS[0,0,:,:]}
    hmax = np.zeros((len(tcs['idx']),len(tcs['time'])),'f')*np.nan

    # start = timeit.default_timer()
    for i in range(0,len(tcs['idx'])):
        start = timeit.default_timer()
        rangle = cbearing(tcs['lat'][i,:],tcs['lon'][i,:])

        for j in range(0,len(tcs['time'])):
            tcw = TCWaves(Vmax=tcs['Vmax'][i,j],Vfm=tcs['Vfm'][i,j],Rmax=tcs['Rmax'][i,j],R34=tcs['R34'][i,j],Lat=tcs['lat'][i,j],Lon=tcs['lon'][i,j],wdname=wdname)
            pwm = tcw.PWModel()
            pwm = tcw.xy_to_lonlat(pwm)
            pwm = tcw.rotate(pwm,rangle[j],wvar="Hs")
            pwm,WW = tcw.blend(pwm,W,wvar="Hs")
            WFHS[i,j,:,:]=np.array(WW['Hs'][:,:])
            hmax[i,j] = float(np.nanmax(pwm['Hs']))
            del pwm,WW,tcw
            print(" +++++ j "+repr(j))

        print(" ========== i "+repr(i))
        wplot(latp,lonp,np.nanmax(WFHS[i,:,:,:],axis=0)*mask,wvar="Hs",flabel=str(tcs['idx'][i]).zfill(3),ftitle=" - m"+str(tcs['idx'][i]).zfill(3))

        # 780s each member
        stop = timeit.default_timer()
        print("Member "+repr(i)+" in "+repr(int(round(stop - start,0)))+" seconds")

    # stop = timeit.default_timer()
    # print('Concluded in '+repr(int(round(stop - start,0)))+' seconds')

    # Probabilities
    qlev_hs = np.array([4.0, 6.0, 9.0, 14.0])
    ind = np.where( np.nanmean(hmax,axis=1)>0. )[0]
    print(" Total of "+repr(len(ind))+" members used")
    hmax = hmax[ind,:]
    WFHS = WFHS[ind,:,:,:]
    WFHS = np.nanmax(WFHS,axis=1)
    PWHS = np.zeros((len(qlev_hs),WFHS.shape[1],WFHS.shape[2]),'f')
    for i in range(0,len(qlev_hs)):
        for j in range(0,WFHS.shape[1]):
            for k in range(0,WFHS.shape[2]):
                aux = np.where(WFHS[:,j,k]>qlev_hs[i])
                if size(aux)>0:
                    PWHS[i,j,k] = float(len(aux[0])/int(WFHS.shape[0]))

    # Plot
    for i in range(0,WFHS.shape[0]):
        wplot(latp,lonp,WFHS[i,:,:]*mask,wvar="Hs",flabel=str(tcs['idx'][i]).zfill(3),ftitle=" - m"+str(tcs['idx'][i]).zfill(3))


    # Probabilities
    plevels=np.arange(10,101,10)
    for i in range(0,len(qlev_hs)):
        fig = plt.figure(figsize=(10, 6))
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.set_extent([np.min(lonp)-0.1, np.max(lonp)+0.1,np.min(latp)-0.1, np.max(latp)+0.1], crs=ccrs.PlateCarree())
        ax.coastlines(resolution='10m', linewidth=1)
        ax.add_feature(cfeature.BORDERS, linewidth=0.5)
        ax.add_feature(cfeature.LAND, facecolor='lightgray')
        ax.add_feature(cfeature.OCEAN, facecolor='lightblue')
        ax.gridlines(draw_labels=True, linestyle='--', color='gray', alpha=0.5)
        cs=ax.contourf(lonp,latp,gaussian_filter(PWHS[i,:,:]*100.,2)*mask,levels=plevels,alpha=0.7,cmap='gist_stern_r',zorder=1,extend="max",transform = ccrs.PlateCarree())
        title = "Prob Hs > "+str(int(qlev_hs[i]))+" m"
        ax.set_title(title); del title
        plt.tight_layout()
        ax = plt.gca(); pos = ax.get_position(); l, b, w, h = pos.bounds; cax = plt.axes([l+0.06, b-0.07, w-0.12, 0.03]) # setup colorbar axes.
        cbar = plt.colorbar(cs,cax=cax, orientation='horizontal', format='%g')
        plt.axes(ax); plt.tight_layout()
        figname = "PHs_above"+str(int(qlev_hs[i]))+"m"
        plt.savefig(figname+".png", dpi=200, facecolor='w', edgecolor='w',
            orientation='portrait', format='png',transparent=False, bbox_inches='tight', pad_inches=0.1)

        plt.close(fig)

