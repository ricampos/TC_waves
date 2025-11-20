#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
TCpwaves.py

VERSION AND LAST UPDATE:
 v1.0  05/07/2025
 v1.1  05/20/2025

PURPOSE:
 Tropical Cyclone Parametric Wave Fields.
 The algorith generates wave field (Hs and Tp) based on inputs of
 Vmax: Maximum sustained wind speed (m/s)
 Vfm: Velocity of forward movement of the cyclone (m/s)
 Rmax: Radius of maximum winds (m)
 R34: Radius of gale (34 kt) winds (m)
 The position Lat/Lon must also be given.

USAGE:
 The configuration file pmodel_config.yaml has the coefficients, parameters, and also provides
  the paths
 The wave diagrams, PModel_WaveDiagrams.nc, is required to run this code (path must be informed 
  in the pmodel_config.yaml file)
 With these two files, one can enter the required information of the 
  cyclone (Vmax, Vfm, Rmax, R34, Lat, Lon, time), for example from IBTrACS, and generate multiple cyclones
 It is also possible to generate a fixed domain (a large area) onto which the cyclones can 
  be interpolated and blended. 

 Auxiliary functions:
  read
  read_mcrall
  read_atcf
  bearing
  cbearing
  ctrack
  wplot
 Class TCWaves
 Functions:
  PWModel
  xy_to_lonlat
  rotate
  blend
  hwplot

 The explanation for each function is contained in the headers, including examples, help(TCpwaves)

 The function PWModel below, regarding the parametric wave model, follows the paper:
  Grossmann-Matheson, G., Young, I.R., Meucci, A., Alves,J.-H., Tamizi, A., 2025.
  A model for the spatial distribution of ocean wave parameters in tropical cyclones. Ocean Engineering.
  https://doi.org/10.1016/j.oceaneng.2024.120091
  The function PWModel is based on the Matlab code given by
  https://doi.org/10.26188/27237156 , and great support provided by Guisela Grossmann-Matheson

OUTPUT:
 Wave information (Hs and Tp) for the given cyclones.

DEPENDENCIES:
 See he imports below.

AUTHOR and DATE:
 05/07/2025: Ricardo M. Campos, first version.
 05/20/2025: Ricardo M. Campos, new auxiliar functions to read NHC tracks and TC data

PERSON OF CONTACT:
 Ricardo M Campos: ricardo.campos@noaa.gov

"""

import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import yaml
from geopy.distance import distance
from geopy.point import Point
from scipy.interpolate import griddata
from scipy.interpolate import RegularGridInterpolator
from datetime import datetime
import pandas as pd
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.ticker
from mpl_toolkits.basemap import cm
colormap = cm.GMT_polar
palette = plt.cm.jet
palette.set_bad('aqua', 1.0)
import warnings; warnings.filterwarnings("ignore")

# Font size and style
sl=13
matplotlib.rcParams.update({'font.size': sl}); plt.rc('font', size=sl) 
matplotlib.rc('xtick', labelsize=sl); matplotlib.rc('ytick', labelsize=sl); matplotlib.rcParams.update({'font.size': sl})


def read(wdname):
    '''
    Read configuration file (.yaml format) and wave diagrams (netcdf format).
    This is used inside PWModel
    '''
    # fname="/home/ricardo/work/noaa/github/TC_waves/pmodel/pmodel_config.yaml"
    with open(wdname, 'r') as file:
        wconfig = yaml.safe_load(file)

    ds = xr.open_dataset(wconfig['pmwdiagrams_path'])

    return wconfig, ds
    del wconfig


def read_mcrall(fdate,fpath):
    '''
    Read Monte Carlo wind speed probability realization files mcrall.
    '''
    # fdate = "2022092400"
    # Time
    if isinstance(fdate, str):
        dt = datetime.strptime(fdate, "%Y%m%d%H")
        t = datetime.strptime(fdate, "%Y%m%d%H").timestamp()
    elif isinstance(fdate, float):
        t = float(fdate)
        dt = datetime.fromtimestamp(t)
    elif isinstance(fdate, datetime):
        dt = fdate
        t = datetime.strptime(fdate, "%Y%m%d%H").timestamp()

    fname = "mcrall_al"+dt.strftime("%m")+dt.strftime("%Y")+"_"+dt.strftime("%m")+dt.strftime("%d")+dt.strftime("%H")+".feather"
    # read data
    tcs = pd.read_feather(fpath+"/"+fname) 
    idx = np.array(tcs['idx'].values)
    flt = np.array(tcs['tau'].values) # forecast lead time
    aVmax = np.array(tcs['vmx'].values)/1.94384
    aRmax = np.array(tcs['rmw'].values)*1852.
    aR34 = np.array(tcs['r34q1'].values+tcs['r34q2'].values+tcs['r34q3'].values+tcs['r34q4'].values)/4.
    aVfm = np.array(tcs['spd'].values)/1.94384
    alat = np.array(tcs['lat'].values)
    alon = np.array(tcs['lon'].values)
    del tcs

    alon[alon<0.]=alon[alon<0.]+360.

    # Reshape and organize
    idx=np.unique(idx)
    flt=np.unique(flt)
    ftime = np.array(flt*3600.+t).astype('double')
    dtime = [datetime.fromtimestamp(t) for t in ftime]
    Vmax=np.reshape(aVmax,(len(idx),len(flt)))
    Rmax=np.reshape(aRmax,(len(idx),len(flt)))
    R34=np.reshape(aR34,(len(idx),len(flt)))
    Vfm=np.reshape(aVfm,(len(idx),len(flt)))
    lat=np.reshape(alat,(len(idx),len(flt)))
    lon=np.reshape(alon,(len(idx),len(flt)))

    result={'idx':idx,'flt':flt,'time':ftime,'dtime':dtime,
    'Vmax':Vmax,'Rmax':Rmax,'R34':R34,'Vfm':Vfm,
    'lat':lat, 'lon':lon}

    return result


def read_atcf(fname):
    '''
    Read Automated Tropical Cyclone Forecast (ATCF) Archive Data Files / Text Files
    https://web.uwm.edu/hurricane-models/models/models.html
    https://ftp.nhc.noaa.gov/atcf/archive/2024/
    https://ftp.nhc.noaa.gov/atcf/archive/README
      wind speed in knots
      pressure in hPa
    '''
    # GFSO: GFS
    # GEFS control: AC00
    # GEFS perturbed members: AP01 02 ..

    models = {}

    with open(fname, 'r') as f:
        lines = f.readlines()
    
    for line in lines:
        parts = [p.strip() for p in line.split(',')]
        if len(parts) < 10:
            continue

        basin = parts[0]
        storm_num = parts[1]
        init_time_str = parts[2]
        model = parts[4]
        fhr = int(parts[5])
        lat_str = parts[6]
        lon_str = parts[7]
        wind = int(parts[8])
        pres = int(parts[9]) if parts[9].isdigit() else None

        # Convert lat/lon from N/S/E/W to float
        lat = float(lat_str[:-1]) / 10.0 * (-1 if lat_str[-1] == 'S' else 1)
        lon = float(lon_str[:-1]) / 10.0 * (-1 if lon_str[-1] == 'W' else 1)

        # Convert init and forecast time
        init_time = pd.to_datetime(init_time_str, format='%Y%m%d%H')
        valid_time = init_time + pd.Timedelta(hours=fhr)

        # Initialize model dict
        if model not in models:
            models[model] = []

        # Append forecast info
        models[model].append((init_time, fhr, valid_time, lat, lon, wind, pres))

    # Convert to DataFrames
    for model in models:
        models[model] = pd.DataFrame(models[model], columns=['init_time', 'fhr', 'time', 'lat', 'lon', 'wind', 'pres'])

    return models
    # ind = np.where(model['AC00']['init_time']=='2024-09-25T06:00:00.000000000')[0]
    # ctrack(np.array(model['AC00']['lat'][ind]),np.array(model['AC00']['lon'][ind]),model['AC00']['time'][ind],np.array(model['AC00']['wind'][ind]),flabel="test_AC",ftitle="test_AC")


def read_ibtracs(fname):
    '''
    Read IBTracks text data
    ibtracs.last3years.list.v04r01.csv
    Only named storms are selected
    IBTRACKS V4 data. It can be downloaded at
    https://www.ncei.noaa.gov/data/international-best-track-archive-for-climate-stewardship-ibtracs/v04r00/access/csv/
    https://www.ncdc.noaa.gov/ibtracs/index.php?name=bib
    '''
    dfibtr = pd.read_csv(fname, header=[0,1])
    inat = np.array(dfibtr.values[:,7]); iname=np.array(dfibtr['NAME'].values[:][:,0]).astype('str')
    icyid = np.array(dfibtr['NUMBER'].values[:][:,0]).astype('int')
    ity = np.array(pd.Series(dfibtr['ISO_TIME'].values.flatten()).str.slice(0,4)).astype('int')
    itm = np.array(pd.Series(dfibtr['ISO_TIME'].values.flatten()).str.slice(5,7)).astype('int')
    itd = np.array(pd.Series(dfibtr['ISO_TIME'].values.flatten()).str.slice(8,10)).astype('int')
    ith = np.array(pd.Series(dfibtr['ISO_TIME'].values.flatten()).str.slice(11,13)).astype('int')
    itime = pd.to_datetime(dfibtr['ISO_TIME'].values[:,0])
    ilat = np.array(dfibtr['LAT'].values[:,0]).astype('float')
    ilon = np.array(dfibtr['LON'].values[:,0]).astype('float')
    iVmax = np.array(dfibtr['USA_WIND'].replace(' ', np.nan).values[:,0]).astype('float')/1.94384 # max wind speed in m/s
    iVfm = np.array(dfibtr['STORM_SPEED'].replace(' ', np.nan).values[:,0]).astype('float')/1.94384 # forward speed in m/s
    iRmax = np.array(dfibtr['USA_RMW'].replace(' ', np.nan).values[:,0]).astype('float')*1852. # in meters
    iR34 = np.array(np.mean(np.array([dfibtr['USA_R34_NE'].replace(' ', np.nan).values[:,0],dfibtr['USA_R34_SE'].replace(' ', np.nan).values[:,0],
        dfibtr['USA_R34_SW'].replace(' ', np.nan).values[:,0],dfibtr['USA_R34_NW'].replace(' ', np.nan).values[:,0]]).astype('float'),axis=0)*1852.)

    print('read_ibtracs: Ibtracks ok'); del dfibtr
    ind=np.where(iname!="UNNAMED")
    ity=np.copy(ity[ind[0]]);itm=np.copy(itm[ind[0]]); itd=np.copy(itd[ind[0]]); ith=np.copy(ith[ind[0]]); itime=np.copy(itime[ind[0]])
    icyid=np.copy(icyid[ind[0]]); iname=np.copy(iname[ind[0]]); ilat=np.copy(ilat[ind[0]]);ilon=np.copy(ilon[ind[0]]);inat=np.copy(inat[ind[0]])
    iVmax=np.copy(iVmax[ind[0]]); iVfm=np.copy(iVfm[ind[0]]); iRmax=np.copy(iRmax[ind[0]]); iR34=np.copy(iR34[ind[0]])
    del ind
    # ilon[ilon<0]=ilon[ilon<0]+360.

    icyid = np.char.add(ity.astype(int).astype(str), np.char.zfill(icyid.astype(int).astype(str), 3))

    wibtr={'ity':ity,'itm':itm,'itd':itd,'ith':ith,'itime':pd.to_datetime(itime),'icyid':np.array(icyid).astype('int'),'iname':iname,'ilat':ilat,'ilon':ilon,'inat':inat,
        'iVmax':iVmax,'iVfm':iVfm,'iRmax':iRmax,'iR34':iR34}

    return wibtr


def bearing(lat1, lon1, lat2, lon2):
    """
    Calculate the propagation direction (going to), based on two pairs of lat/lons
    """
    lat1_rad = np.deg2rad(lat1)
    lat2_rad = np.deg2rad(lat2)
    dlon_rad = np.deg2rad(lon2 - lon1)
    x = np.sin(dlon_rad) * np.cos(lat2_rad)
    y = np.cos(lat1_rad) * np.sin(lat2_rad) - np.sin(lat1_rad) * np.cos(lat2_rad) * np.cos(dlon_rad)
    bearing_rad = np.arctan2(x, y)
    bearing_deg = (np.rad2deg(bearing_rad) + 360) % 360  # normalize to [0, 360)
    return bearing_deg

def cbearing(alat,alon):
    """
    Generate the final propagation direction (going to) for given arrays of lat and lon
    """
    alon[alon>180]=alon[alon>180]-360.
    bdir1=np.zeros(len(alat),'f')*np.nan
    bdir2=np.zeros(len(alat),'f')*np.nan
    for i in range(0,len(alat)-1):
        bdir1[i] = bearing(alat[i],alon[i],alat[i+1],alon[i+1])

    for i in range(1,len(alat)):
        bdir2[i] = bearing(alat[i-1],alon[i-1],alat[i],alon[i])

    bdir1[-1]=bdir2[-1]; bdir2[0]=bdir1[0]

    rad1 = np.deg2rad(bdir1); rad2 = np.deg2rad(bdir2)
    sin_avg = (np.sin(rad1) + np.sin(rad2)) / 2
    cos_avg = (np.cos(rad1) + np.cos(rad2)) / 2
    avg_dir_rad = np.arctan2(sin_avg, cos_avg)
    avg_dir_deg = (np.rad2deg(avg_dir_rad) + 360) % 360
    return avg_dir_deg



# ---- Visualization Functions ----

def ctrack(alat,alon,atime,aVmax,flabel="no",ftitle="no",style=1):
    """
    Plot cyclone tracks and save the figure, based on arrays of lat, lon, time, and aVmax (in knots)
    """

    alat = np.atleast_2d(alat); alon = np.atleast_2d(alon)
    aVmax = np.atleast_2d(aVmax)

    alon[alon>180]=alon[alon>180]-360.

    # aVmax must be in knots
    fig = plt.figure(figsize=(10, 6))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent([np.min(alon)-8., np.max(alon)+8.,np.min(alat)-5., np.max(alat)+5.], crs=ccrs.PlateCarree())
    ax.coastlines(resolution='10m', linewidth=1)
    ax.add_feature(cfeature.BORDERS, linewidth=0.5)
    ax.add_feature(cfeature.LAND, facecolor='lightgray')
    ax.add_feature(cfeature.OCEAN, facecolor='lightblue')
    ax.gridlines(draw_labels=True, linestyle='--', color='gray', alpha=0.5)

    if style==1:

        for i in range(0,alat.shape[0]):

            ax.plot(alon[i,:], alat[i,:], color='white', linewidth=2,zorder=1)
            ax.plot(alon[i,:], alat[i,:], color='silver', linewidth=1,zorder=1)

            for j in range(0,alat.shape[1]):

                if atime[j].hour == 0:
                    # https://www.nhc.noaa.gov/aboutsshws.php
                    if aVmax[i,j]<34:
                        ax.plot(alon[i,j], alat[i,j], marker='o', color='navy', linewidth=1,zorder=2)
                    if aVmax[i,j]>=34 and aVmax[i,j]<48:
                        ax.plot(alon[i,j], alat[i,j], marker='o', color='darkgreen', linewidth=1,zorder=2)
                    if aVmax[i,j]>=48 and aVmax[i,j]<64:
                        ax.plot(alon[i,j], alat[i,j], marker='o', color='gold', linewidth=1,zorder=2)
                    if aVmax[i,j]>=64 and aVmax[i,j]<83:
                        ax.plot(alon[i,j], alat[i,j], marker='o', color='orange', linewidth=1,zorder=4)
                    if aVmax[i,j]>=83 and aVmax[i,j]<96:
                        ax.plot(alon[i,j], alat[i,j], marker='o', color='darkorange', linewidth=1,zorder=4)
                    if aVmax[i,j]>=96 and aVmax[i,j]<113:
                        ax.plot(alon[i,j], alat[i,j], marker='o', color='red', linewidth=1,zorder=5)
                    if aVmax[i,j]>=113 and aVmax[i,j]<137:
                        ax.plot(alon[i,j], alat[i,j], marker='o', color='firebrick', linewidth=1,zorder=5)
                    if aVmax[i,j]>=137:
                        ax.plot(alon[i,j], alat[i,j], marker='o', color='darkred', linewidth=1,zorder=5)

                    if i==0:
                        ax.text(alon[i,j] + 0.2, alat[i,j] + 0.2, str(atime[j].strftime('%b%d')), fontsize=9, zorder=6)
                else:
                    if aVmax[i,j]<34:
                        ax.plot(alon[i,j], alat[i,j], marker='.', color='navy', linewidth=2,zorder=2)
                    if aVmax[i,j]>=34 and aVmax[i,j]<48:
                        ax.plot(alon[i,j], alat[i,j], marker='.', color='darkgreen', linewidth=2,zorder=2)
                    if aVmax[i,j]>=48 and aVmax[i,j]<64:
                        ax.plot(alon[i,j], alat[i,j], marker='.', color='gold', linewidth=2,zorder=2)
                    if aVmax[i,j]>=64 and aVmax[i,j]<83:
                        ax.plot(alon[i,j], alat[i,j], marker='.', color='orange', linewidth=2,zorder=4)
                    if aVmax[i,j]>=83 and aVmax[i,j]<96:
                        ax.plot(alon[i,j], alat[i,j], marker='.', color='darkorange', linewidth=2,zorder=4)
                    if aVmax[i,j]>=96 and aVmax[i,j]<113:
                        ax.plot(alon[i,j], alat[i,j], marker='.', color='red', linewidth=2,zorder=5)
                    if aVmax[i,j]>=113 and aVmax[i,j]<137:
                        ax.plot(alon[i,j], alat[i,j], marker='.', color='firebrick', linewidth=2,zorder=5)
                    if aVmax[i,j]>=137:
                        ax.plot(alon[i,j], alat[i,j], marker='.', color='darkred', linewidth=2,zorder=5)

            print("  - ctrack : track "+repr(i)+" ok")

    else:

        for i in range(0,alat.shape[0]):
            ax.plot(alon[i,:], alat[i,:], color='firebrick', linewidth=1,zorder=1)
            print("  - ctrack : track "+repr(i)+" ok")

    if ftitle=="no":
        ax.set_title("Cyclone Track")
    else:
        ax.set_title("Cyclone Track "+ftitle)

    plt.tight_layout()
    if flabel=="no":
        figname="CycloneTrack.png"
    else:
        figname="CycloneTrack_"+flabel+".png"

    plt.savefig(figname, dpi=300, facecolor='w', edgecolor='w',
        orientation='portrait', format='png',transparent=False, bbox_inches='tight', pad_inches=0.1)

    plt.close(fig)


def wplot(lat,lon,W,wvar="Hs",flabel="no",ftitle="no"):
    """
    Plot and save figure, single time and fixed domain
    """
    lon[lon<0.]=lon[lon<0.]+360
    levels=np.array(np.arange(0,19,1)).astype('int')
    if str(wvar).upper()=="HS":       
        btitle=str(wvar)+" (m)"
    elif str(wvar).upper()=="TP":
        btitle=str(wvar)+" (s)"

    fig = plt.figure(figsize=(10, 6))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent([lon.min()-0.1,lon.max()+0.1,lat.min()-0.1,lat.max()+0.1], crs=ccrs.PlateCarree())
    ax.coastlines(resolution='10m', linewidth=1)
    ax.add_feature(cfeature.BORDERS, linewidth=0.5)
    ax.add_feature(cfeature.STATES, linewidth=0.5)
    ax.add_feature(cfeature.LAND, facecolor='lightgray')
    ax.add_feature(cfeature.OCEAN, facecolor='white')
    ax.gridlines(draw_labels=True, linestyle='--', color='gray', alpha=0.5)
    cs=ax.contourf(lon,lat,W,levels,cmap=palette,extend="max", zorder=1,transform = ccrs.PlateCarree())
    if ftitle=="no":
        ax.set_title(btitle)
    else:
        ax.set_title(btitle+ftitle)

    plt.tight_layout()
    ax = plt.gca(); pos = ax.get_position()
    l, b, w, h = pos.bounds
    cax = plt.axes([l+0.07, b-0.07, w-0.12, 0.025]) # setup colorbar axes.
    cbar=plt.colorbar(cs,cax=cax, orientation='horizontal'); cbar.ax.tick_params(labelsize=10)
    plt.axes(ax); plt.tight_layout()
    if flabel=="no":
        figname="W_"+str(wvar)+".png"
    else:
        figname="W_"+str(wvar)+"_"+flabel+".png"

    plt.savefig(figname, dpi=200, facecolor='w', edgecolor='w',
        orientation='portrait', format='png',transparent=False, bbox_inches='tight', pad_inches=0.1)

    plt.close(fig)



class TCWaves:

    def __init__(self,Vmax=None,Vfm=None,Rmax=None,R34=None,Lat=None,Lon=None,wdname=None):

        """
        Initialize TCWaves object with required parameters.
        """
        required_inputs = {"Vmax": Vmax,"Vfm": Vfm,"Rmax": Rmax,"R34": R34,"Lat": Lat,"Lon": Lon,"wdname": wdname}
        for name, value in required_inputs.items():
            if value is None:
                raise ValueError(f"Missing required input {name}.")

        self.Vmax = float(Vmax); self.Vfm = float(Vfm)
        self.Rmax = float(Rmax); self.R34 = float(R34)
        self.Lat = float(Lat); self.Lon = float(Lon)
        self.wdname = str(wdname)
        self.wvar = np.array(['Hs','Tp','Uwnd','Vwnd']).astype('str')


    def PWModel(self):
        '''
        Main program, where the wave fields are generated for the vortexes.
        Based on the Matlab code provided by Guisela Grossmann-Matheson
          https://doi.org/10.26188/27237156
        '''
        self.Vmax = min(np.maximum(self.Vmax, 17), 78)
        self.Vfm = min(self.Vfm, 15)
        self.Rmax = min(np.maximum(self.Rmax, 15000), 60000)
        self.R34 = min(np.maximum(self.R34, 200000), 400000)

        wconfig, ds = read(self.wdname)

        # Parameterised Fetch Equation
        F_P32 = (
            wconfig['pmcoeff_a'] * self.Vmax**3 +
            wconfig['pmcoeff_b'] * self.Vmax**2 +
            wconfig['pmcoeff_c'] * self.Vfm**2 +
            wconfig['pmcoeff_d'] * self.Vmax**2 * self.Vfm +
            wconfig['pmcoeff_e'] * self.Vmax * self.Vfm**2 +
            wconfig['pmcoeff_f'] * self.Vmax * self.Vfm +
            wconfig['pmcoeff_g'] * self.Vmax +
            wconfig['pmcoeff_h'] * self.Vfm +
            wconfig['pmcoeff_i'] ) * np.exp(wconfig['pmcoeff_C'] * self.Vfm)

        # Correction factors for Rmax (lambda) and for R34 (gamma)
        clambda = 0.85 * np.log10(self.R34 / 30000.0) + 1.0
        cgamma = 0.65 * np.log10(self.R34 / 300000.0) + 1.0
        # Final Model Fetch 
        F = F_P32*clambda*cgamma/1000.
        # Calculate max value of Hs
        alpha = 0.89 # calibration factor
        Hs_max = alpha * (0.0016 * (9.81 * F * 1000.)**0.5 * self.Vmax) / 9.81

        # Spatial distribution using Wave Diagrams (previously run)
        Z = ds['Z'].values[:] # Nondimensional significant wave height, Hs/Hs(max)
        T = ds['Tp'].values[:] # Peak Wave Period (Tp)
        U = ds['upeak'].values[:] # Zonal component peak direction
        V = ds['vpeak'].values[:] # Meridional component peak direction

        # Interpolation of Wave Diagrams to the exact given conditions
        xi = np.array([[self.Vmax, self.Vfm, self.R34/1000., self.Rmax/1000.]])
        interp = RegularGridInterpolator((ds['Vmax'].values, ds['Vfm'].values, ds['R34'].values, ds['R'].values), Z, method='linear', bounds_error=False, fill_value=None)
        ZZ = interp(xi)[0]; del interp, Z
        interp = RegularGridInterpolator((ds['Vmax'].values, ds['Vfm'].values, ds['R34'].values, ds['R'].values), T, method='linear', bounds_error=False, fill_value=None)
        TP = interp(xi)[0]; del interp, T
        interp = RegularGridInterpolator((ds['Vmax'].values, ds['Vfm'].values, ds['R34'].values, ds['R'].values), U, method='linear', bounds_error=False, fill_value=None)
        UU = interp(xi)[0]; del interp, U
        interp = RegularGridInterpolator((ds['Vmax'].values, ds['Vfm'].values, ds['R34'].values, ds['R'].values), V, method='linear', bounds_error=False, fill_value=None)
        VV = interp(xi)[0]; del interp, V

        # Flip if South Hemisphere
        if self.Lat < 0:
           ZZ = np.fliplr(ZZ.T)
           TP = np.fliplr(TP.T)
           UU = np.fliplr(UU.T)
           VV = np.fliplr(VV.T)

        # Aply Hsmax to nondimensional field
        HS=ZZ*Hs_max 

        result={'xrm':ds['xrm'].values[:],'yrm':ds['yrm'].values[:],
        'Hs':np.array(HS),'Tp':np.array(TP),
        'Uwnd':np.array(UU), 'Vwnd':np.array(VV)}

        return result
        del HS, TP, UU, VV


    def xy_to_lonlat(self,pwm):
        """
        Create the lat lon for the vortexes based on x and y arrays (meters).
        """
        X, Y = np.meshgrid(pwm['xrm'], pwm['yrm'])
        X_flat = X.flatten(); Y_flat = Y.flatten()
        origin = Point(self.Lat, self.Lon)
        lats = []; lons = []
        for dx, dy in zip(X_flat, Y_flat):
            bearing = (np.degrees(np.arctan2(dx, dy)) + 360) % 360
            destination = distance(kilometers=np.hypot(dx, dy)).destination(origin, bearing)
            lats.append(destination.latitude)
            lons.append(destination.longitude)

        pwm['lat'] = np.array(lats).reshape(X.shape)
        pwm['lon'] = np.array(lons).reshape(X.shape)

        return pwm


    def rotate(self,pwm,rangle,wvar="all"):
        """
        Rotate the vortex based on the propagation heading (rangle).
        """
        lon_flat = pwm['lon'].ravel()
        lat_flat = pwm['lat'].ravel()
        lon0 = np.mean(pwm['lon'])
        lat0 = np.mean(pwm['lat'])
        theta = -rangle * np.pi / 180  # degrees to radians, negative for clockwise
        # Approximate lon/lat as Cartesian for small area (equirectangular projection)
        x = (lon_flat - lon0) * np.cos(np.deg2rad(lat0))
        y = lat_flat - lat0
        # Apply 2D rotation
        x_rot = x * np.cos(theta) - y * np.sin(theta)
        y_rot = x * np.sin(theta) + y * np.cos(theta)
        # Convert back to lon/lat
        lon_rot = x_rot / np.cos(np.deg2rad(lat0)) + lon0
        lat_rot = y_rot + lat0
        # Construct rotated points array
        rotated_points = np.column_stack((lon_rot, lat_rot))

        mdist = np.nanmin([np.abs(pwm['yrm'][-1]-np.nanmean(pwm['yrm'])), np.abs(pwm['xrm'][-1]-np.nanmean(pwm['xrm']))])
        X, Y = np.meshgrid(pwm['yrm'], pwm['xrm'])
        gdist = np.sqrt(X**2 + Y**2)
        gdist[gdist>mdist]=np.nan; gdist=gdist/gdist

        # Interpolate onto the original lat/lon grid
        if wvar=="all":
            for i in range(0,len(self.wvar)):
                value = griddata(rotated_points, np.array(pwm[self.wvar[i]]*gdist).ravel(), (pwm['lon'], pwm['lat']), method='linear')
                pwm[self.wvar[i]]=value
                del value

        else:
            if str(wvar).upper()=="HS" or str(wvar).upper()=="SWH":
                value = griddata(rotated_points, np.array(pwm['Hs']*gdist).ravel(), (pwm['lon'], pwm['lat']), method='linear')
                pwm['Hs']=value
                del value

            elif str(wvar).upper()=="TP":
                value = griddata(rotated_points, np.array(pwm['Tp']*gdist).ravel(), (pwm['lon'], pwm['lat']), method='linear')
                pwm['Tp']=value
                del value

        return pwm


    def blend(self,pwm,W,wvar="all"):
        """
        Interpolate and blend the vortexes onto a fixed given grid.
        """

        auxlon = np.array(pwm['lon']); auxlon[auxlon<0]=auxlon[auxlon<0]+360.
        # Interpolate and collocate into the final grid and domain
        tcgrid = np.column_stack((pwm['lat'].ravel(), auxlon.ravel()))
        flat_grid, flon_grid = np.meshgrid(W['lat'], W['lon'], indexing='ij')
        target_grid = np.column_stack((flat_grid.ravel(), flon_grid.ravel()))
        if wvar=="all":
            WFHs = griddata(tcgrid, pwm['Hs'].ravel(), target_grid, method='linear', fill_value=np.nan).reshape(W['Hs'].shape)
            WFTp = griddata(tcgrid, pwm['Tp'].ravel(), target_grid, method='linear', fill_value=np.nan).reshape(W['Tp'].shape)
            WFUwnd = griddata(tcgrid, pwm['Uwnd'].ravel(), target_grid, method='linear', fill_value=np.nan).reshape(W['Uwnd'].shape)
            WFVwnd = griddata(tcgrid, pwm['Vwnd'].ravel(), target_grid, method='linear', fill_value=np.nan).reshape(W['Vwnd'].shape)
            WW={'lat':W['lat'],'lon':W['lon'],'Hs':WFHs, 'Tp':WFTp, 'Uwnd':WFUwnd, 'Vwnd':WFVwnd}
        else:
            if str(wvar).upper()=="HS" or str(wvar).upper()=="SWH":
                WFHs = griddata(tcgrid, pwm['Hs'].ravel(), target_grid, method='linear', fill_value=np.nan).reshape(W['Hs'].shape)
                WW={'lat':W['lat'],'lon':W['lon'],'Hs':WFHs}
            elif str(wvar).upper()=="TP":
                WFTp = griddata(tcgrid, pwm['Tp'].ravel(), target_grid, method='linear', fill_value=np.nan).reshape(W['Tp'].shape)
                WW={'lat':W['lat'],'lon':W['lon'],'Tp':WFTp}

        return pwm,WW

    def hwplot(self,pwm,WW,wvar="Hs",flabel="no",ftitle="no"):
        """
        Plot and save figures (Hs and Tp), for both the cyclones (following) and the large domain (fixed).
        """
        # Plot following the Cyclone ---
        # Hs
        lmax = np.maximum(int(np.ceil(np.nanmax(pwm['Hs']))),8)
        levels=np.array(np.arange(0,lmax+1,1)).astype('int')
        plt.figure(figsize=(6,6))
        # ax = plt.axes(projection=ccrs.PlateCarree())
        ax = plt.axes(projection=ccrs.Mercator())
        ax.set_extent([pwm['lon'].min()-0.1,pwm['lon'].max()+0.1,pwm['lat'].min()-0.1,pwm['lat'].max()+0.1], crs=ccrs.PlateCarree())
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='grey', alpha=0.5, linestyle='--', zorder=3)
        gl.xlabel_style = {'size': 9, 'color': 'k','rotation':0}; gl.ylabel_style = {'size': 9, 'color': 'k','rotation':0}
        cs=ax.contourf(pwm['lon'],pwm['lat'],pwm['Hs'],levels,cmap=palette,extend="max", zorder=1,transform = ccrs.PlateCarree())
        ax.add_feature(cartopy.feature.OCEAN,facecolor=("white"))
        ax.add_feature(cartopy.feature.LAND,facecolor=("lightgrey"), edgecolor='grey',linewidth=0.5, zorder=2)
        ax.add_feature(cartopy.feature.BORDERS, edgecolor='grey', linestyle='-',linewidth=0.5, alpha=1, zorder=3)
        ax.coastlines(resolution='50m', color='grey',linewidth=0.5, linestyle='-', alpha=1, zorder=4)
        if ftitle=="no":
            ax.set_title("Hs (m)")
        else:
            ax.set_title("Hs (m) "+ftitle)

        plt.tight_layout()
        ax = plt.gca(); pos = ax.get_position()
        l, b, w, h = pos.bounds
        cax = plt.axes([l+0.07, b-0.07, w-0.12, 0.025]) # setup colorbar axes.
        cbar=plt.colorbar(cs,cax=cax, orientation='horizontal'); cbar.ax.tick_params(labelsize=10)
        # tick_locator = ticker.MaxNLocator(nbins=7); cbar.locator = tick_locator; cbar.update_ticks()
        plt.axes(ax); plt.tight_layout()
        if flabel=="no":
            figname="TCw_Hs.png"
        else:
            figname="TCw_Hs_"+flabel+".png"

        plt.savefig(figname, dpi=200, facecolor='w', edgecolor='w',
            orientation='portrait', format='png',transparent=False, bbox_inches='tight', pad_inches=0.1)

        if str(wvar).upper()!="HS":
            # Tp
            lmax = np.maximum(int(np.ceil(np.nanmax(pwm['Tp']))),8)
            levels=np.array(np.arange(0,lmax+1,1)).astype('int')
            cs=ax.contourf(pwm['lon'],pwm['lat'],pwm['Tp'],levels,cmap=palette,extend="max", zorder=1,transform = ccrs.PlateCarree())
            if ftitle=="no":
                ax.set_title("Tp (s)")
            else:
                ax.set_title("Tp (s) "+ftitle)

            plt.tight_layout()
            cbar=plt.colorbar(cs,cax=cax, orientation='horizontal'); cbar.ax.tick_params(labelsize=10)
            plt.axes(ax); plt.tight_layout()
            if flabel=="no":
                figname="TCw_Tp.png"
            else:
                figname="TCw_Tp_"+flabel+".png"

            plt.savefig(figname, dpi=200, facecolor='w', edgecolor='w',
                orientation='portrait', format='png',transparent=False, bbox_inches='tight', pad_inches=0.1)

            plt.close('all'); del ax

        # Plot Fixed Domain ---
        # Hs
        lmax = np.maximum(int(np.ceil(np.nanmax(WW['Hs']))),8)
        levels=np.array(np.arange(0,lmax+1,1)).astype('int')
        lonp=WW['lon']; lonp[lonp<0]=lonp[lonp<0]+360
        plt.figure(figsize=(9,5.5))
        ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=lonp.mean()))
        ax.set_extent([lonp.min(),lonp.max()+0.1,WW['lat'].min()-0.1,WW['lat'].max()+0.1], crs=ccrs.PlateCarree())
        gl = ax.gridlines(crs=ccrs.PlateCarree(),  xlocs=range(-180,180, 20), draw_labels=True, linewidth=0.5, color='grey', alpha=0.5, linestyle='--')
        gl.xlabel_style = {'size': 9, 'color': 'k','rotation':0}; gl.ylabel_style = {'size': 9, 'color': 'k','rotation':0}
        cs=ax.contourf(lonp,WW['lat'],WW['Hs'],levels=levels,cmap='jet',zorder=1,extend="max",transform = ccrs.PlateCarree())
        ax.add_feature(cartopy.feature.OCEAN,facecolor=("white"))
        ax.add_feature(cartopy.feature.LAND,facecolor=("lightgrey"), edgecolor='grey',linewidth=0.5, zorder=2)
        ax.add_feature(cartopy.feature.BORDERS, edgecolor='grey', linestyle='-',linewidth=0.5, alpha=1, zorder=3)
        ax.coastlines(resolution='110m', color='dimgrey',linewidth=0.5, linestyle='-', alpha=1, zorder=4)
        if ftitle=="no":
            ax.set_title("Hs (m)")
        else:
            ax.set_title("Hs (m) "+ftitle)

        plt.tight_layout()
        ax = plt.gca(); pos = ax.get_position(); l, b, w, h = pos.bounds; cax = plt.axes([l+0.06, b-0.07, w-0.12, 0.03]) # setup colorbar axes.
        cbar = plt.colorbar(cs,cax=cax, orientation='horizontal', format='%g')
        plt.axes(ax); plt.tight_layout()
        if flabel=="no":
            figname="Wfield_Hs.png"
        else:
            figname="Wfield_Hs_"+flabel+".png"

        plt.savefig(figname, dpi=200, facecolor='w', edgecolor='w',
            orientation='portrait', format='png',transparent=False, bbox_inches='tight', pad_inches=0.1)

        if str(wvar).upper()!="HS":
            # Tp
            lmax = np.maximum(int(np.ceil(np.nanmax(WW['Tp']))),8)
            levels=np.array(np.arange(0,lmax+1,1)).astype('int')
            cs=ax.contourf(lonp,WW['lat'],WW['Tp'],levels=levels,cmap='jet',zorder=1,extend="max",transform = ccrs.PlateCarree())
            plt.tight_layout()
            cbar = plt.colorbar(cs,cax=cax, orientation='horizontal', format='%g')
            if ftitle=="no":
                ax.set_title("Tp (s)")
            else:
                ax.set_title("Tp (s) "+ftitle)

            plt.axes(ax); plt.tight_layout()
            if flabel=="no":
                figname="Wfield_Tp.png"
            else:
                figname="Wfield_Tp_"+flabel+".png"

            plt.savefig(figname, dpi=200, facecolor='w', edgecolor='w',
                orientation='portrait', format='png',transparent=False, bbox_inches='tight', pad_inches=0.1)

            plt.close('all'); del ax


