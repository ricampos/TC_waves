#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
 Support conversion of matlab to python scripts.
 Read .mat files (wave diagrams), organize arrays, and save in netCDF4 format
"""

import numpy as np
from matplotlib.mlab import *
import scipy.io
import netCDF4 as nc
import yaml
fnetcdf="NETCDF4"

if __name__ == "__main__":

    pmatpath="/home/ricardo/work/noaa/analysis/TC_Waves/modeling/PModel/PWModel_Matlab_package/Matlab_example_for_Paper/wave_diagrams"

    print(" ")
    with open('/home/ricardo/work/noaa/github/TC_waves/pmodel/pmodel_config.yaml', 'r') as file:
        wconfig = yaml.safe_load(file)

    print(" Reading yaml configuration, ok")

    # 378 model runs with the spatial distribution of wave parameters
    # delta-pressure, Vfm, R34, Rmax
    Hs=np.zeros( (len(wconfig['Dp_runs']), len(wconfig['Vfm_runs']), len(wconfig['R34_runs']), len(wconfig['R_runs']),281,281), 'f')*np.nan
    T=np.zeros( (len(wconfig['Dp_runs']), len(wconfig['Vfm_runs']), len(wconfig['R34_runs']), len(wconfig['R_runs']),281,281), 'f')*np.nan
    Tp=np.zeros( (len(wconfig['Dp_runs']), len(wconfig['Vfm_runs']), len(wconfig['R34_runs']), len(wconfig['R_runs']),281,281), 'f')*np.nan
    Z=np.zeros( (len(wconfig['Dp_runs']), len(wconfig['Vfm_runs']), len(wconfig['R34_runs']), len(wconfig['R_runs']),281,281), 'f')*np.nan
    upeak=np.zeros( (len(wconfig['Dp_runs']), len(wconfig['Vfm_runs']), len(wconfig['R34_runs']), len(wconfig['R_runs']),281,281), 'f')*np.nan
    vpeak=np.zeros( (len(wconfig['Dp_runs']), len(wconfig['Vfm_runs']), len(wconfig['R34_runs']), len(wconfig['R_runs']),281,281), 'f')*np.nan

    c=0
    for i in range(0,len(wconfig['Dp_runs'])):
        for j in range(0,len(wconfig['Vfm_runs'])):
            for k in range(0,len(wconfig['R34_runs'])):
                for l in range(0,len(wconfig['R_runs'])):

                    # Load .mat file
                    data = scipy.io.loadmat(pmatpath+"/wave_"+str(int(wconfig['Dp_runs'][i]))+"_"+str(int(wconfig['Vfm_runs'][j]*10))+"_"+str(wconfig['R34_runs'][k])+"_"+str(wconfig['R_runs'][l])+".mat")

                    Hs[i,j,k,l,:,:]=data['Hs']
                    T[i,j,k,l,:,:]=data['T']
                    Tp[i,j,k,l,:,:]=data['Tp']
                    Z[i,j,k,l,:,:]=data['Z'] # nondimensional significant wave height
                    upeak[i,j,k,l,:,:]=data['upeak'] # zonal component peak direction
                    vpeak[i,j,k,l,:,:]=data['vpeak'] # meridional component peak direction
                    if c==0:
                        xrm = np.round(data['xrm'][:,0]*30.) # distance in x direction,
                        yrm = np.round(data['yrm'][:,0]*30.) # distance in y direction

                    del data
                    c=c+1

    print(" Wave Diagrams, reading and organizing, ok. Total of "+repr(c)+" diagrams.")  

    # Save netCDF4 output file
    ncfile = nc.Dataset('PModel_WaveDiagrams.nc', "w", format=fnetcdf)
    ncfile.history="PModel (Grossmann-Matheson et al., 2025) spatial distributions of wave parameters in Tropical Cyclones."
    # create  dimensions
    ncfile.createDimension('Delta_pressure', len(wconfig['Dp_runs']) )
    ncfile.createDimension('Vfm', len(wconfig['Vfm_runs']) )
    ncfile.createDimension('R34', len(wconfig['R34_runs']) )
    ncfile.createDimension('R', len(wconfig['R_runs']) )
    ncfile.createDimension('grid', len(xrm) )
    # create variables
    vDelta_pressure = ncfile.createVariable('Delta_pressure',np.float32,('Delta_pressure'))
    vVmax = ncfile.createVariable('Vmax',np.float32,('Delta_pressure'))
    vVfm = ncfile.createVariable('Vfm',np.float32,('Vfm'))
    vR34 = ncfile.createVariable('R34',np.float32,('R34'))
    vR = ncfile.createVariable('R',np.float32,('R'))
    # results
    vHs = ncfile.createVariable('Hs',np.float32,('Delta_pressure','Vfm','R34','R','grid','grid'))
    vT = ncfile.createVariable('T',np.float32,('Delta_pressure','Vfm','R34','R','grid','grid'))
    vTp = ncfile.createVariable('Tp',np.float32,('Delta_pressure','Vfm','R34','R','grid','grid'))
    vZ = ncfile.createVariable('Z',np.float32,('Delta_pressure','Vfm','R34','R','grid','grid'))
    vupeak = ncfile.createVariable('upeak',np.float32,('Delta_pressure','Vfm','R34','R','grid','grid'))
    vvpeak = ncfile.createVariable('vpeak',np.float32,('Delta_pressure','Vfm','R34','R','grid','grid'))
    vxrm = ncfile.createVariable('xrm',np.float32,('grid'))
    vyrm = ncfile.createVariable('yrm',np.float32,('grid'))
    # Adding long names to the variables
    vDelta_pressure.long_name = "Central pressure drop (pambiant-p0)"
    vVmax.long_name = "Maximum sustained wind speed"
    vVfm.long_name = "Velocity of forward movement"
    vR34.long_name = "Radius to gales (34 kt)"
    vR.long_name = "Radius to maximum winds"
    vHs.long_name = "Significant Wave Height (Hs)"
    vT.long_name = "Nondimensional Wave Period (T)"
    vTp.long_name = "Peak Wave Period (Tp)"
    vZ.long_name = "Nondimensional significant wave height, Hs/Hs(max)"
    vupeak.long_name = "Zonal component peak direction"
    vvpeak.long_name = "Meridional component peak direction"
    vxrm.long_name = "distance in x direction (lon)"
    vyrm.long_name = "distance in y direction (lat)"
    # Assign units
    vDelta_pressure.units = 'hPa'
    vVmax.units = 'm/s'
    vVfm.units = 'm/s'
    vR34.units='km'; vR.units='km'
    vHs.units='m'; 	vZ.units='nondimensional, 0 to 1'
    vT.units='nondimensional, 0 to 1'; vTp.units='s'
    vxrm.units='km'
    vyrm.units='km'

    # Allocate Data
    vDelta_pressure[:] = np.array(wconfig['Dp_runs']).astype('float')[:]
    vVmax[:] = np.array(wconfig['Vmax_runs']).astype('float')[:]
    vVfm[:] = np.array(wconfig['Vfm_runs']).astype('float')[:]
    vR34[:] = np.array(wconfig['R34_runs']).astype('float')[:]
    vR[:] = np.array(wconfig['R_runs']).astype('float')[:]
    vHs[:,:,:,:,:,:] = Hs[:,:,:,:,:,:]
    vT[:,:,:,:,:,:] = T[:,:,:,:,:,:]
    vTp[:,:,:,:,:,:] = Tp[:,:,:,:,:,:]
    vZ[:,:,:,:,:,:] = Z[:,:,:,:,:,:]
    vupeak[:,:,:,:,:,:] = upeak[:,:,:,:,:,:]
    vvpeak[:,:,:,:,:,:] = vpeak[:,:,:,:,:,:]
    vxrm[:] = xrm[:]
    vyrm[:] = yrm[:]
    #
    ncfile.close()
    print(" Output netcdf file saved, ok")


