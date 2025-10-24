#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from matplotlib.mlab import *
import pandas as pd
from pylab import *
import os
import netCDF4 as nc
import pyresample
import time
import timeit
from calendar import timegm
import sys
import warnings; warnings.filterwarnings("ignore")
# netcdf format
fnetcdf="NETCDF4"


# weight function for pyresample
def wf(pdist):
    a=(1 - pdist / (dlim+1))
    return (abs(a)+a)/2


if __name__ == "__main__":

    # INPUTS 
    # power of initial array 10**pia (size) that will be used to allocate satellite data (faster than append)
    pia=10
    # Maximum distance (m) for pyresample weighted average, See: https://doi.org/10.3390/rs15082203
    dlim=35000.
    # Maximum temporal distance (s) for pyresample weighted average
    maxti=1800.
    # Directory where AODN altimeter data is saved, downloaded using wfetchsatellite_AODN_Altimeter.sh
    dirs='/work/noaa/marine/ricardo.campos/data/AODN/altimeter'

    # Satellite missions available at AODN dataset, pick one as this code runs one satellite at a time!
    if len(sys.argv) <= 2 :
        s=int(sys.argv[1]) # argument satellite ID for satellite mission selection. 
        # s=0 is JASON3, s=1 is JASON2 etc. See list below in sdname
    else:
        s=int(sys.argv[1])
        nseg=int(sys.argv[2])
        seg=int(sys.argv[3])

    sdname=np.array(['JASON3','JASON2','CRYOSAT2','JASON1','HY2','HY2B','SARAL','SENTINEL3A','ENVISAT','ERS1','ERS2','GEOSAT','GFO','TOPEX','SENTINEL3B','CFOSAT','SENTINEL6A'])
    sname=np.array(['JASON-3','JASON-2','CRYOSAT-2','JASON-1','HY-2','HY-2B','SARAL','SENTINEL-3A','ENVISAT','ERS-1','ERS-2','GEOSAT','GFO','TOPEX','SENTINEL-3B','CFOSAT','SENTINEL-6A'])
    # Ongoing sat missions:
    # CFOSAT, SARAL, CRYOSAT2, HY2B, JASON3, SENTINEL3A, SENTINEL3B, SENTINEL6A
    # 0, 2, 5, 6, 7, 14, 15, 16

    # Quality Control parameters
    max_swh_rms = 1.5  # Max RMS of the band significant wave height
    max_sig0_rms = 0.8 # Max RMS of the backscatter coefficient
    max_swh_qc = 2.0 # Max SWH Ku band quality control
    hsmax=20.; wspmax=90.
    min_swh_numval = np.array([17,17,17,17,17,17,17,17,17,17,17,-9999,3,7,17,-9999,17])
    # ---------

    start = timeit.default_timer()

    print(" Processing segment "+str(seg)+" of total "+str(nseg)+". SAT "+sdname[s])

    # Read buoys Reference
    df = pd.read_csv('Data_REF.txt', sep='\t')
    lines_to_read = len(df) // nseg
    if seg<nseg:
        df = df.iloc[(seg-1)*lines_to_read:(seg)*lines_to_read, :]
    else:
        df = df.iloc[(seg-1)*lines_to_read::, :]

    blat = np.array(df['lat'][:]); blon = np.array(df['lon'][:])
    ftime = np.array( (pd.to_datetime(df['time'][:], format='%Y%m%d%H%M') - pd.Timestamp('1970-01-01')) // pd.Timedelta('1s') ).astype('double')
    btime = np.array( (pd.to_datetime(df['buoy_time'][:], format='%Y%m%d%H%M') - pd.Timestamp('1970-01-01')) // pd.Timedelta('1s') ).astype('double')
    # ---------------

    # Date interval for satellite data allocation
    adatemin=btime.min()-3600.; adatemax=btime.max()+3600.
    adatemin=(adatemin-float(timegm( time.strptime('1985010100', '%Y%m%d%H') )))/24./3600.
    adatemax=(adatemax-float(timegm( time.strptime('1985010100', '%Y%m%d%H') )))/24./3600.


    # Read Sat Data -----------------
    #  domain of interest
    auxlat=np.array(np.arange(np.floor(blat.min())-1.,np.ceil(blat.max())+1.,1)).astype('int')
    ablon=np.copy(blon); ablon[ablon<0.]=ablon[ablon<0.]+360.
    auxlon=np.array(np.arange(np.floor(ablon.min())-1.,np.ceil(ablon.max())+1.,1)).astype('int')

    # Read and allocate satellite data into arrays
    ast=np.double(np.zeros((10**pia),'d')); aslat=np.zeros((10**pia),'f'); aslon=np.zeros((10**pia),'f');
    ahskcal=np.zeros((10**pia),'f')
    awndcal=np.zeros((10**pia),'f'); asig0knstd=np.zeros((10**pia),'f');
    aswhknobs=np.zeros((10**pia),'f'); aswhknstd=np.zeros((10**pia),'f'); aswhkqc=np.zeros((10**pia),'f')
    ii=0
    for j in auxlat:
        for k in auxlon:

            if j>=0:
                hem='N'
            else:
                hem='S'

            try: 
                fu=nc.Dataset(dirs+'/'+sdname[s]+'/IMOS_SRS-Surface-Waves_MW_'+sname[s]+'_FV02_'+str(np.abs(j)).zfill(3)+hem+'-'+str(k).zfill(3)+'E-DM00.nc')
            except:
                print(dirs+'/'+sdname[s]+'/IMOS_SRS-Surface-Waves_MW_'+sname[s]+'_FV02_'+str(np.abs(j)).zfill(3)+hem+'-'+str(k).zfill(3)+'E-DM00.nc does not exist'); vai=0
            else:
                st=np.double(fu.variables['TIME'][:])
                if np.size(st)>10:
                    slat=fu.variables['LATITUDE'][:]
                    slon=fu.variables['LONGITUDE'][:]
                    wndcal=fu.variables['WSPD_CAL'][:]
                    try: 
                        hskcal=fu.variables['SWH_KU_CAL'][:]
                        sig0knstd=fu.variables['SIG0_KU_std_dev'][:]
                        swhknobs=fu.variables['SWH_KU_num_obs'][:]
                        swhknstd=fu.variables['SWH_KU_std_dev'][:]
                        swhkqc=fu.variables['SWH_KU_quality_control'][:]
                    except:
                        print(' error reading KU, pick KA')
                        hskcal=fu.variables['SWH_KA_CAL'][:]
                        sig0knstd=fu.variables['SIG0_KA_std_dev'][:]
                        swhknobs=fu.variables['SWH_KA_num_obs'][:]
                        swhknstd=fu.variables['SWH_KA_std_dev'][:]
                        swhkqc=fu.variables['SWH_KA_quality_control'][:]

                    if ii+np.size(st) <= ast.shape[0] :
                        if (st.shape[0]==wndcal.shape[0]) & (slat.shape[0]==slon.shape[0]) & (wndcal.shape[0]==hskcal.shape[0]) :    
                            ast[ii:ii+st.shape[0]]=np.array(st).astype('double')
                            aslat[ii:ii+st.shape[0]]=np.array(slat).astype('float')
                            aslon[ii:ii+st.shape[0]]=np.array(slon).astype('float')
                            ahskcal[ii:ii+st.shape[0]]=np.array(hskcal).astype('float')
                            awndcal[ii:ii+st.shape[0]]=np.array(wndcal).astype('float')
                            asig0knstd[ii:ii+st.shape[0]]=np.array(sig0knstd).astype('float')
                            aswhknobs[ii:ii+st.shape[0]]=np.array(swhknobs).astype('float')
                            aswhknstd[ii:ii+st.shape[0]]=np.array(swhknstd).astype('float')
                            aswhkqc[ii:ii+st.shape[0]]=np.array(swhkqc).astype('float')
                            ii=ii+st.shape[0]

                    else:
                        sys.exit('Small array to allocate the satellite data! Increase the power of initial array (pia)')

                    del st,slat,slon,hskcal,wndcal,sig0knstd,swhknobs,swhknstd,swhkqc
                    fu.close(); del fu

    print(' Done reading and allocating satellite data '+sdname[s])
    del ii

    # Quality Control Check ----
    indq = np.where( (aswhknstd<=max_swh_rms) & (asig0knstd<=max_sig0_rms) & (aswhknobs>=min_swh_numval[s]) & (aswhkqc<=max_swh_qc) & (ahskcal>0.3) & (ahskcal<hsmax) & (awndcal>0.3) & (awndcal<wspmax) & (ast>=adatemin) & (ast<=adatemax) )
    del asig0knstd,aswhknobs,aswhknstd,aswhkqc,adatemin,adatemax
    print(' '); print(" Collocation ")
    if np.size(indq)>10:
        ast=np.double(np.copy(ast[indq[0]]))
        ast=np.double(np.copy(ast)*24.*3600.+float(timegm( time.strptime('1985010100', '%Y%m%d%H') )))
        aslat=np.copy(aslat[indq[0]]); aslon=np.copy(aslon[indq[0]])
        ahskcal=np.copy(ahskcal[indq[0]])
        awndcal=np.copy(awndcal[indq[0]])
        del indq

        # Collocation -------------------
        fhskcal=np.zeros((btime.shape[0]),'f')*np.nan; fwndcal=np.zeros((btime.shape[0]),'f')*np.nan
        fhskcaln=np.zeros((btime.shape[0]),'f')*np.nan; fwndcaln=np.zeros((btime.shape[0]),'f')*np.nan

        prlon=np.copy(aslon); prlon[prlon>180.]=prlon[prlon>180.]-360.

        ubtime=np.unique(btime)
        for t in range(0,ubtime.shape[0]):

            indt1 = np.where( abs(ast[:]-ubtime[t]) < maxti )
            indt2 = np.where( btime[:]==ubtime[t] )[0]

            if np.size(indt1)>0:
                indt1=indt1[0]

                # Orig space
                orig_def = pyresample.geometry.SwathDefinition(lons=prlon[indt1], lats=aslat[indt1])
                # Target
                targ_def = pyresample.geometry.SwathDefinition(lons=blon[indt2],lats=blat[indt2])

                # By distance function wf
                fhskcal[indt2] = pyresample.kd_tree.resample_custom(orig_def,ahskcal[indt1],targ_def,radius_of_influence=dlim,weight_funcs=wf,fill_value=0,nprocs=0)
                fwndcal[indt2] = pyresample.kd_tree.resample_custom(orig_def,awndcal[indt1],targ_def,radius_of_influence=dlim,weight_funcs=wf,fill_value=0,nprocs=0)
                # nearest
                fhskcaln[indt2] = pyresample.kd_tree.resample_nearest(orig_def,ahskcal[indt1],targ_def,radius_of_influence=dlim,fill_value=0,nprocs=0)
                fwndcaln[indt2] = pyresample.kd_tree.resample_nearest(orig_def,awndcal[indt1],targ_def,radius_of_influence=dlim,fill_value=0,nprocs=0)

                del indt1, indt2, orig_def, targ_def

        print(sdname[s]+' Collocation Done')

        fhskcal[(fhskcal < 0.3) | (fhskcal > hsmax)] = -999.999
        fhskcal[np.isnan(fhskcal)==True]=-999.999
        fwndcal[(fwndcal < 1.) | (fwndcal > wspmax)] = -999.999
        fwndcal[np.isnan(fwndcal)==True]=-999.999
        fhskcaln[(fhskcaln < 0.3) | (fhskcaln > hsmax)] = -999.999
        fhskcaln[np.isnan(fhskcaln)==True]=-999.999
        fwndcaln[(fwndcaln < 1.) | (fwndcaln > wspmax)] = -999.999
        fwndcaln[np.isnan(fwndcaln)==True]=-999.999

        # Save final data
        df = pd.DataFrame({
            'time': pd.to_datetime(ftime, unit='s').strftime('%Y%m%d%H%M'),
            'buoy_time': pd.to_datetime(btime, unit='s').strftime('%Y%m%d%H%M'),
            'lat': np.round(blat,5),
            'lon': np.round(blon,5),
            'hs_avr': np.round(fhskcal,4),
            'wnd_avr': np.round(fwndcal,4),
            'hs_nrst': np.round(fhskcaln,4),
            'wnd_nrst': np.round(fwndcaln,4)
        })

        fname="Data_REF_"+sdname[s]+"_"+str(seg).zfill(2)+"_"+str(nseg).zfill(2)+".txt"
        if seg==1:
            df.to_csv(fname, sep='\t', index=False, header=True)
        else:
            df.to_csv(fname, sep='\t', index=False, header=False)

    stop = timeit.default_timer()
    print('Concluded in '+repr(int(round(stop - start,0)))+' seconds')

