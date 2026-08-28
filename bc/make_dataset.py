#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Build dataset for ML modeling
Total of 129 variables with ~32k matchups for training, distributed over 4 years and 84 TCs.
"""

import numpy as np
import pandas as pd
import netCDF4 as nc
import warnings; warnings.filterwarnings("ignore")
fnetcdf="NETCDF4"

def read_ibtracs(fname):
    '''
    Read IBTracks text data
    Only named storms (TS) are selected
    IBTRACKS V4 data. It can be downloaded at
    https://www.ncei.noaa.gov/data/international-best-track-archive-for-climate-stewardship-ibtracs/v04r00/access/csv/
    https://www.ncdc.noaa.gov/ibtracs/index.php?name=bib
    '''
    dfibtr = pd.read_csv(fname, header=[0,1])
    inat = np.array(dfibtr.values[:,7]); iname=np.array(dfibtr['NAME'].values[:][:,0]).astype('str')
    icyid = np.array(dfibtr['NUMBER'].values[:][:,0]).astype('int')
    ity = np.array(pd.Series(dfibtr['ISO_TIME'].values.flatten()).str.slice(0,4)).astype('int')
    itime = pd.to_datetime(dfibtr['ISO_TIME'].values[:,0])
    ilat = np.array(dfibtr['LAT'].values[:,0]).astype('float')
    ilon = np.array(dfibtr['LON'].values[:,0]).astype('float')
    iVmax = np.array(dfibtr['USA_WIND'].replace(' ', np.nan).values[:,0]).astype('float')/1.94384 # max wind speed in m/s
    iVfm = np.array(dfibtr['STORM_SPEED'].replace(' ', np.nan).values[:,0]).astype('float')/1.94384 # forward speed in m/s
    iDir = np.array(dfibtr['STORM_DIR'].replace(' ', np.nan).values[:,0]).astype('float')

    iRmax = np.array(dfibtr['USA_RMW'].replace(' ', np.nan).values[:,0]).astype('float')*1852. # in meters

    iR34_NE = np.array(dfibtr['USA_R34_NE'].replace(' ', np.nan).values[:,0]).astype('float')*1852.; iR34_NE[np.isnan(iR34_NE)==True]=0.
    iR34_NW = np.array(dfibtr['USA_R34_NW'].replace(' ', np.nan).values[:,0]).astype('float')*1852.; iR34_NW[np.isnan(iR34_NW)==True]=0.
    iR34_SW = np.array(dfibtr['USA_R34_SW'].replace(' ', np.nan).values[:,0]).astype('float')*1852.; iR34_SW[np.isnan(iR34_SW)==True]=0.
    iR34_SE = np.array(dfibtr['USA_R34_SE'].replace(' ', np.nan).values[:,0]).astype('float')*1852.; iR34_SE[np.isnan(iR34_SE)==True]=0.
    iR34 = np.mean(np.array([iR34_NE,iR34_NW,iR34_SW,iR34_SE]),axis=0)

    iR50_NE = np.array(dfibtr['USA_R50_NE'].replace(' ', np.nan).values[:,0]).astype('float')*1852.; iR50_NE[np.isnan(iR50_NE)==True]=0.
    iR50_NW = np.array(dfibtr['USA_R50_NW'].replace(' ', np.nan).values[:,0]).astype('float')*1852.; iR50_NW[np.isnan(iR50_NW)==True]=0.
    iR50_SW = np.array(dfibtr['USA_R50_SW'].replace(' ', np.nan).values[:,0]).astype('float')*1852.; iR50_SW[np.isnan(iR50_SW)==True]=0.
    iR50_SE = np.array(dfibtr['USA_R50_SE'].replace(' ', np.nan).values[:,0]).astype('float')*1852.; iR50_SE[np.isnan(iR50_SE)==True]=0.
    iR50 = np.mean(np.array([iR50_NE,iR50_NW,iR50_SW,iR50_SE]),axis=0)

    iR64_NE = np.array(dfibtr['USA_R64_NE'].replace(' ', np.nan).values[:,0]).astype('float')*1852.; iR64_NE[np.isnan(iR64_NE)==True]=0.
    iR64_NW = np.array(dfibtr['USA_R64_NW'].replace(' ', np.nan).values[:,0]).astype('float')*1852.; iR64_NW[np.isnan(iR64_NW)==True]=0.
    iR64_SW = np.array(dfibtr['USA_R64_SW'].replace(' ', np.nan).values[:,0]).astype('float')*1852.; iR64_SW[np.isnan(iR64_SW)==True]=0.
    iR64_SE = np.array(dfibtr['USA_R64_SE'].replace(' ', np.nan).values[:,0]).astype('float')*1852.; iR64_SE[np.isnan(iR64_SE)==True]=0.
    iR64 = np.mean(np.array([iR64_NE,iR64_NW,iR64_SW,iR64_SE]),axis=0)

    print('read_ibtracs: Ibtracks ok'); del dfibtr
    ind=np.where((iname!="UNNAMED") & (inat=='TS'))
    ity=np.copy(ity[ind[0]]); itime=np.copy(itime[ind[0]])
    icyid=np.copy(icyid[ind[0]]); iname=np.copy(iname[ind[0]]); ilat=np.copy(ilat[ind[0]]);ilon=np.copy(ilon[ind[0]]);inat=np.copy(inat[ind[0]])
    iVmax=np.copy(iVmax[ind[0]]); iVfm=np.copy(iVfm[ind[0]]); iDir=np.copy(iDir[ind[0]])
    iRmax=np.copy(iRmax[ind[0]])
    iR34=np.copy(iR34[ind[0]]); iR34_NE=np.copy(iR34_NE[ind[0]]); iR34_NW=np.copy(iR34_NW[ind[0]]); iR34_SW=np.copy(iR34_SW[ind[0]]); iR34_SE=np.copy(iR34_SE[ind[0]])
    iR50=np.copy(iR50[ind[0]]); iR50_NE=np.copy(iR50_NE[ind[0]]); iR50_NW=np.copy(iR50_NW[ind[0]]); iR50_SW=np.copy(iR50_SW[ind[0]]); iR50_SE=np.copy(iR50_SE[ind[0]])
    iR64=np.copy(iR64[ind[0]]); iR64_NE=np.copy(iR64_NE[ind[0]]); iR64_NW=np.copy(iR64_NW[ind[0]]); iR64_SW=np.copy(iR64_SW[ind[0]]); iR64_SE=np.copy(iR64_SE[ind[0]])

    del ind
    # ilon[ilon<0]=ilon[ilon<0]+360.

    icyid = np.char.add(ity.astype(int).astype(str), np.char.zfill(icyid.astype(int).astype(str), 3))

    wibtr={'itime':pd.to_datetime(itime),'icyid':np.array(icyid).astype('int'),'iname':iname,'ilat':ilat,'ilon':ilon,'inat':inat,
        'iVmax':iVmax,'iVfm':iVfm,'iDir':iDir,'iRmax':iRmax,
        'iR34':iR34, 'iR34_NE':iR34_NE, 'iR34_NW':iR34_NW, 'iR34_SE':iR34_SE, 'iR34_SW':iR34_SW, 
        'iR50':iR50, 'iR50_NE':iR50_NE, 'iR50_NW':iR50_NW, 'iR50_SE':iR50_SE, 'iR50_SW':iR50_SW,
        'iR64':iR64, 'iR64_NE':iR64_NE, 'iR64_NW':iR64_NW, 'iR64_SE':iR64_SE, 'iR64_SW':iR64_SW}

    return wibtr


def ot_to_epoch_seconds(ot):
    """Convert integer times like 202408041500 (YYYYMMDDHHMM) to seconds since 1970-01-01."""
    dt = pd.to_datetime(np.asarray(ot).astype(str), format="%Y%m%d%H%M")
    epoch = (dt - pd.Timestamp("1970-01-01")) // pd.Timedelta("1s")
    return np.asarray(epoch, dtype="int64")
 

if __name__ == "__main__":

    # Read Obs and Model
    dirdata="/home/ricardo/work/noaa/analysis/TC_Waves/4postproc/data/"

    # HAFS
    dh = pd.read_csv(dirdata+"Data_HAFS_fromObs.txt", sep='\t')
    hhs = np.array(dh['HTSGW_surface']); hhs[hhs<0.]=np.nan
    hwnd = np.sqrt( np.array(dh['UGRD_surface'])**2 + np.array(dh['VGRD_surface'])**2 )
    hwnd[hwnd<0.]=np.nan; hwnd[hwnd>100.]=np.nan
    del dh

    # GEFS
    for i in range(0,31):
        dg = pd.read_csv(dirdata+"Data_GEFS_fromObs_m"+str(i).zfill(2)+".txt", sep='\t')
        if i==0:
            ghs = np.array([dg['HTSGW_surface']])
            gwnd = np.array([np.sqrt( np.array(dg['UGRD_surface'])**2 + np.array(dg['VGRD_surface'])**2 )])
        else:
            ghs = np.append(ghs, np.array([dg['HTSGW_surface']]), axis=0)
            gwnd = np.append(gwnd, np.array([np.sqrt( np.array(dg['UGRD_surface'])**2 + np.array(dg['VGRD_surface'])**2 )]), axis=0)

    ghs_c = np.copy(ghs[0,:]); ghs_em = np.nanmean(ghs,axis=0)
    gwnd_c = np.copy(gwnd[0,:]); gwnd_em = np.nanmean(gwnd,axis=0)
    ghs_c[ghs_c<0.]=np.nan; ghs_em[ghs_em<0.]=np.nan
    gwnd_c[gwnd_c<0.]=np.nan; gwnd_c[gwnd_c>100.]=np.nan;
    gwnd_em[gwnd_em<0.]=np.nan; gwnd_em[gwnd_em>100.]=np.nan
    del ghs, gwnd

    # PModel
    df = pd.read_csv(dirdata+"Data_Obs_PModel_Default.txt", sep='\t')
    ot = df.time.values
    ot_dt = pd.to_datetime(ot.astype(str), format='%Y%m%d%H%M') # datetime
    ohs = np.array(df['obs_hs']); ownd = np.array(df['obs_wnd'])
    mhs = np.array(df['pm_hs']); mtp = np.array(df['pm_tp'])
    mlat =  np.array(df['lat']); mlon = np.array(df['lon']); mlon[mlon>180.]=mlon[mlon>180.]-360.
    mwnd = np.sqrt( np.array(df['pm_uwnd'])**2 + np.array(df['pm_vwnd'])**2 )
    oid = np.array(df['id']).astype('str')
    cmap = np.array(df['cmap']); csec = np.array(df['csec'])
    ccdist = np.array(df['ccdist']); ccangle = np.array(df['ccangle'])
    cid = np.array(df['cid'])

    # Select quality data, obs and models. Exclude some platforms
    ind = np.where( (ohs>0.2) & (ghs_em>0.2) & (mhs>0.2) & (hhs>0.2) & (cmap==5) &
        (np.char.startswith(oid, 'DWSD')!=True) &
        (np.char.startswith(oid, 'MICROSWIFT')!=True) &
        (np.char.startswith(oid, 'WSRA')!=True) &
        (np.char.startswith(oid, 'SAILDRONE')!=True) )[0]

    ot = ot[ind]; ot_dt = ot_dt[ind]
    mlat = mlat[ind]; mlon = mlon[ind]
    ohs = ohs[ind]; ownd = ownd[ind]
    mhs = mhs[ind]; mtp = mtp[ind]; mwnd = mwnd[ind]
    oid = oid[ind]; cid = cid[ind]
    csec = csec[ind]; ccdist = ccdist[ind]; ccangle = ccangle[ind]
    hhs = hhs[ind]; hwnd = hwnd[ind]
    ghs_c = ghs_c[ind]; ghs_em = ghs_em[ind]
    gwnd_c = gwnd_c[ind]; gwnd_em = gwnd_em[ind]

    del ind

    # Prepare for cross validation: Select events with min of 10 obs/model matchups available.
    fcid=[]
    for i in np.unique(cid):
        ind=np.where((cid==i)&(ohs>1.0))[0]
        if len(ind)>=10:
            fcid=np.append(fcid,int(i))

        del ind

    fcid = np.array(fcid).astype('int') # Total 84 cyclones with valid obs
    ind = np.where(np.isin(cid, fcid)==True)[0]
    ot = ot[ind]; ot_dt = ot_dt[ind]
    mlat = mlat[ind]; mlon = mlon[ind]
    ohs = ohs[ind]; ownd = ownd[ind]
    mhs = mhs[ind]; mtp = mtp[ind]; mwnd = mwnd[ind]
    oid = oid[ind]; cid = cid[ind]
    csec = csec[ind]; ccdist = ccdist[ind]*1000.; ccangle = ccangle[ind]
    hhs = hhs[ind]; hwnd = hwnd[ind]
    ghs_c = ghs_c[ind]; ghs_em = ghs_em[ind]
    gwnd_c = gwnd_c[ind]; gwnd_em = gwnd_em[ind]
    del ind

    # Read IBtracks
    fcnames='/home/ricardo/work/noaa/analysis/TC_Waves/1preproc/cyclonemap/cyclone_names.txt'
    with open(fcnames) as f:
        cnames = [line.strip() for line in f if line.strip()]

    fname='/home/ricardo/work/noaa/analysis/TC_Waves/1preproc/cyclonemap/ibtracs.last4years.list.v04r01.csv'
    wibtr = read_ibtracs(fname)

    # Prepare allocation/organizing
    tlag = np.array([0,3,6,12,24])
    jlday=np.zeros((len(ot_dt)),'f')*np.nan

    clat=np.zeros((len(tlag),len(ot_dt)),'f')*np.nan
    clon=np.zeros((len(tlag),len(ot_dt)),'f')*np.nan
    vmx=np.zeros((len(tlag),len(ot_dt)),'f')*np.nan
    rmw=np.zeros((len(tlag),len(ot_dt)),'f')*np.nan
    vfm=np.zeros((len(tlag),len(ot_dt)),'f')*np.nan

    cdir_sin=np.zeros((len(tlag),len(ot_dt)),'f')*np.nan
    cdir_cos=np.zeros((len(tlag),len(ot_dt)),'f')*np.nan

    r34=np.zeros((len(tlag),len(ot_dt)),'f')*np.nan
    r34_ne=np.zeros((len(tlag),len(ot_dt)),'f')*np.nan
    r34_nw=np.zeros((len(tlag),len(ot_dt)),'f')*np.nan
    r34_se=np.zeros((len(tlag),len(ot_dt)),'f')*np.nan
    r34_sw=np.zeros((len(tlag),len(ot_dt)),'f')*np.nan

    r50=np.zeros((len(tlag),len(ot_dt)),'f')*np.nan
    r50_ne=np.zeros((len(tlag),len(ot_dt)),'f')*np.nan
    r50_nw=np.zeros((len(tlag),len(ot_dt)),'f')*np.nan
    r50_se=np.zeros((len(tlag),len(ot_dt)),'f')*np.nan
    r50_sw=np.zeros((len(tlag),len(ot_dt)),'f')*np.nan

    r64=np.zeros((len(tlag),len(ot_dt)),'f')*np.nan
    r64_ne=np.zeros((len(tlag),len(ot_dt)),'f')*np.nan
    r64_nw=np.zeros((len(tlag),len(ot_dt)),'f')*np.nan
    r64_se=np.zeros((len(tlag),len(ot_dt)),'f')*np.nan
    r64_sw=np.zeros((len(tlag),len(ot_dt)),'f')*np.nan

    # Include variables with time lag
    for i in range(0,len(ot_dt)):
        jlday[i] = int(pd.Timestamp(ot_dt[i]).dayofyear)
        for j in range(0,len(tlag)):
            ind=np.where( (wibtr['itime']==(ot_dt[i]-pd.Timedelta(hours=tlag[j])).round('3h')) & (wibtr['iname']==cnames[cid[i]]) )
            if size(ind)>0:
                ind=ind[0][0]
                
                clat[j,i] = wibtr['ilat'][ind]; clon[j,i] = wibtr['ilon'][ind]
                vmx[j,i] = wibtr['iVmax'][ind]; rmw[j,i] = wibtr['iRmax'][ind]
                vfm[j,i] = wibtr['iVfm'][ind]
                dir_rad = np.deg2rad(wibtr['iDir'][ind])
                cdir_sin[j,i] = np.sin(dir_rad); cdir_cos[j,i] = np.cos(dir_rad)

                r34[j,i] = wibtr['iR34'][ind]
                r34_ne[j,i] = wibtr['iR34_NE'][ind]
                r34_nw[j,i] = wibtr['iR34_NW'][ind]
                r34_se[j,i] = wibtr['iR34_SE'][ind]
                r34_sw[j,i] = wibtr['iR34_SW'][ind]

                r50[j,i] = wibtr['iR50'][ind]
                r50_ne[j,i] = wibtr['iR50_NE'][ind]
                r50_nw[j,i] = wibtr['iR50_NW'][ind]
                r50_se[j,i] = wibtr['iR50_SE'][ind]
                r50_sw[j,i] = wibtr['iR50_SW'][ind]

                r64[j,i] = wibtr['iR64'][ind]
                r64_ne[j,i] = wibtr['iR64_NE'][ind]
                r64_nw[j,i] = wibtr['iR64_NW'][ind]
                r64_se[j,i] = wibtr['iR64_SE'][ind]
                r64_sw[j,i] = wibtr['iR64_SW'][ind]

                del ind, dir_rad

        print(repr(i))

    # Save final dataset. 32582 obs/model matchups. 131 variables (129 for the ML).
    ncfile = nc.Dataset('TC_ModelObs_CyclInfo.nc', "w", format=fnetcdf) 
    ncfile.history='TC dataset prepared for ML modeling. QC Obs and Model (PModel, GEFS, and HAFS), plus cyclone info from IBTracks.'
    ncfile.info='tlag is the time lag (hours, backwards) to select cyclone data.'
    # create  dimensions.
    ncfile.createDimension('time', len(ot))
    ncfile.createDimension('tlag', len(tlag))
    # time arrays
    vot = ncfile.createVariable('time',np.dtype('float64'),('time'))
    vjlday = ncfile.createVariable('jlday',np.dtype('int32'),('time'))
    vtlag = ncfile.createVariable('tlag',np.dtype('int32'),('tlag'))
    # Obs and Model data
    void = ncfile.createVariable('obs_ID',np.dtype('a25'),('time'))
    vcid = ncfile.createVariable('cid',np.dtype('int32'),('time'))
    vmlat = ncfile.createVariable('mlat',np.dtype('float32'),('time'))
    vmlon = ncfile.createVariable('mlon',np.dtype('float32'),('time'))
    vcsec = ncfile.createVariable('csec',np.dtype('int32'),('time'))
    vccdist = ncfile.createVariable('ccdist',np.dtype('float32'),('time'))
    vccangle = ncfile.createVariable('ccangle',np.dtype('float32'),('time'))
    vohs= ncfile.createVariable('ohs',np.dtype('float32'),('time'))
    vownd= ncfile.createVariable('ownd',np.dtype('float32'),('time'))
    # PModel
    vmhs= ncfile.createVariable('mhs',np.dtype('float32'),('time'))
    vmtp= ncfile.createVariable('mtp',np.dtype('float32'),('time'))
    vmwnd= ncfile.createVariable('mwnd',np.dtype('float32'),('time'))
    # HAFS
    vhhs = ncfile.createVariable('hhs',np.dtype('float32'),('time'))
    vhwnd = ncfile.createVariable('hwnd',np.dtype('float32'),('time'))
    # GEFS
    vghs_c = ncfile.createVariable('ghs_c',np.dtype('float32'),('time'))
    vghs_em = ncfile.createVariable('ghs_em',np.dtype('float32'),('time'))
    vgwnd_c = ncfile.createVariable('gwnd_c',np.dtype('float32'),('time'))
    vgwnd_em = ncfile.createVariable('gwnd_em',np.dtype('float32'),('time'))
    # Cyclone Info
    vclat = ncfile.createVariable('clat',np.dtype('float32'),('tlag','time'))
    vclon = ncfile.createVariable('clon',np.dtype('float32'),('tlag','time'))
    vvmx = ncfile.createVariable('vmx',np.dtype('float32'),('tlag','time'))
    vrmw = ncfile.createVariable('rmw',np.dtype('float32'),('tlag','time'))
    vvfm = ncfile.createVariable('vfm',np.dtype('float32'),('tlag','time'))
    vcdir_sin = ncfile.createVariable('cdir_sin',np.dtype('float32'),('tlag','time'))
    vcdir_cos = ncfile.createVariable('cdir_cos',np.dtype('float32'),('tlag','time'))
    vr34 = ncfile.createVariable('r34',np.dtype('float32'),('tlag','time'))
    vr34_ne = ncfile.createVariable('r34_ne',np.dtype('float32'),('tlag','time'))
    vr34_nw = ncfile.createVariable('r34_nw',np.dtype('float32'),('tlag','time'))
    vr34_se = ncfile.createVariable('r34_se',np.dtype('float32'),('tlag','time'))
    vr34_sw = ncfile.createVariable('r34_sw',np.dtype('float32'),('tlag','time'))
    vr50 = ncfile.createVariable('r50',np.dtype('float32'),('tlag','time'))
    vr50_ne = ncfile.createVariable('r50_ne',np.dtype('float32'),('tlag','time'))
    vr50_nw = ncfile.createVariable('r50_nw',np.dtype('float32'),('tlag','time'))
    vr50_se = ncfile.createVariable('r50_se',np.dtype('float32'),('tlag','time'))
    vr50_sw = ncfile.createVariable('r50_sw',np.dtype('float32'),('tlag','time'))
    vr64 = ncfile.createVariable('r64',np.dtype('float32'),('tlag','time'))   
    vr64_ne = ncfile.createVariable('r64_ne',np.dtype('float32'),('tlag','time')) 
    vr64_nw = ncfile.createVariable('r64_nw',np.dtype('float32'),('tlag','time')) 
    vr64_se = ncfile.createVariable('r64_se',np.dtype('float32'),('tlag','time')) 
    vr64_sw = ncfile.createVariable('r64_sw',np.dtype('float32'),('tlag','time')) 
    # Units
    vot.units = "seconds since 1970-01-01 00:00:00.0 0:00"
    vmlat.units = 'degrees_north' ; vmlon.units = 'degrees_east'
    vclat.units = 'degrees_north' ; vclon.units = 'degrees_east'
    vr34.units = "m"; vvmx.units = "m/s"
    vohs.units = "m"; vownd.units = "m/s"
    vmtp.units = 's'; vmwnd.units = 'm/s'
    vcsec.units = '1'; vccdist.units = 'm'; vccangle.units = 'degree'
    vhhs.units = 'm'; vhwnd.units = 'm/s'
    vghs_c.units = 'm'; vghs_em.units = 'm'
    vgwnd_c.units = 'm/s'; vgwnd_em.units = 'm/s'
    vrmw.units = 'm'; vvfm.units = 'm/s'
    vcdir_sin.units = '1'; vcdir_cos.units = '1'
    vr34_ne.units = 'm'; vr34_nw.units = 'm'; vr34_se.units = 'm'; vr34_sw.units = 'm'
    vr50.units = 'm'; vr50_ne.units = 'm'; vr50_nw.units = 'm'; vr50_se.units = 'm'; vr50_sw.units = 'm'
    vr64.units = 'm'; vr64_ne.units = 'm'; vr64_nw.units = 'm'; vr64_se.units = 'm'; vr64_sw.units = 'm'
    vjlday.units = '1'; vtlag.units = 'hours'
    # Allocate Data
    vot[:] = np.array(ot_to_epoch_seconds(ot)).astype('float64')[:]; vjlday[:] = jlday[:]; vtlag[:] = tlag[:]
    # Obs and Model data
    void[:] = oid[:]; vcid[:] = cid[:]
    vmlat[:] = mlat[:]; vmlon[:] = mlon[:]
    vcsec[:] = csec[:]; vccdist[:] = ccdist[:]; vccangle[:] = ccangle[:]
    vohs[:] = ohs[:]; vownd[:] = ownd[:]
    vmhs[:] = mhs[:]; vmtp[:] = mtp[:]; vmwnd[:] = mwnd[:]
    vhhs[:] = hhs[:]; vhwnd[:] = hwnd[:]
    vghs_c[:] = ghs_c[:]; vghs_em[:] = ghs_em[:]
    vgwnd_c[:] = gwnd_c[:]; vgwnd_em[:] = gwnd_em[:]
    # Cyclone Info
    vclat[:,:] = clat[:,:]; vclon[:,:] = clon[:,:]
    vvmx[:,:] = vmx[:,:]; vrmw[:,:] = rmw[:,:]; vvfm[:,:] = vfm[:,:]
    vcdir_sin[:,:] = cdir_sin[:,:]; vcdir_cos[:,:] = cdir_cos[:,:]
    vr34[:,:] = r34[:,:]; vr34_ne[:,:] = r34_ne[:,:]; vr34_nw[:,:] = r34_nw[:,:]; vr34_se[:,:] = r34_se[:,:]; vr34_sw[:,:] = r34_sw[:,:]
    vr50[:,:] = r50[:,:]; vr50_ne[:,:] = r50_ne[:,:]; vr50_nw[:,:] = r50_nw[:,:]; vr50_se[:,:] = r50_se[:,:]; vr50_sw[:,:] = r50_sw[:,:]
    vr64[:,:] = r64[:,:]; vr64_ne[:,:] = r64_ne[:,:]; vr64_nw[:,:] = r64_nw[:,:]; vr64_se[:,:] = r64_se[:,:]; vr64_sw[:,:] = r64_sw[:,:]
    ncfile.close()
    print('netcdf ok ')

    # List of TC names
    selected_names = np.array(cnames)[np.array(cid, dtype=int)]
    np.savetxt("selected_cyclones.txt", np.unique(selected_names), fmt="%s")

