# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 18:34:51 2024

@author: GnegyE
"""

import rdata
import xarray as xr
import glob, os
import numpy as np
import pandas as pd
from netCDF4 import Dataset
import datetime
import pymannkendall as mk


#%% read from RData file

file = '/Users/GnegyE/Desktop/CanKrig_189912_201909.RData'
parsed = rdata.parser.parse_file(file)
converted = rdata.conversion.convert(parsed)

lats = converted['lats'].values
lons = converted['lons'].values
times = converted['yms']
precip = converted['CanKrig']

times_dt = pd.to_datetime(times, format='%Y%m')

precip_da = xr.DataArray(precip,
                  dims=('lats', 'lons','time'),
                  coords={'lats':  (('lats', 'lons'), lats), 'lons':  (('lats', 'lons'), lons), 'time': times_dt})

          
#%% select time period of interest and base period

#inclusive years
startyear = 1948
endyear = 2018

startyear_base = 1961
endyear_base = 1990

#%% annual 
precip_select = precip_da.sel(time=slice(str(startyear),str(endyear)))
precip_base = precip_da.sel(time=slice(str(startyear_base),str(endyear_base)))

#TODO: fix sum of NaNs being 0
precip_yr = precip_select.resample(time="AS").sum() #annual sums

#same but for base period
precip_yr_base = precip_base.resample(time="AS").sum()
precip_mean_base = precip_yr_base.mean(dim='time')

#puts cankrig data in same format as cangrd
precip_wrt_base = (precip_yr - precip_mean_base)/precip_mean_base

#%% same but for seasonal 

start_date = pd.Timestamp(str(startyear-1)+'-12-01') #get dec of previous year for first DJF
end_date = pd.Timestamp(str(endyear)+'-11-01') #stop at nov of final year 

start_date_base = pd.Timestamp(str(startyear_base-1)+'-12-01')
end_date_base = pd.Timestamp(str(endyear_base)+'-11-01')

precip_select_seas = precip_da.sel(time=slice(start_date,end_date))
precip_base_seas = precip_da.sel(time=slice(start_date_base,end_date_base))

#TODO: fix sum of NaNs being 0
precip_seas = precip_select_seas.resample(time="QS-DEC").sum()
precip_seas_base = precip_base_seas.resample(time="QS-DEC").sum()

precip_seas_gr = precip_seas.groupby('time.season')
precip_seas_means_base = precip_seas_base.groupby('time.season').mean('time')

precip_wrt_base_DJF = (precip_seas_gr['DJF']-precip_seas_means_base.sel(season='DJF'))/precip_seas_means_base.sel(season='DJF')
precip_wrt_base_MAM = (precip_seas_gr['MAM']-precip_seas_means_base.sel(season='MAM'))/precip_seas_means_base.sel(season='MAM')
precip_wrt_base_JJA = (precip_seas_gr['JJA']-precip_seas_means_base.sel(season='JJA'))/precip_seas_means_base.sel(season='JJA')
precip_wrt_base_SON = (precip_seas_gr['SON']-precip_seas_means_base.sel(season='SON'))/precip_seas_means_base.sel(season='SON')

#%% read in cangrd data

points_file = #path to "CANGRD_points_LL.csv" file on shared drive 
data_filepath = #path to "precip_anomalies" fodler on shared drive 
  
x_len = 95
y_len = 125

df = pd.read_csv(points_file, sep=",", header=None, names=["id1", "id2", "lat", "lon"])

lats_cangrd = df['lat'].values.reshape((x_len, y_len))
lons_cangrd = df['lon'].values.reshape((x_len, y_len))

#lil function to read cangrid data 
def read_cangrd(file):
    #file in format YYYYDD (DD is month or 13 for ann, 14-17 for seas - see cangrd readme)
    #open the GRD file
    df = pd.read_csv(data_filepath + 'p' + file + '.grd', header=None,skiprows=[0,1,2,3,4],names=['pr'])
    #fix spacing issues 
    df = df['pr'].str.split(expand=True).stack().values
    data = df.reshape((x_len, y_len)).astype(float)
    data[data > 10e30] = np.nan
    return data

#inclusive years
startyear_cangrd = 1948
endyear_cangrd = 2012

years_cangrd = [datetime.datetime(year, 1, 1) for year in range(startyear_cangrd, endyear_cangrd + 1)]

def load_cangrd(seas_ID):
    cangrd_data = np.empty((x_len, y_len, len(years_cangrd)))
    
    for i in range(len(years_cangrd)):
        year = years_cangrd[i].year
        print(year)
        file = str(year) + seas_ID #13 for "annual" - see readme for cangrd
        cangrd_data[:, :, i] = read_cangrd(file) 

    cangrd_da = xr.DataArray(cangrd_data,
                      dims=('lats', 'lons','time'),
                      coords={'lats':  (('lats', 'lons'), lats_cangrd), 'lons':  (('lats', 'lons'), lons_cangrd), 'time': years_cangrd})


    return(cangrd_da)

cangrd_da_ann = load_cangrd("13")
cangrd_da_djf = load_cangrd("14")
cangrd_da_mam = load_cangrd("15")
cangrd_da_jja = load_cangrd("16")
cangrd_da_son = load_cangrd("17")


#%%

def get_trendline(pr):
    mk_samp = mk.original_test(pr.values)
    
    # for absolute trends
    #x = np.linspace(0,len(pr.time),len(pr.time)+1)
    #trendline = mk_samp.slope * x + mk_samp.intercept
    #perc_change = (trendline[-1] - trendline[0])/(abs(trendline[0])) * 100
    
    # for anomolies
    perc_change = mk_samp.slope * len(pr.time) * 100  
    
    return perc_change  

#%%

# pr_da should already be in yearly (seasonal) format (lat,lon,years)
def get_trends(pr_da):
    
    #make empty array to fill
    perc_change = np.empty(pr_da.lats.shape)
        
    for i in range(pr_da.lats.shape[0]):
        print(str(i) + "/" + str(pr_da.lats.shape[0]))
        
        for j in range(pr_da.lats.shape[1]):
            
            #all timesteps (of selected range) at one lat/lon point
            precip_i = pr_da.isel(lats=i,lons=j)
            
            #if its out of the canada terrestrial boundary
            if np.all(np.isnan(precip_i)) or np.all(precip_i==0):
                perc_change[i,j] = np.nan

            else:
                perc_change[i, j] = get_trendline(precip_i)
                          
    return(perc_change)

#%% anomolies
cankrig_perc_change_ann = get_trends(precip_wrt_base)
cankrig_perc_change_DJF = get_trends(precip_wrt_base_DJF)
cankrig_perc_change_MAM = get_trends(precip_wrt_base_MAM)
cankrig_perc_change_JJA = get_trends(precip_wrt_base_JJA)
cankrig_perc_change_SON = get_trends(precip_wrt_base_SON)

#%% raw precip

#cankrig_perc_change_ann = get_trends(precip_yr)
cankrig_perc_change_DJF = get_trends(precip_seas_gr['DJF'])
cankrig_perc_change_MAM = get_trends(precip_seas_gr['MAM'])
cankrig_perc_change_JJA = get_trends(precip_seas_gr['JJA'])
cankrig_perc_change_SON = get_trends(precip_seas_gr['SON'])


#%% cangrd

cangrd_perc_change_ann = get_trends(cangrd_da_ann)
cangrd_perc_change_djf = get_trends(cangrd_da_djf)
cangrd_perc_change_mam = get_trends(cangrd_da_mam)
cangrd_perc_change_jja = get_trends(cangrd_da_jja)
cangrd_perc_change_son = get_trends(cangrd_da_son)


#%% save it (cankrig) as a netcdf file

#file = "C:/Users/GnegyE/Desktop/precip_trends/precip_trends_" + str(startyear) + "-" + str(endyear) + "_cankrig_anomolies.nc"
file = "C:/Users/GnegyE/Desktop/precip_trends/precip_trends_" + str(startyear) + "-" + str(endyear) + "_cankrig_abs.nc"

with Dataset(file,'w') as nc_file:

    lat_dim = nc_file.createDimension('lat',np.shape(lats)[0])
    lon_dim = nc_file.createDimension('lon',np.shape(lons)[1])

    ann_trend_var = nc_file.createVariable('ann_precip_trend','f8',('lat','lon'))
    djf_trend_var = nc_file.createVariable('djf_precip_trend','f8',('lat','lon'))
    mam_trend_var = nc_file.createVariable('mam_precip_trend','f8',('lat','lon'))
    jja_trend_var = nc_file.createVariable('jja_precip_trend','f8',('lat','lon'))
    son_trend_var = nc_file.createVariable('son_precip_trend','f8',('lat','lon'))

    lat_var = nc_file.createVariable('lat','f8',('lat','lon'))
    lon_var = nc_file.createVariable('lon','f8',('lat','lon'))
    
    ann_trend_var.Description = str(startyear) + "-" + str(endyear) +' CanKrig annual precipitation trend'
    djf_trend_var.Description = str(startyear) + "-" + str(endyear) +' CanKrig DJF precipitation trend'
    mam_trend_var.Description = str(startyear) + "-" + str(endyear) +' CanKrig MAM precipitation trend'
    jja_trend_var.Description = str(startyear) + "-" + str(endyear) +' CanKrig JJA precipitation trend'
    son_trend_var.Description = str(startyear) + "-" + str(endyear) +' CanKrig SON precipitation trend'

    ann_trend_var[:,:] = cankrig_perc_change_ann[:,:]
    djf_trend_var[:,:] = cankrig_perc_change_DJF[:,:]
    mam_trend_var[:,:] = cankrig_perc_change_MAM[:,:]
    jja_trend_var[:,:] = cankrig_perc_change_JJA[:,:]
    son_trend_var[:,:] = cankrig_perc_change_SON[:,:]
    
    lat_var[:,:] = lats[:,:]
    lon_var[:,:] = lons[:,:]

#%% save cangrd as a netcdf file

file = "C:/Users/GnegyE/Desktop/precip_trends/precip_trends_" + str(startyear_cangrd) + "-" + str(endyear_cangrd) + "_cangrd_anomolies.nc"

with Dataset(file,'w') as nc_file:

    lat_dim = nc_file.createDimension('lat',np.shape(lats_cangrd)[0])
    lon_dim = nc_file.createDimension('lon',np.shape(lons_cangrd)[1])

    ann_trend_var = nc_file.createVariable('ann_precip_trend','f8',('lat','lon'))
    djf_trend_var = nc_file.createVariable('djf_precip_trend','f8',('lat','lon'))
    mam_trend_var = nc_file.createVariable('mam_precip_trend','f8',('lat','lon'))
    jja_trend_var = nc_file.createVariable('jja_precip_trend','f8',('lat','lon'))
    son_trend_var = nc_file.createVariable('son_precip_trend','f8',('lat','lon'))

    lat_var = nc_file.createVariable('lat','f8',('lat','lon'))
    lon_var = nc_file.createVariable('lon','f8',('lat','lon'))
    
    ann_trend_var.Description = str(startyear) + "-" + str(endyear) +' CanGRD annual precipitation trend'
    djf_trend_var.Description = str(startyear) + "-" + str(endyear) +' CanGRD DJF precipitation trend'
    mam_trend_var.Description = str(startyear) + "-" + str(endyear) +' CanGRD MAM precipitation trend'
    jja_trend_var.Description = str(startyear) + "-" + str(endyear) +' CanGRD JJA precipitation trend'
    son_trend_var.Description = str(startyear) + "-" + str(endyear) +' CanGRD SON precipitation trend'

    ann_trend_var[:,:] = cangrd_perc_change_ann[:,:]
    djf_trend_var[:,:] = cangrd_perc_change_djf[:,:]
    mam_trend_var[:,:] = cangrd_perc_change_mam[:,:]
    jja_trend_var[:,:] = cangrd_perc_change_jja[:,:]
    son_trend_var[:,:] = cangrd_perc_change_son[:,:]
    
    lat_var[:,:] = lats_cangrd[:,:]
    lon_var[:,:] = lons_cangrd[:,:]
