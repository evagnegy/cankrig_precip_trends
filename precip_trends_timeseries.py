# -*- coding: utf-8 -*-
"""
Created on Sat May 11 09:12:27 2024

@author: GnegyE
"""

import rdata
import xarray as xr
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import datetime



#%% read from RData file

file = '/Users/GnegyE/Desktop/CanKrig_189912_201909.RData'
parsed = rdata.parser.parse_file(file)
converted = rdata.conversion.convert(parsed)

lats = converted['lats'].values
lons = converted['lons'].values
time = converted['yms']
precip = converted['CanKrig']

times = pd.to_datetime(time, format='%Y%m')

precip_da = xr.DataArray(precip,
                  dims=('lats', 'lons','time'),
                  coords={'lats':  (('lats', 'lons'), lats), 'lons':  (('lats', 'lons'), lons), 'time': times})

          
#%% select time period of interest and base period

#inclusive years
startyear = 1948
endyear = 2018

startyear_base = 1961
endyear_base = 1990

precip_select = precip_da.sel(time=slice(str(startyear),str(endyear)))
precip_base = precip_da.sel(time=slice(str(startyear_base),str(endyear_base)))

#%% annual spatial averages
precip_yr = precip_select.resample(time="AS").sum() #annual sums
area_weights = np.cos(np.deg2rad(precip_da['lats']))
precip_yr_weighted = precip_yr.weighted(area_weights).mean(dim=("lons","lats")) #weighted avg

#same but for base period
precip_yr_base = precip_base.resample(time="AS").sum()
precip_yr_weighted_base = precip_yr_base.weighted(area_weights).mean(dim=("lons","lats"))
precip_avg_weighted_base = precip_yr_weighted_base.mean() #average 30 year period

#puts cankrig data in same format as cangrd
precip_wrt_base = (precip_yr_weighted - precip_avg_weighted_base)/precip_avg_weighted_base

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

cangrd_data = np.empty((x_len, y_len, len(years_cangrd)))

for i in range(len(years_cangrd)):
    year = years_cangrd[i].year
    print(year)
    file = str(year) + '13' #13 for "annual" - see readme for cangrd
    cangrd_data[:, :, i] = read_cangrd(file) 

#%%

cangrd_da = xr.DataArray(cangrd_data,
                  dims=('lats', 'lons','time'),
                  coords={'lats':  (('lats', 'lons'), lats_cangrd), 'lons':  (('lats', 'lons'), lons_cangrd), 'time': years_cangrd})


#%% annual spatial averages - cangrd
area_weights = np.cos(np.deg2rad(cangrd_da['lats']))
precip_wrt_base_cangrd = cangrd_da.weighted(area_weights).mean(dim=("lons","lats")) #get weighted avg

#%%

plt.plot(precip_wrt_base_cangrd.time,precip_wrt_base_cangrd*100,label='CanGRD')
plt.plot(precip_wrt_base.time,precip_wrt_base*100,label='CanKrig')

plt.legend()
plt.ylabel('Percent, %')
plt.title('Annual precipitation anomolies, Canada')
