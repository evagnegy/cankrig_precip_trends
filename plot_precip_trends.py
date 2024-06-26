# -*- coding: utf-8 -*-
"""
Created on Sat May 11 09:12:27 2024

@author: GnegyE
"""


import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from netCDF4 import Dataset
import cartopy.feature as cf
import matplotlib as mpl
import matplotlib.colors as pltcol
import pylab as pl
import matplotlib.colorbar as cb

#%% read from netcdf files

#file = 'C:/Users/GnegyE/Desktop/precip_trends/nc_files/precip_trends_1948-2018_cankrig_anomolies.nc'
#file = 'C:/Users/GnegyE/Desktop/precip_trends/nc_files/precip_trends_1948-2018_cankrig_abs.nc'
file = 'C:/Users/GnegyE/Desktop/precip_trends/nc_files/precip_trends_1948-2012_cangrd_anomolies.nc'

#update based on cankrig or cangrd
startyear = 1948
endyear=2012 #cangrd=2012, cankrig=2018


nc = Dataset(file,'r')

ann_precip_trend = nc.variables['ann_precip_trend']
djf_precip_trend = nc.variables['djf_precip_trend']
mam_precip_trend = nc.variables['mam_precip_trend']
jja_precip_trend = nc.variables['jja_precip_trend']
son_precip_trend = nc.variables['son_precip_trend']

lats = nc.variables['lat']
lons = nc.variables['lon']

#%%

#made colorbar from image on this page: 
#https://www.researchgate.net/figure/Change-in-average-precipitation-based-on-multi-model-mean-projections-for-2081-2100_fig1_342989239    
colors_pr = ['#a05323', '#b46b28', '#ce853e' , '#f6cd84', '#f4e19e', '#fdf7ba', \
             '#d0e8c6', '#a3d4aa', '#50bfa0', '#4b86a3', '#255e7f', '#123859']


#how many colors are wanted if using step colorbar (as opposed to gradient, see N below)
#the -2 is if you want arrows at the min/max with different colors
color_count = len(colors_pr)-2

#optional but cleans up the labels
lim = [-.75,.75] #data limits (unit not in %)

ticks = np.linspace(lim[0],lim[1],color_count+1) #location where I want ticks
ticklabels = [str(round(i*100)) for i in (ticks)] #makes the values be a % rather than decimal


# makes the colorbar
# NOTE: if you wanted the colorbar to be a gradient, you would change N to be really high (e.g. 1000) and remove [1:-1] indexing
cmap_pr = pltcol.LinearSegmentedColormap.from_list("custom", colors_pr[1:-1],N=color_count) #[1:-1] is to not include the min/max arrow colors
#uncomment this line to see a gradient instead of step
cmap_pr.set_over(colors_pr[-1]) #add the max arrow color
cmap_pr.set_under(colors_pr[0]) #add the min arrow color

# this is a function that plots a colorbar, where the input is a colormap (cmap)
def plot_cbar(cmap_pr,lim,label,norm=1,ticklabels=1,ticks=1):
    
    pl.figure(figsize=(20, 1.5),dpi=250)
    pl.gca().set_visible(False)
    cax = pl.axes([0.1, 0.2, 0.8, 0.6])
    
    if norm==1:
        norm = pltcol.Normalize(vmin=lim[0], vmax=lim[1])
    if isinstance(ticks,int):
        ticks=None
    else:
        ticks=ticks
    cb.ColorbarBase(cax,cmap=cmap_pr,norm=norm,orientation="horizontal", extend='both',ticks=ticks)
    
    if ticklabels!=1:
        cax.set_xticks(ticks,ticklabels)
    
    cax.xaxis.set_tick_params(size=12,width=4)
    pl.xticks(fontsize=38)
    pl.xlabel(label,fontsize=38)

    plt.show()

#see if you like the colorbar
plot_cbar(cmap_pr,lim,'Precip Change (%)',ticks=ticks,ticklabels=ticklabels)

#%%
#plot the data 


cmap=cmap_pr
vmin=lim[0]*100
vmax=lim[1]*100
extend='both'

def make_perc_plots(data, seas):
    fig = plt.figure(figsize=(10,10),dpi=200)
    ax = fig.add_subplot(1,1,1, projection=ccrs.RotatedPole(pole_latitude=42.5,pole_longitude=83))
    
    
    plt.pcolormesh(lons,lats,data,transform=ccrs.PlateCarree(),cmap=cmap,vmin=vmin,vmax=vmax)
    #plt.scatter(lons,lats,c=precip,s=0.3,transform=ccrs.PlateCarree(),cmap=cmap,vmin=vmin,vmax=vmax)
    
    
    states_provinces = cf.NaturalEarthFeature(category='cultural',name='admin_1_states_provinces_lines',scale='50m',facecolor='none')
    
    ax.coastlines(resolution='50m')
    ax.add_feature(cf.BORDERS)
    ax.add_feature(states_provinces)
    
    ax.set_extent([-135,-55,40,85],crs=ccrs.PlateCarree())
    #ax.set_extent([-142,-133,55,65],crs=ccrs.PlateCarree())
    
    cbar_ax = fig.add_axes([0.2,0.15,0.62,0.02])
    
    fig.colorbar(mpl.cm.ScalarMappable(cmap=cmap, norm=mpl.colors.Normalize(vmin=vmin, vmax=vmax)),cax=cbar_ax,orientation='horizontal',extend=extend,ticks=ticks*100)
    cbar_ax.tick_params(labelsize=16)
    cbar_ax.set_xlabel('Percentage Change (%)',size=18)
    
    fig.suptitle('CanGRD ' + str(startyear) +'-' + str(endyear) + ": " + seas + "\n(relative to 1961-1990 mean)",fontsize=20,y=0.88)
    
    output_dir = '/Users/GnegyE/Desktop/precip_trends/figures/'
    plt.savefig(output_dir + seas + '_cangrd_anomolies.png',bbox_inches='tight')

make_perc_plots(ann_precip_trend, "Annual")
make_perc_plots(djf_precip_trend, "DJF")
make_perc_plots(mam_precip_trend, "MAM")
make_perc_plots(jja_precip_trend, "JJA")
make_perc_plots(son_precip_trend, "SON")
