# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 20:27:48 2022

@author: Lauren.Dares
"""

import pandas as pd
import numpy as np
import os
from osgeo import gdal
import geopandas as geo
from glob import glob
import rioxarray as rx
import rasterstats as rs
import rasterio as rio
from rasterio.enums import Resampling
from datetime import datetime as dt

def read_file(file):
    with rio.open(file) as src:
        return(src.read(1))

def select_months(file, ch1, ch2, dates):
    firstday = file[ch1:ch2]
    if firstday in dates:
        return file
    
def make_null(array, value):
    array = array.astype('float')
    array[array == value] = np.NAN
    return array

#%% yearly means for sst (summer)
def sst_means(year):
    var='Sea Surface Temperature'
    folder='Mapped/Rasters/'
    data_dir = 'OB.DAAC Data/{var}/{folder}'.format(var=var, folder=folder) # Or sys.argv[1]
    file_list = glob(os.path.join(data_dir, '*.tif'))
    
    smr_md = ['0501', '0601', '0701', '0801', '0901']
    dates = [str(year)+x for x in smr_md]
    smr_list = [select_months(i, 142, 150, dates) for i in file_list]
    smr_list1 = list(filter(None.__ne__, smr_list))

    # Read all data as a list of numpy arrays 
    array_list = [read_file(x) for x in smr_list1]
    array_list1 = [make_null(x, -32767) for x in array_list]

    # Perform averaging
    array_out = np.nanmean(array_list1, axis=0)
    array_out1 = array_out*0.004999999888241291

    # Get metadata from one of the input files
    with rio.open(smr_list1[0]) as src:
        meta = src.meta
    
    meta.update(dtype=rio.float32)
    
    # Write output file
    with rio.open('OB.DAAC Data/Sea Surface Temperature/Mapped/Rasters/Mean_SST/{var}_smr_{yr}_mean1.tif'.format(var=var, yr=year), 'w', **meta) as dst:
        dst.write(array_out1.astype(rio.float32), 1)

[sst_means(year) for year in range(2002, 2023)]

#%% yearly means for sst (winter)
def means(year, var, scale_factor, season, ch_start, ch_end):
    folder='Mapped/Rasters/'
    data_dir = 'OB.DAAC Data/{var}/{folder}'.format(var=var, folder=folder) # Or sys.argv[1]
    file_list = glob(os.path.join(data_dir, '*.tif'))
    
    if season == 'win':
        dates = [str(year)+'1001', str(year)+'1101',str(year)+'1201',
             str(year+1)+'0101', str(year+1)+'0201', str(year+1)+'0301', str(year+1)+'0401']
    else:
        smr_md = ['0501', '0601', '0701', '0801', '0901']
        dates = [str(year)+x for x in smr_md]
    smr_list = [select_months(i, ch_start, ch_end, dates) for i in file_list]
    smr_list1 = list(filter(None.__ne__, smr_list))

    # Read all data as a list of numpy arrays 
    array_list = [read_file(x) for x in smr_list1]
    array_list1 = [make_null(x, -32767) for x in array_list]

    # Perform averaging
    array_out = np.nanmean(array_list1, axis=0)
    array_out1 = array_out*scale_factor

    # Get metadata from one of the input files
    with rio.open(smr_list1[0]) as src:
        meta = src.meta
    
    meta.update(dtype=rio.float32)
    
    # Write output file
    with rio.open('OB.DAAC Data/Sea Surface Temperature/Mapped/Rasters/Mean_SST/{var}_{season}_{yr}_mean2.tif'.format(var=var, yr=year, season=season), 'w', **meta) as dst:
        dst.write(array_out1.astype(rio.float32), 1)

#%%
sst_scale_factor = 0.004999999888241291
[means(year, 'Sea Surface Temperature', sst_scale_factor, 'win', 143, 151) for year in range(2002, 2023)]

cha_scale_factor = 1
[means(year, 'Chlorophyll-A', cha_scale_factor, 'win', 136, 144) for year in range(2002, 2023)]
[means(year, 'Chlorophyll-A', cha_scale_factor, 'smr', 133, 141) for year in range(2002, 2023)]
#%%
[print(year) for year in range(2002, 2023)]
#%%SST Summer for entire study period
years = range(2002, 2023)
smr_md = ['0501', '0601', '0701', '0801', '0901']

data_dir = 'OB.DAAC Data/{var}/{folder}'.format(var=var, folder=folder) # Or sys.argv[1]
file_list = glob(os.path.join(data_dir, '*.tif'))

dates = [str(y)+x for y in years for x in smr_md]
smr_list = [select_months(i, 143, 151, dates) for i in file_list]
smr_list1 = list(filter(None.__ne__, smr_list))

array_list = [read_file(x) for x in smr_list1]
array_list1 = [make_null(x, -32767) for x in array_list]

# Perform averaging
array_out = np.nanmean(array_list1, axis=0)
array_out1 = array_out*0.004999999888241291

# Get metadata from one of the input files
with rio.open(smr_list1[0]) as src:
    meta = src.meta

meta.update(dtype=rio.float32)

# Write output file
with rio.open('OB.DAAC Data/Sea Surface Temperature/Mapped/Rasters/Mean_SST/{var}_smr_allyrs_mean1.tif'.format(var=var), 'w', **meta) as dst:
    dst.write(array_out1.astype(rio.float32), 1)

#%%

var='Sea Surface Temperature'
folder='/Mapped/Rasters/Mean_SST'

data_dir = 'OB.DAAC Data/{var}{folder}'.format(var=var, folder=folder) # Or sys.argv[1]
file_list = glob(os.path.join(data_dir, '*.tif'))
        
dates = ['20220501', '20220601', '20220701', '20220801', '20220901']
smr_list = [select_months(i, 143, 151, dates) for i in file_list]
smr_list1 = list(filter(None.__ne__, smr_list))

# Read all data as a list of numpy arrays 
array_list = [read_file(x) for x in smr_list1]
array_list1 = [make_null(x, -32767) for x in array_list]

# Perform averaging
array_out = np.nanmean(array_list1, axis=0)
array_out1 = array_out*0.004999999888241291

# Get metadata from one of the input files
with rio.open(smr_list1[0]) as src:
    meta = src.meta

meta.update(dtype=rio.float32)

# Write output file
with rio.open('OB.DAAC Data/Sea Surface Temperature/Mapped/Rasters/Mean_SST/{var}_smr_{yr}_mean.tif'.format(var=var, yr=year), 'w', **meta) as dst:
    dst.write(array_out1.astype(rio.float32), 1)

#%%SST Summer for entire study period
years = range(2002, 2023)
smr_md = ['0501', '0601', '0701', '0801', '0901']

data_dir = 'OB.DAAC Data/{var}/{folder}'.format(var=var, folder=folder) # Or sys.argv[1]
file_list = glob(os.path.join(data_dir, '*.tif'))

dates = [str(y)+x for y in years for x in smr_md]
smr_list = [select_months(i, 143, 151, dates) for i in file_list]
smr_list1 = list(filter(None.__ne__, smr_list))

array_list = [read_file(x) for x in smr_list1]
array_list1 = [make_null(x, -32767) for x in array_list]

# Perform averaging
array_out = np.nanmean(array_list1, axis=0)
array_out1 = array_out*0.004999999888241291

# Get metadata from one of the input files
with rio.open(smr_list1[0]) as src:
    meta = src.meta

meta.update(dtype=rio.float32)

# Write output file
with rio.open('OB.DAAC Data/Sea Surface Temperature/Mapped/Rasters/Mean_SST/{var}_smr_allyrs_mean.tif'.format(var=var), 'w', **meta) as dst:
    dst.write(array_out1.astype(rio.float32), 1)


#%% SST Winter
dates = ['20220101', '20220201', '20220301', '20220401', '20221001', '20221101', '20221201']
win_list = [select_months(i, 143, 151, dates) for i in file_list]
win_list1 = list(filter(None.__ne__, win_list))

# Read all data as a list of numpy arrays 
array_list = [read_file(x) for x in win_list1]
array_list1 = [make_null(x, -32767) for x in array_list]

# Perform averaging
array_out = np.nanmean(array_list1, axis=0)
array_out1 = array_out*0.004999999888241291

# Get metadata from one of the input files
with rio.open(win_list1[0]) as src:
    meta = src.meta

meta.update(dtype=rio.float32)

# Write output file
with rio.open('OB.DAAC Data/Sea Surface Temperature/Mapped/Rasters/Mean_SST/{var}_win_{yr}_mean.tif'.format(var=var, yr=year), 'w', **meta) as dst:
    dst.write(array_out1.astype(rio.float32), 1)
    

#%%SST Winter for entire study period
years = range(2002, 2023)
win_md = ['0101', '0201', '0301', '0401', '1001', '1101', '1201']

data_dir = 'OB.DAAC Data/{var}/{folder}'.format(var=var, folder=folder) # Or sys.argv[1]
file_list = glob(os.path.join(data_dir, '*.tif'))

dates = [str(y)+x for y in years for x in win_md]
win_list = [select_months(i, 143, 151, dates) for i in file_list]
win_list1 = list(filter(None.__ne__, win_list))

array_list = [read_file(x) for x in win_list1]
array_list1 = [make_null(x, -32767) for x in array_list]

# Perform averaging
array_out = np.nanmean(array_list1, axis=0)
array_out1 = array_out*0.004999999888241291

# Get metadata from one of the input files
with rio.open(smr_list1[0]) as src:
    meta = src.meta

meta.update(dtype=rio.float32)

# Write output file
with rio.open('OB.DAAC Data/Sea Surface Temperature/Mapped/Rasters/Mean_SST/{var}_win_allyrs_mean.tif'.format(var=var), 'w', **meta) as dst:
    dst.write(array_out1.astype(rio.float32), 1)


#%%Do for ChA
var='Chlorophyll-A'
folder='/Mapped/Rasters'
year = '2022'

data_dir = 'OB.DAAC Data/{var}/{folder}'.format(var=var, folder=folder) # Or sys.argv[1]
file_list = glob(os.path.join(data_dir, '*.tif'))

cha_dates = ['20220501', '20220601', '20220701', '20220801', '20220901']

smr_cha_list = [select_months(i, 133, 141, cha_dates) for i in file_list]
smr_cha_list1 = list(filter(None.__ne__, smr_cha_list))

array_list = [read_file(x) for x in smr_cha_list1]

array_list1 = [make_null(x, -32767) for x in array_list]

# Perform averaging, multiply by correction factor
array_out = np.nanmean(array_list1, axis=0)
array_out1 = array_out*1

# Get metadata from one of the input files
with rio.open(smr_cha_list1[0]) as src:
    meta = src.meta

meta.update(dtype=rio.float32)

# Write output file
with rio.open('OB.DAAC Data/Chlorophyll-A/Mapped/Rasters/Mean_ChA/cha_smr_{yr}_mean.tif'.format(yr=year), 'w', **meta) as dst:
    dst.write(array_out1.astype(rio.float32), 1)
   
#array_mst = (array_out1-7.168162)/7.064486
#with rasterio.open('C:/Users/Lauren.Dares/OneDrive - Ocean Wise Conservation Association/Documents/OB.DAAC Data/Chlorophyll-A/Mapped/Rasters/Mean ChA/cha_smr_2019_mst.tif', 'w', **meta) as dst:
#    dst.write(array_mst.astype(rasterio.float32), 1)

#%%ChA Summer for entire study period
years = range(2002, 2023)
smr_md = ['0501', '0601', '0701', '0801', '0901']

data_dir = 'OB.DAAC Data/{var}/{folder}'.format(var=var, folder=folder) # Or sys.argv[1]
file_list = glob(os.path.join(data_dir, '*.tif'))

dates = [str(y)+x for y in years for x in smr_md]
smr_list = [select_months(i, 133, 141, dates) for i in file_list]
smr_list1 = list(filter(None.__ne__, smr_list))

array_list = [read_file(x) for x in smr_list1]
array_list1 = [make_null(x, -32767) for x in array_list]

# Perform averaging
array_out = np.nanmean(array_list1, axis=0)
array_out1 = array_out*1

# Get metadata from one of the input files
with rio.open(smr_list1[0]) as src:
    meta = src.meta

meta.update(dtype=rio.float32)

# Write output file
with rio.open('OB.DAAC Data/{var}/Mapped/Rasters/Mean_ChA/{var}_smr_allyrs_mean.tif'.format(var=var), 'w', **meta) as dst:
    dst.write(array_out1.astype(rio.float32), 1)

#%% ChA Winter
cha_dates = ['20211001', '20211101', '20211201', '20220101', '20220201','20220301', '20220401']

win_cha_list = [select_months(i, 133, 141, cha_dates) for i in file_list]
win_cha_list1 = list(filter(None.__ne__, win_cha_list))

array_list = [read_file(x) for x in win_cha_list1]
array_list1 = [make_null(x, -32767) for x in array_list]

# Perform averaging
array_out = np.nanmean(array_list1, axis=0)
array_out1 = array_out*1

# Get metadata from one of the input files
with rasterio.open(win_cha_list1[0]) as src:
    meta = src.meta

meta.update(dtype=rasterio.float32)

# Write output file
with rasterio.open('OB.DAAC Data/Chlorophyll-A/Mapped/Rasters/Mean_ChA/cha_win_{yr}_mean.tif'.format(yr=year), 'w', **meta) as dst:
    dst.write(array_out1.astype(rasterio.float32), 1)


#array_mst = (array_out1-6.079538)/7.128631
#with rasterio.open('C:/Users/Lauren.Dares/OneDrive - Ocean Wise Conservation Association/Documents/OB.DAAC Data/Chlorophyll-A/Mapped/Rasters/Mean ChA/cha_win_2019_mst1.tif', 'w', **meta) as dst:
#    dst.write(array_mst.astype(rasterio.float32), 1)

#%%ChA Winter for entire study period
years = range(2002, 2023)
win_md = ['0101', '0201', '0301', '0401', '1001', '1101', '1201']

data_dir = 'OB.DAAC Data/{var}/{folder}'.format(var=var, folder=folder) # Or sys.argv[1]
file_list = glob(os.path.join(data_dir, '*.tif'))

dates = [str(y)+x for y in years for x in win_md]
win_list = [select_months(i, 133, 141, dates) for i in file_list]
win_list1 = list(filter(None.__ne__, win_list))

array_list = [read_file(x) for x in win_list1]
array_list1 = [make_null(x, -32767) for x in array_list]

# Perform averaging
array_out = np.nanmean(array_list1, axis=0)
array_out1 = array_out*1

# Get metadata from one of the input files
with rio.open(win_list1[0]) as src:
    meta = src.meta

meta.update(dtype=rio.float32)

# Write output file
with rio.open('OB.DAAC Data/{var}/Mapped/Rasters/Mean_ChA/{var}_win_allyrs_mean.tif'.format(var=var), 'w', **meta) as dst:
    dst.write(array_out1.astype(rio.float32), 1)


#%%Resample depth rasters to 4 km x 4 km
depth = rx.open_rasterio('env_means/winter/gebco_2021_resampled_4km.tif').squeeze()
depth = np.where(depth>=0, np.nan, depth)
# Get metadata from one of the input files
with rio.open('env_means/winter/gebco_2021_resampled_4km.tif') as src:
    meta = src.meta

meta.update(dtype=rasterio.float32)

# Write output file
with rasterio.open('gebco_2021_lessthanzero.tif', 'w', **meta) as dst:
    dst.write(depth.astype(rasterio.float32), 1)

ref = gdal.Open("env_means/summer/Chlorophyll-A_smr_allyrs_mean.tif", 0)
ref_trans = ref.GetGeoTransform()
x_res = ref_trans[1]
y_res = -ref_trans[5]

kwargs = {"format": "GTiff", "xRes":x_res, "yRes":y_res}

depth1 = gdal.Warp('env_means/winter/gebco_2021_resampled_4km_python1.tif', 'GEBCO 2021 Bathymetry/GEBCO_15_Nov_2021_6fcf1f1d3513/gebco_2021_n57.0_s45.0_w-139.0_e-120.0.tif', **kwargs)
