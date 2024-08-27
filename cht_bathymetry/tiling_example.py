# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 13:48:55 2022

@author: ormondt
"""
import glob
import cht

from cht.bathymetry.tiling import TiledDataset

name      = "new_jersey_delaware_coned_2015"
long_name = "2015 USGS CoNED Topobathymetric Model: New Jersey and Delaware (1888 - 2014)"
database_path = "c:\\work\\delftdashboard\\data\\bathymetry"

url = "https://chs.coast.noaa.gov/htdata/raster2/elevation/NewJersey_Delaware_Coned_Topobathy_DEM_2015_5040/"
skip_file = "skiplist.txt"
skip_file = None

attributes = {}

ds = TiledDataset(name, database_path, attrs=attributes)

file_list = glob.glob("*.tif")

ds.build(file_list,
         make_base_tiles=True,
         make_high_level_tiles=True,
         url=url,
         crs=crs,
         skip_file=skip_file)
