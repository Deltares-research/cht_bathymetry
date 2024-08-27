# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 13:48:55 2022

@author: ormondt
"""
import glob

from cht.bathymetry.tiling import TiledDataset

name      = "gebco22"

# Path where tiles will be stored
database_path = "c:\\work\\delftdashboard\\data\\bathymetry"

# Download from https://www.gebco.net/data_and_products/gridded_bathymetry_data/
file_list = glob.glob("c:\\work\\projects\\delftdashboard\\gebco22\\*.tif")

# These will be written to yml file
attributes = {}
attributes["long_name"]                = "GEBCO 2022"
attributes["source"]                   = "BODC"
attributes["url"]                      = "https://www.gebco.net/data_and_products/gridded_bathymetry_data/"
attributes["vertical_reference_level"] = "MSL"
attributes["vertical_units"]           = "m"
attributes["difference_with_msl"]      = 0.0
attributes["created_by"]               = "Maarten van Ormondt"
attributes["disclaimer"]               = "These data are made available in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE."

ds = TiledDataset()

ds.build(name,
         file_list,
         database_path,
         make_base_tiles=True,
         make_high_level_tiles=True,
         attrs=attributes)
