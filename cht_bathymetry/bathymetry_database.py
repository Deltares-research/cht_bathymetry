# -*- coding: utf-8 -*-
"""
Created on Sun Apr 25 10:58:08 2021

@author: Maarten van Ormondt
"""

import os
import numpy as np
import yaml
import toml
from pyproj import Transformer
from cht_utils.misc_tools import interp2

from .netcdf_tiles_v1 import BathymetryDatasetNetCDFTilesV1
from .netcdf_tiles_v2 import BathymetryDatasetNetCDFTilesV2
from .tiled_web_map import BathymetryDatasetTiledWebMap

class BathymetryDatabase:
    """
    The main Bathymetry Database class
    
    :param pth: Path name where bathymetry tiles will be cached.
    :type pth: string            
    """
    
    def __init__(self, pth):
        if pth:
            self.path    = pth
            self.dataset = []
            self.read()
            self.initialized = True
        else:
            self.initialized = False
    
    def initialize(self, pth):
        if not self.initialized:
            self.path    = pth
            self.dataset = []
            self.read()
        self.initialized = True
       
    def read(self):
        """
        Reads meta-data of all datasets in the database. 
        """
        
        # Read in database
        # yml_file = os.path.join(self.path, "bathymetry.yml")
        # datasets = yaml2dict(yml_file)
        tml_file = os.path.join(self.path, "bathymetry.tml")
        datasets = toml.load(tml_file)

        for d in datasets["dataset"]:

            name = d["name"]

            if "path" in d:
                path = d["path"]
            else:
                path = os.path.join(self.path, name)

            # Read the meta data for this dataset
            fname = os.path.join(path, name + ".tml")

            if os.path.exists(fname):
                metadata = toml.load(fname)
                dataset_format = metadata["format"]
            else:
                print("Could not find metadata file for dataset " + name + " ! Skipping dataset.")
                continue

            if dataset_format == "netcdf_tiles_v1":
                dataset = BathymetryDatasetNetCDFTilesV1(name, path)
            elif dataset_format == "netcdf_tiles_v2":
                dataset = BathymetryDatasetNetCDFTilesV2(name, path)
            elif dataset_format == "tiled_web_map":
                dataset = BathymetryDatasetTiledWebMap(name, path)

            dataset.database = self    
            
            self.dataset.append(dataset)


    # def get_crs(self, dataset_name):
    #     """
    #     Returns coordinate reference system (CRS) of dataset

    #     :param dataset_name: Name of requested bathymetry dataset.
    #     :type dataset_name: str
    #     :return: CRS  
    #     :rtype: str
    #     """
        
    #     # Read in data from database
    #     # Find corresponding dataset dataset
    #     for d in self.dataset:
    #         if d.name == dataset_name:
    #             return d.crs

    def get_bathymetry_on_points(self, xz, yz, dxmin, crs, bathymetry_list, method="linear"):
        zz = self.get_bathymetry_on_grid(xz, yz, crs, bathymetry_list, method=method, coords="points", dxmin=dxmin)
        return zz

    def get_bathymetry_on_grid(self, xz, yz, crs, bathymetry_list, method="linear", coords="grid", dxmin=1.0e6):

        if xz.ndim == 2:
            # xy and yz are a grid
            zz = np.full(xz.shape, np.nan)
            dx = np.sqrt((xz[0,1] - xz[0,0])**2 + (yz[0,1] - yz[0,0])**2)
            dy = np.sqrt((xz[1,0] - xz[0,0])**2 + (yz[1,0] - yz[0,0])**2)
        else:
            if coords == "grid":
                zz = np.full((len(yz), len(xz)), np.nan)
                dx = xz[1] - xz[0]
                dy = yz[1] - yz[0]
                xz, yz = np.meshgrid(xz, yz)
            else:    
                zz = np.full(xz.shape, np.nan)

        # Determine resolution to get bathy data
        if coords == "grid":
            # Resolution follow from grid
            if crs.is_geographic:
                dx = min(111111.0 * dx,
                        111111.0 * dy * np.cos(np.pi * np.max(np.abs(yz)) / 180.0))
            else:
                dx = min(dx, dy)
        else:
            dx = dxmin        

        # Loop through bathymetry datasets
        for ibathy, bathymetry in enumerate(bathymetry_list):

            dataset_name = bathymetry["name"] # name
            zmin         = bathymetry["zmin"]
            zmax         = bathymetry["zmax"]

            dataset = self.get_dataset(dataset_name)

            transformer = Transformer.from_crs(crs,
                                               dataset.crs,
                                               always_xy=True)

            if np.isnan(zz).any():
                xzb, yzb = transformer.transform(xz, yz)
                xmin = np.nanmin(np.nanmin(xzb))
                xmax = np.nanmax(np.nanmax(xzb))
                ymin = np.nanmin(np.nanmin(yzb))
                ymax = np.nanmax(np.nanmax(yzb))
                ddx = 0.05 * (xmax - xmin)
                ddy = 0.05 * (ymax - ymin)
                xl = [xmin - ddx, xmax + ddx]
                yl = [ymin - ddy, ymax + ddy]

                # Get DEM data
                xb, yb, zb = dataset.get_data(xl,
                                              yl,
                                              max_cell_size=dx)

                # If zb equal np.nan, then there is not data
                if not np.isnan(zb).all():
                    zb[np.where(zb < zmin)] = np.nan
                    zb[np.where(zb > zmax)] = np.nan
                    zz1 = interp2(xb, yb, zb, xzb, yzb, method=method)
                    isn = np.where(np.isnan(zz))
                    zz[isn] = zz1[isn]

        return zz

    def get_dataset(self, name):
        for dataset in self.dataset:
            if dataset.name == name:
                return dataset
        return None

    def dataset_names(self, source=None):
        short_name_list = []
        long_name_list = []
        source_name_list = []
        for dataset in self.dataset:
            ok = False
            if source:
                if dataset.source == source:
                    ok = True
            else:
                ok = True
            if ok:
                short_name_list.append(dataset.name)
                long_name_list.append(dataset.long_name)
                source_name_list.append(dataset.source)
        return short_name_list, long_name_list, source_name_list

    def sources(self):

        sources = []
        source_names = []

        for dataset in self.dataset:
            source = dataset.source
            if source in source_names:
                # Existing source
                for src in sources:
                    if src.name == source:
                        src.dataset.append(dataset)
            else:
                # New source
                src = BathymetrySource(source)
                src.dataset.append(dataset)
                sources.append(src)
                source_names.append(source)

        return source_names, sources

class BathymetrySource:    
    def __init__(self, name):        
        self.name    = name
        self.dataset = []

def dict2yaml(file_name, dct, sort_keys=False):
    yaml_string = yaml.dump(dct, sort_keys=sort_keys)    
    file = open(file_name, "w")  
    file.write(yaml_string)
    file.close()

def yaml2dict(file_name):
    file = open(file_name,"r")
    dct = yaml.load(file, Loader=yaml.FullLoader)
    return dct

bathymetry_database = BathymetryDatabase(None)
