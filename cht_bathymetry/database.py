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
import boto3
from botocore import UNSIGNED
from botocore.client import Config

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
    
    def __init__(self,
                 path=None,
                 s3_bucket=None,
                 s3_key=None,
                 s3_region=None):
        self.path    = None
        self.dataset = []
        self.initialized = False
        self.s3_client = None
        self.s3_bucket = s3_bucket
        self.s3_key = s3_key
        self.s3_region = s3_region
        if path is not None:
            self.initialize(path)
    
    def initialize(self, path):
        self.path    = path
        self.dataset = []
        self.read()
        self.initialized = True
       
    def read(self):
        """
        Reads meta-data of all datasets in the database. 
        """

        if self.path is None:
            print("Path to bathymetry database not set !")
            return
        
        # Check if the path exists. If not, create it.
        if not os.path.exists(self.path):
            os.makedirs(self.path)

        # Read in database
        tml_file = os.path.join(self.path, "bathymetry.tml")

        if not os.path.exists(tml_file):
            print("Warning! Bathymetry database file not found: " + tml_file)
            return

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
            elif os.path.exists(os.path.join(path, "metadata.tml")):
                metadata = toml.load(os.path.join(path, "metadata.tml"))
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

    def write(self):
        """
        """        
        pass
        # # Write in database
        # tml_file = os.path.join(self.path, "bathymetry.tml")
        # datasets = {"dataset": []}
        # for dataset in self.dataset:
        #     datasets["dataset"].append(dataset.get_metadata())
        # dict2yaml(tml_file, datasets, sort_keys=False)

    def check_online_database(self):
        if self.s3_client is None:
            self.s3_client = boto3.client('s3', config=Config(signature_version=UNSIGNED))
        if self.s3_bucket is None:
            return
        # First download a copy of bathymetry.tml and call it bathymetry_s3.tml
        key = f"{self.s3_key}/bathymetry.tml"
        filename = os.path.join(self.path, "bathymetry_s3.tml")
        print("Updating bathymetry database ...")
        try:
            self.s3_client.download_file(Bucket=self.s3_bucket,     # assign bucket name
                                         Key=key,           # key is the file name
                                         Filename=filename) # storage file path
        except:
            # Download failed
            print(f"Failed to download {key} from {self.s3_bucket}. Database will not be updated.")
            return

        # Read bathymetry_s3.tml
        short_name_list, long_name_list, source_name_list = self.dataset_names()
        datasets_s3 = toml.load(filename)
        datasets_added = False
        added_names = []
        # Loop through s3 datasets, and check whether they exist in the local database.
        # If so, check if the metadata also exists. If not, make local folder and download the metadata.
        # Additionally, check if available_tiles.nc in s3 and not in local database, download it.
        for d in datasets_s3["dataset"]:
            # Get list of existing datasets
            s3_name = d["name"]
            if s3_name not in short_name_list:
                # Dataset not in local database
                print(f"Adding bathymetry dataset {s3_name} to local database ...")
                # Create folder and download metadata
                path = os.path.join(self.path, s3_name)
                os.makedirs(path, exist_ok=True)
                key = f"{self.s3_key}/{s3_name}/metadata.tml"
                filename = os.path.join(path, "metadata.tml")
                # Download metadata
                try:
                    self.s3_client.download_file(Bucket=self.s3_bucket, # assign bucket name
                                                Key=key,               # key is the file name
                                                Filename=filename)     # storage file path
                except Exception as e:
                    print(e)
                    print(f"Failed to download {key}. Skipping dataset.")
                    continue
                # Read metadata.tml
                metadata = toml.load(filename)
                if "available_tiles" in metadata:
                    if metadata["available_tiles"]:
                        # Download available_tiles.nc
                        key = f"{self.s3_key}/{s3_name}/available_tiles.nc"
                        filename = os.path.join(path, "available_tiles.nc")
                        try:
                            self.s3_client.download_file(Bucket=self.s3_bucket, # assign bucket name
                                                        Key=key,               # key is the file name
                                                        Filename=filename)
                        except Exception as e:
                            print(f"Failed to download {key}. Skipping dataset.")
                            continue    
                key = f"{self.s3_key}/{s3_name}/index.html"
                filename = os.path.join(path, "index.html")
                try:
                    self.s3_client.download_file(Bucket=self.s3_bucket, # assign bucket name
                                                Key=key,               # key is the file name
                                                Filename=filename)
                except Exception as e:
                    print(f"index.html could not be downloaded")
                # Necessary data has been downloaded    
                datasets_added = True
                added_names.append(s3_name)
        # Write new local bathymetry.tml
        if datasets_added:
            d = {}
            d["dataset"] = []
            for name in short_name_list:
                d["dataset"].append({"name": name})
            for name in added_names:
                d["dataset"].append({"name": name})
            # Now write the new bathymetry.tml
            with open(os.path.join(self.path, "bathymetry.tml"), "w") as tml:
                toml.dump(d, tml)            
            # Read the database again
            self.dataset = []
            self.read()
        else:
            print("No new datasets were added to the local database.")

    def get_bathymetry_on_points(self, xz, yz, dxmin, crs, bathymetry_list, method="linear"):
        zz = self.get_bathymetry_on_grid(xz, yz, crs, bathymetry_list, method=method, coords="points", dxmin=dxmin)
        return zz

    def get_bathymetry_on_grid(self, xz, yz, crs, bathymetry_list,
                               method="linear",
                               coords="grid",
                               dxmin=1.0e6,
                               waitbox=None):

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
                                              max_cell_size=dx,
                                              waitbox=waitbox)

                # If zb equal np.nan, then there is no data
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
