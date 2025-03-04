# -*- coding: utf-8 -*-
"""
Created on Sun Apr 25 10:58:08 2021

@author: Maarten van Ormondt
"""
import xarray as xr
from pathlib import Path
import rasterio
import time

from hydromt import DataCatalog

from .dataset import BathymetryDataset

class BathymetryDatasetCOG(BathymetryDataset):
    """
    Bathymetry dataset class 

    """

    def __init__(self, name, path):
        super().__init__()
        
        self.name              = name
        self.path              = path
        self.local_path        = path
        self.read_metadata()
        self.data              = xr.Dataset()
        self.path              = Path(self.local_path) / self.filename
        if self.crs is None:
            with rasterio.open(self.path) as dataset:
                self.crs = dataset.crs
            
    def get_data(self, xl, yl, max_cell_size, waitbox=None):
        """
        Reads data from database. Returns xarray dataset in same coordinate system (3857) as dataset. Resolution is determined by max_cell_size.
        """

        bbox=(xl[0], xl[1], yl[0], yl[1])

        time0 = time.time()

        dc = DataCatalog()
        da_elv = dc.get_rasterdataset(
            self.path,
            bbox=bbox,
            buffer=10,
            variables=["elevtn"],
            zoom=(max_cell_size, "meter"),
        )
        time1 = time.time()

        print(f"Time to get data array: {(time1-time0):.2f} s")

        time0 = time.time()
        x = da_elv.x.values[:]
        y = da_elv.y.values[:]
        z = da_elv.values[:]
        time1 = time.time()

        print(f"Time to read data to numpy arrays: {(time1-time0):.2f} s")

        return x, y, z
    
