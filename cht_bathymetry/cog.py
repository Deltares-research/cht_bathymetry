# -*- coding: utf-8 -*-
"""
Created on Sun Apr 25 10:58:08 2021

@author: Maarten van Ormondt
"""
import os
import xarray as xr
from pathlib import Path
import rasterio
import numpy as np
# from rasterio.enums import Resampling
# from rasterio.windows import from_bounds
import rioxarray

# from hydromt import DataCatalog

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
        # if self.crs is None:
        #     if self.path.exists():
        #         with rasterio.open(self.path) as dataset:
        #             self.crs = dataset.crs
        #     else:
        #         print("Warning: CRS not defined for dataset {self.name} as the file does not exist (yet).")   
            
    def get_data(self, xl, yl, max_cell_size=1000.0, waitbox=None):
        """
        Reads data from database. Returns xarray dataset in same coordinate system (3857) as dataset. Resolution is determined by max_cell_size.
        """

        if not self.path.exists():
            if hasattr(self, "s3_key") and hasattr(self, "s3_bucket"):
                # Download first !
                self.download()

        # First find appropriate overview level based on max pixel size
        with rasterio.open(self.path) as src:
            overview_level = get_appropriate_overview_level(src, max_cell_size)

        rds = rioxarray.open_rasterio(self.path,
                                      masked=False,
                                      overview_level=overview_level)

        data = rds.rio.clip_box(
            minx=xl[0],
            miny=yl[0],
            maxx=xl[1],
            maxy=yl[1],
        )
        x = data.x.values[:]
        y = data.y.values[:]
        z = data.values[0,:,:]
        z[z == rds.rio.nodata] = np.nan

        rds.close()

        return x, y, z

    def download(self):

        # Download the COG file
        print(f"Downloading {self.filename} from S3")
        print("This may take a while...")

        key = f"{self.s3_key}/{self.filename}"
        filename = os.path.join(self.local_path, self.filename)
        try:
            self.database.s3_client.download_file(Bucket=self.s3_bucket, # assign bucket name
                                                  Key=key,               # key is the file name
                                                  Filename=filename)

            print("Downloading done.")

        except Exception as e:
            print(f"Failed to download {key}. Skipping dataset.")


def get_appropriate_overview_level(src, max_pixel_size):
    """
    Given a rasterio dataset `src` and a desired `max_pixel_size`, 
    determine the appropriate overview level (zoom level) that fits 
    the maximum resolution allowed by `max_pixel_size`.
    """
    # Get the original resolution (pixel size) in terms of x and y
    original_resolution = src.res  # Tuple of (x_resolution, y_resolution)
    if src.crs.is_geographic:
        original_resolution = original_resolution[0] * 111000, original_resolution[1] * 111000  # Convert to meters
    # Get the overviews for the dataset
    overview_levels = src.overviews(1)  # Overview levels for the first band (if multi-band, you can adjust this)
    
    # If there are no overviews, return 0 (native resolution)
    if not overview_levels:
        return 0
    
    # Calculate the resolution for each overview by multiplying the original resolution by the overview factor
    resolutions = [(original_resolution[0] * factor, original_resolution[1] * factor) for factor in overview_levels]
    
    # Find the highest overview level that is smaller than or equal to the max_pixel_size
    selected_overview = 0
    for i, (x_res, y_res) in enumerate(resolutions):
        if x_res <= max_pixel_size and y_res <= max_pixel_size:
            selected_overview = i
        else:
            break

    return selected_overview
