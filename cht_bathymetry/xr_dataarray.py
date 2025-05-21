# -*- coding: utf-8 -*-
"""
This module defines the BathymetryDatasetCOG class, which represents a cloud-optimized GeoTIFF (COG) dataset for bathymetry data. 
It provides methods to initialize the dataset, read data from the dataset, and download the dataset from an S3 bucket.

Classes:
    BathymetryDatasetCOG: A class for handling cloud-optimized GeoTIFF bathymetry datasets.

Functions:
    get_appropriate_overview_level(src: rasterio.io.DatasetReader, max_pixel_size: float) -> int:
        Determines the appropriate overview level for a rasterio dataset based on the maximum pixel size.

Usage:
    from .cog import BathymetryDatasetCOG
"""

import os
import xarray as xr
from pathlib import Path
import rasterio
import numpy as np
import rioxarray

from .dataset import BathymetryDataset


class BathymetryDatasetDataArray(BathymetryDataset):
    """
    XR DataArray dataset class
    """

    def __init__(self, name: str, path: str):
        """
        Initialize the BathymetryDatasetCOG class.

        Parameters:
        name (str): The name of the dataset.
        path (str): The path to the dataset (dummy value).
        """
        super().__init__()

        self.name: str = name
        self.data: xr.DataArray = xr.DataArray()

    def get_data(
        self,
        xl: list[float],
        yl: list[float],
        max_cell_size: float = 1000.0,
        waitbox: None = None,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Reads data from the database. Returns arrays x, y, z in the same coordinate system as the dataset. 
        Resolution is determined by max_cell_size.

        Parameters:
        xl (list[float]): List of x coordinates (longitude).
        yl (list[float]): List of y coordinates (latitude).
        max_cell_size (float): Maximum cell size for the resolution. Default is 1000.0.
        waitbox (None): Placeholder for a waitbox object. Default is None.

        Returns:
        tuple[np.ndarray, np.ndarray, np.ndarray]: Returns three numpy arrays representing x, y, and z coordinates.
        """

        rds = self.data

        # Check if bounding box covers the bounds of the dataset
        if xl[1] < rds.rio.bounds()[0] or xl[0] > rds.rio.bounds()[2] or yl[1] < rds.rio.bounds()[1] or yl[0] > rds.rio.bounds()[3]:
            # print("Bounding box is outside the dataset bounds.")
            return np.nan, np.nan, np.nan

        data = rds.rio.clip_box(
            minx=xl[0],
            miny=yl[0],
            maxx=xl[1],
            maxy=yl[1],
        )
        x = data.x.values[:]
        y = data.y.values[:]
        z = data.values[0, :, :]
        z[z == rds.rio.nodata] = np.nan

        rds.close()

        return x, y, z

