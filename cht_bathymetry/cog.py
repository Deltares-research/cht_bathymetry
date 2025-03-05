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


class BathymetryDatasetCOG(BathymetryDataset):
    """
    Cloud-optimized GeoTiFF dataset class
    """

    def __init__(self, name: str, path: str):
        """
        Initialize the BathymetryDatasetCOG class.

        Parameters:
        name (str): The name of the dataset.
        path (str): The path to the dataset.
        """
        super().__init__()

        self.name: str = name
        self.path: str = path
        self.local_path: str = path
        self.read_metadata()
        self.data: xr.Dataset = xr.Dataset()
        self.path: Path = Path(self.local_path) / self.filename

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

        if not self.path.exists():
            if hasattr(self, "s3_key") and hasattr(self, "s3_bucket"):
                # Download first !
                self.download()

        # First find appropriate overview level based on max pixel size
        with rasterio.open(self.path) as src:
            overview_level = get_appropriate_overview_level(src, max_cell_size)

        rds = rioxarray.open_rasterio(
            self.path, masked=False, overview_level=overview_level
        )

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

    def download(self) -> None:
        """
        Download the COG file from S3.
        """
        print(f"Downloading {self.filename} from S3")
        print("This may take a while...")

        key = f"{self.s3_key}/{self.filename}"
        filename = os.path.join(self.local_path, self.filename)
        try:
            self.database.s3_client.download_file(
                Bucket=self.s3_bucket,  # assign bucket name
                Key=key,  # key is the file name
                Filename=filename,
            )

            print("Downloading done.")

        except Exception as e:
            print(f"Failed to download {key}. Skipping dataset.")


def get_appropriate_overview_level(
    src: rasterio.io.DatasetReader, max_pixel_size: float
) -> int:
    """
    Given a rasterio dataset `src` and a desired `max_pixel_size`,
    determine the appropriate overview level (zoom level) that fits
    the maximum resolution allowed by `max_pixel_size`.

    Parameters:
    src (rasterio.io.DatasetReader): The rasterio dataset reader object.
    max_pixel_size (float): The maximum pixel size for the resolution.

    Returns:
    int: The appropriate overview level.
    """
    # Get the original resolution (pixel size) in terms of x and y
    original_resolution = src.res  # Tuple of (x_resolution, y_resolution)
    if src.crs.is_geographic:
        original_resolution = (
            original_resolution[0] * 111000,
            original_resolution[1] * 111000,
        )  # Convert to meters
    # Get the overviews for the dataset
    overview_levels = src.overviews(
        1
    )  # Overview levels for the first band (if multi-band, you can adjust this)

    # If there are no overviews, return 0 (native resolution)
    if not overview_levels:
        return 0

    # Calculate the resolution for each overview by multiplying the original resolution by the overview factor
    resolutions = [
        (original_resolution[0] * factor, original_resolution[1] * factor)
        for factor in overview_levels
    ]

    # Find the highest overview level that is smaller than or equal to the max_pixel_size
    selected_overview = 0
    for i, (x_res, y_res) in enumerate(resolutions):
        if x_res <= max_pixel_size and y_res <= max_pixel_size:
            selected_overview = i
        else:
            break

    return selected_overview
