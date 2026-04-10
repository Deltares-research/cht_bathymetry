"""
Cloud-Optimized GeoTIFF (COG) bathymetry dataset.

Provides :class:`BathymetryDatasetCOG` for reading depth data from local or
S3-hosted COG files, selecting an appropriate internal overview level based on
the requested resolution.
"""

import os
from pathlib import Path

import numpy as np
import rasterio
import rioxarray
import xarray as xr

from .dataset import BathymetryDataset


class BathymetryDatasetCOG(BathymetryDataset):
    """
    Bathymetry dataset backed by a Cloud-Optimized GeoTIFF file.

    The file may reside locally or be fetched on first use from an S3 bucket
    when ``s3_bucket`` and ``s3_key`` metadata attributes are present.
    """

    def __init__(self, name: str, path: str) -> None:
        """
        Initialise a COG dataset.

        Parameters
        ----------
        name : str
            Short identifier for the dataset (also used as the filename stem
            and for metadata lookup).
        path : str
            Directory that contains the COG file and its metadata TOML.
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
        Read depth values from the COG for a bounding box.

        Selects the most detailed internal overview whose pixel size is still
        no larger than *max_cell_size* (in metres).  Downloads the file from
        S3 first if it is not yet available locally.

        Parameters
        ----------
        xl : list[float]
            ``[x_min, x_max]`` bounds in the dataset's coordinate system.
        yl : list[float]
            ``[y_min, y_max]`` bounds in the dataset's coordinate system.
        max_cell_size : float, optional
            Maximum pixel size in metres used to choose the overview level.
            Default is ``1000.0``.
        waitbox : None, optional
            Reserved for a progress-dialog object; currently unused.

        Returns
        -------
        x : np.ndarray
            1-D array of x coordinates of the returned grid.
        y : np.ndarray
            1-D array of y coordinates of the returned grid.
        z : np.ndarray
            2-D depth array (NaN where no data).  Returns scalar ``np.nan``
            for all three values when the bounding box lies outside the
            dataset extent.
        """

        if not self.path.exists():
            if hasattr(self, "s3_key") and hasattr(self, "s3_bucket"):
                # Download first !
                self.download()

        # First find appropriate overview level based on max pixel size
        with rasterio.open(self.path) as src:
            overview_level, ok = get_appropriate_overview_level(src, max_cell_size)

        if ok:
            rds = rioxarray.open_rasterio(
                self.path, masked=False, overview_level=overview_level
            )
        else:
            rds = rioxarray.open_rasterio(self.path, masked=False)

        # Check if bounding box covers the bounds of the dataset
        if (
            xl[1] < rds.rio.bounds()[0]
            or xl[0] > rds.rio.bounds()[2]
            or yl[1] < rds.rio.bounds()[1]
            or yl[0] > rds.rio.bounds()[3]
        ):
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

    def download(self) -> None:
        """
        Download the COG file from the configured S3 bucket.

        Uses ``self.database.s3_client`` together with ``self.s3_bucket``,
        ``self.s3_key``, and ``self.filename`` to construct the S3 object key
        and target local path.
        """
        print(f"Downloading {self.filename} from S3")
        print("This may take a while...")

        key = f"{self.s3_key}/{self.filename}"
        filename = os.path.join(self.local_path, self.filename)
        try:
            self.database.s3_client.download_file(
                Bucket=self.s3_bucket,
                Key=key,
                Filename=filename,
            )

            print("Downloading done.")

        except Exception:
            print(f"Failed to download {key}. Skipping dataset.")


def get_appropriate_overview_level(
    src: rasterio.io.DatasetReader, max_pixel_size: float
) -> tuple[int | None, bool]:
    """
    Choose the best internal overview level for a target pixel size.

    Iterates over the overview pyramid of *src* and returns the index of the
    highest-resolution overview whose pixel size is still no larger than
    *max_pixel_size*.

    Parameters
    ----------
    src : rasterio.io.DatasetReader
        An open rasterio dataset.
    max_pixel_size : float
        Maximum acceptable pixel size in metres.

    Returns
    -------
    selected_overview : int or None
        Zero-based overview index, or ``None`` when no suitable overview
        exists.
    ok : bool
        ``True`` when a suitable overview was found, ``False`` otherwise.
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
        return None, False

    # Calculate the resolution for each overview by multiplying the original resolution by the overview factor
    resolutions = [
        (original_resolution[0] * factor, original_resolution[1] * factor)
        for factor in overview_levels
    ]

    if resolutions[0][0] > max_pixel_size and resolutions[0][1] > max_pixel_size:
        # Even the lowest overview is larger than max_pixel_size
        return None, False

    # Find the highest overview level that is smaller than or equal to the max_pixel_size
    selected_overview = 0
    for i, (x_res, y_res) in enumerate(resolutions):
        if x_res <= max_pixel_size and y_res <= max_pixel_size:
            selected_overview = i
        else:
            break

    return selected_overview, True
