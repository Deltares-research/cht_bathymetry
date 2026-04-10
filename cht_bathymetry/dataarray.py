"""
xarray DataArray bathymetry dataset.

Provides :class:`BathymetryDatasetDataArray` for wrapping an in-memory
:class:`xarray.DataArray` (e.g. loaded from a GeoTIFF) as a bathymetry
dataset that can be queried through the standard database interface.
"""

import numpy as np
import xarray as xr

from .dataset import BathymetryDataset


class BathymetryDatasetDataArray(BathymetryDataset):
    """
    Bathymetry dataset backed by an :class:`xarray.DataArray`.

    Intended for datasets that are already loaded into memory (e.g. opened
    with :func:`rioxarray.open_rasterio`) rather than read tile-by-tile from
    disk.
    """

    def __init__(self, da: xr.DataArray, name: str) -> None:
        """
        Initialise a DataArray dataset.

        Parameters
        ----------
        da : xr.DataArray
            Source DataArray with spatial coordinates and a ``rio`` accessor.
        name : str
            Short identifier for the dataset.
        """
        super().__init__()

        self.name: str = name
        self.data = da

    def get_data(
        self,
        xl: list[float],
        yl: list[float],
        max_cell_size: float = 1000.0,
        waitbox: None = None,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Return depth values clipped to a bounding box.

        Parameters
        ----------
        xl : list[float]
            ``[x_min, x_max]`` bounds in the dataset's coordinate system.
        yl : list[float]
            ``[y_min, y_max]`` bounds in the dataset's coordinate system.
        max_cell_size : float, optional
            Maximum cell size in metres (unused here — resolution is fixed by
            the DataArray).  Default is ``1000.0``.
        waitbox : None, optional
            Reserved for a progress-dialog object; currently unused.

        Returns
        -------
        x : np.ndarray
            1-D array of x coordinates.
        y : np.ndarray
            1-D array of y coordinates.
        z : np.ndarray
            2-D depth array (NaN where no data).  Returns scalar ``np.nan``
            for all three values when the bounding box lies outside the
            dataset extent.
        """

        rds = self.data

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
        z = data.values[:, :]
        # Squeeze z to remove singleton dimensions
        z = np.squeeze(z)

        z[z == rds.rio.nodata] = np.nan

        return x, y, z

    def get_bbox(self) -> list[float]:
        """
        Return the bounding box of the dataset.

        Returns
        -------
        list[float]
            ``[minx, miny, maxx, maxy]`` in the dataset's native CRS.
        """
        return [
            self.data.rio.bounds()[0],
            self.data.rio.bounds()[1],
            self.data.rio.bounds()[2],
            self.data.rio.bounds()[3],
        ]

    def get_lon_lat_range(self, **kwargs) -> tuple[list[float] | None, list[float] | None]:
        """
        Return the longitude and latitude range of the dataset.

        Reprojects the bounding box to geographic coordinates when the
        dataset's CRS is projected.

        Returns
        -------
        lon_range : list[float] or None
            ``[lon_min, lon_max]``, or ``None`` when the bounding box cannot
            be determined.
        lat_range : list[float] or None
            ``[lat_min, lat_max]``, or ``None`` when the bounding box cannot
            be determined.
        """
        # Get bbox in WGS 84
        # Convert bounds to lon/lat
        bbox = self.get_bbox(**kwargs)
        if bbox is None:
            return None, None
        else:
            if self.crs.is_geographic:
                lon_range = [bbox[0], bbox[2]]
                lat_range = [bbox[1], bbox[3]]
            else:
                # Convert to lon/lat range
                lon_range, lat_range = self.crs.transform_bounds(bbox[0:2], bbox[2:4])
        return lon_range, lat_range
