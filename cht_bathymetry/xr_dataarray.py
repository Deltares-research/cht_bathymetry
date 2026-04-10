"""
xarray DataArray bathymetry dataset (alternate constructor).

Provides a variant of the DataArray-backed bathymetry dataset that accepts a
``(name, path)`` constructor signature matching the other dataset types,
storing the DataArray separately via ``self.data``.

Note: for the primary DataArray-based dataset (constructed directly from an
:class:`xarray.DataArray` object) see :mod:`cht_bathymetry.dataarray`.
"""

import numpy as np
import xarray as xr

from .dataset import BathymetryDataset


class BathymetryDatasetDataArray(BathymetryDataset):
    """
    Bathymetry dataset backed by an :class:`xarray.DataArray`.

    Constructed from a ``(name, path)`` pair for interface consistency with
    tile-based dataset types.  The actual DataArray must be assigned to
    ``self.data`` after construction.
    """

    def __init__(self, name: str, path: str) -> None:
        """
        Initialise the dataset with a name and placeholder path.

        Parameters
        ----------
        name : str
            Short dataset identifier.
        path : str
            Path string (accepted for interface consistency; not used to load
            data directly in this class).
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
        Return depth values clipped to a bounding box.

        Parameters
        ----------
        xl : list[float]
            ``[x_min, x_max]`` bounds in the dataset's coordinate system.
        yl : list[float]
            ``[y_min, y_max]`` bounds in the dataset's coordinate system.
        max_cell_size : float, optional
            Maximum cell size in metres (unused — resolution is fixed by the
            DataArray).  Default is ``1000.0``.
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
        z = data.values[0, :, :]
        z[z == rds.rio.nodata] = np.nan

        rds.close()

        return x, y, z
