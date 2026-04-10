"""
cht_bathymetry: Bathymetry database and dataset access library.

Provides :class:`BathymetryDatabase` for managing collections of bathymetry
datasets in various formats (NetCDF tiles, COG, tiled web maps, xarray
DataArrays) and interpolating depth values onto arbitrary grids or point sets.
"""

from .database import BathymetryDatabase  # noqa: F401
