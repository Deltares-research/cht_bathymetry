"""
Base class for all bathymetry dataset types.

Defines the :class:`BathymetryDataset` interface that concrete dataset
implementations (NetCDF tiles, COG, tiled web map, DataArray) inherit from.
Handles metadata loading from TOML files and CRS resolution.
"""

import os
from typing import Any, List, Optional

import toml
from pyproj import CRS


class BathymetryDataset:
    """
    Abstract base class for a bathymetry dataset.

    Subclasses must implement :meth:`get_data`, :meth:`get_bbox`, and
    :meth:`get_lon_lat_range`.
    """

    def __init__(self) -> None:
        """
        Initialise default attributes shared by all dataset types.
        """
        self.database: Any = None
        self.name: str = ""
        self.long_name: str = ""
        self.source: str = ""
        self.data_format: str = ""  # netcdftiles, geotiff, etc
        self.nr_zoom_levels: int = 0
        self.zoom_level: List[int] = []
        self.coordinate_system: List[str] = []
        self.use_cache: bool = True
        self.remote_path: str = ""
        self.path: str = ""
        self.local_path: str = ""
        self.vertical_units: str = "m"
        self.vertical_reference_level_name: str = "MSL"
        self.vertical_reference_level_difference_with_MSL: float = 0.0
        self.crs: Optional[CRS] = None

    def read_metadata(self) -> None:
        """
        Read a TOML metadata file and populate instance attributes.

        Looks first for ``<local_path>/<name>.tml``, then falls back to
        ``<local_path>/metadata.tml``.  Every key in the TOML file is set as
        an attribute via :func:`setattr`.  The ``longname`` key is accepted as
        an alias for ``long_name`` for backwards compatibility.

        Raises
        ------
        FileNotFoundError
            If neither metadata file exists at ``local_path``.
        """
        tml_file = os.path.join(self.local_path, self.name + ".tml")
        if not os.path.exists(tml_file):
            tml_file = os.path.join(self.local_path, "metadata.tml")
        tml = toml.load(tml_file)
        for key in tml:
            setattr(self, key, tml[key])
        # Long name for backwards compatibility
        if "longname" in tml:
            self.long_name = tml["longname"]
        # Make sure there is always a long_name
        if self.long_name == "":
            self.long_name = self.name

        if "coord_ref_sys_name" in tml:
            self.crs = CRS(tml["coord_ref_sys_name"])

    def get_data(self) -> None:
        """
        Return depth data for a bounding box.

        Placeholder — subclasses must override this method.
        """
        pass

    def get_bbox(self, **kwargs) -> None:
        """
        Return the spatial bounding box of the dataset.

        Placeholder — subclasses must override this method.
        """
        pass

    def get_lon_lat_range(self, **kwargs) -> None:
        """
        Return the longitude and latitude range of the dataset.

        Placeholder — subclasses must override this method.
        """
        pass
