# -*- coding: utf-8 -*-
"""
This module defines the BathymetryDataset class, which represents a bathymetry dataset.
It provides methods to read metadata and get data from the dataset.

Classes:
    BathymetryDataset: A class for handling bathymetry datasets.

Usage:
    from .dataset import BathymetryDataset
"""

import os
import toml
from pyproj import CRS
from typing import Any, List, Optional

class BathymetryDataset:
    """
    Bathymetry dataset class
    """

    def __init__(self):
        """
        Initialize the BathymetryDataset class.
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
        Read metadata file and set attributes.

        Raises:
        FileNotFoundError: If the metadata file does not exist.
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
        Placeholder method to get data from the dataset.
        """
        pass

    def get_bbox(self, **kwargs) -> None:
        """
        Placeholder method to get the bounding box of the dataset.
        """
        pass