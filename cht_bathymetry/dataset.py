# -*- coding: utf-8 -*-
"""
Created on Sun Apr 25 10:58:08 2021

@author: Maarten van Ormondt
"""
import os
import toml
from pyproj import CRS

class BathymetryDataset:
    """
    Bathymetry dataset class 
    """
    def __init__(self):        
        self.database          = None
        self.name              = ""
        self.long_name         = ""
        self.source            = ""
        self.data_format       = "" # netcdftiles, geotiff, etc
        self.nr_zoom_levels    = 0
        self.zoom_level        = []
        self.coordinate_system = []
        self.use_cache         = True
        self.remote_path       = ""
        self.path              = ""
        self.local_path        = ""
        self.vertical_units    = "m"
        self.vertical_reference_level_name = "MSL"
        self.vertical_reference_level_difference_with_MSL = 0.0
        self.crs               = None

    def read_metadata(self):
        # Read metadata file
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
            
    def get_data(self):
        pass

    def get_bbox(self, **kwargs):
        pass