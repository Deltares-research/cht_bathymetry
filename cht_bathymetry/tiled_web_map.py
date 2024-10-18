# -*- coding: utf-8 -*-
"""
Created on Sun Apr 25 10:58:08 2021

@author: Maarten van Ormondt
"""

from cht_tiling import TiledWebMap

from .dataset import BathymetryDataset

class BathymetryDatasetTiledWebMap(BathymetryDataset):
    """
    Bathymetry dataset class 

    :ivar name: initial value: ''
    :ivar nr_zoom_levels: initial value: 0
    """

    def __init__(self, name, path):
        super().__init__()
        
        self.name              = name
        self.path              = path
        self.local_path        = path
        self.read_metadata()
        self.data = TiledWebMap(self.local_path, name, parameter="elevation")
            
    def get_data(self, xl, yl, max_cell_size):
        """
        Reads data from database. Returns x, y, z in same coordinate system (3857) as dataset. Resolution is determined by max_cell_size.
        """
        x, y, z = self.data.get_data(xl, yl, max_cell_size)
        
        return x, y, z
