# -*- coding: utf-8 -*-
"""
Created on Sun Apr 25 10:58:08 2021

@author: Maarten van Ormondt
"""

import os
# import xml.etree.ElementTree as ET
import numpy as np
# import math
import urllib
import netCDF4 as nc
import yaml
from shutil import copyfile
# from pyproj import CRS
# from pyproj import Transformer
# from cht_utils.misc_tools import interp2

from .dataset import BathymetryDataset

class ZoomLevel:
    def __init__(self):        
        self.dx = 0.0
        self.dy = 0.0
        self.i_available = []
        self.j_available = []

class BathymetryDatasetNetCDFTilesV2(BathymetryDataset):
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
        self.use_cache         = True
        self.read_metadata()
        self.read_tile_structure()

    def read_tile_structure(self):
        # Read netcdf file with dimensions
        nc_file = os.path.join(self.local_path, self.name + ".nc")
        ds   = nc.Dataset(nc_file)
        self.pixels_in_tile           = ds['tile_size_x'][0]
        self.nr_zoom_levels           = ds.dimensions['zoom_levels'].size
        for izoom in range(self.nr_zoom_levels):
            zl = ZoomLevel()
            zl.x0 = ds['x0'][izoom]
            zl.y0 = ds['y0'][izoom]
            zl.dx = ds['grid_size_x'][izoom]
            zl.dy = ds['grid_size_y'][izoom]
            zl.nr_tiles_x = ds['nr_tiles_x'][izoom]
            zl.nr_tiles_y = ds['nr_tiles_y'][izoom]
            self.zoom_level.append(zl)
            
    def get_data(self, xl, yl, max_cell_size, waitbox=None):
        """
        Reads data from database. Returns x, y, z in same coordinate system as dataset. Resolution is determined by max_cell_size.
        """

        izoom = 0
        
        just_get_tiles = False
        iopendap       = False
                  
        x = np.nan
        y = np.nan
        z = np.nan

        cell_size_x = np.array([])
        cell_size_y = np.array([])
        
        for zl in self.zoom_level:
            cell_size_x = np.append(cell_size_x, zl.dx)
            cell_size_y = np.append(cell_size_y, zl.dy)

            # Should multiply with unit here...
        
        if izoom == 0:
            # Find zoom level based on resolution
            if self.crs.is_geographic:
                cell_size_x = cell_size_x*111111
                cell_size_y = cell_size_y*111111
            # Find first level with cell size greater than max cell size    
            ilev1 = find_last(cell_size_x<=max_cell_size)
            if ilev1 is None:
                # All levels have cell sizes smaller than max cell size
                ilev1 = 0
            ilev2 = ilev1
        elif izoom == -1:
            # Get all the data from each zoom level !
            just_get_tiles = True
            ilev1 = 0
            ilev2 = self.nr_zoom_levels
        else:
            ilev1 = izoom
            ilev2 = izoom
        
        for ilev in range(ilev1,ilev2+1):
            
            x0  = float(self.zoom_level[ilev].x0)   # lower-left corner x
            y0  = float(self.zoom_level[ilev].y0)   # lower-left corner y
            dx  = float(self.zoom_level[ilev].dx)   # cell size x
            dy  = float(self.zoom_level[ilev].dy)   # cell size y
            nx  = self.pixels_in_tile        # number of pixels in tile x
            if self.name == "rws_vaklodingen":
                ny = 625
            else:
                ny  = self.pixels_in_tile        # number of pixels in tile y
            nnx = self.zoom_level[ilev].nr_tiles_x  # number of tiles x
            nny = self.zoom_level[ilev].nr_tiles_y  # number of tiles y

            iav = self.zoom_level[ilev].i_available
            if not np.any(iav):
                ncfile = os.path.join(self.local_path, self.name + ".nc")
                ds  = nc.Dataset(ncfile)
                iav = ds["available_zl" + str(ilev + 1).zfill(2)][:].data
                self.zoom_level[ilev].i_available = iav

#            vert_unit = self.vertical_coordinate_system.unit
            
            tile_size_x = dx*nx
            tile_size_y = dy*ny
            
             # Directories and names
            name   = self.name;
            levdir = 'zl' + str(ilev + 1).zfill(2)
            
            iopendap = False
            ipdrive = False

            if self.url[0:4] == 'http':
                # Tiles stored on OpenDAP server
                iopendap  = True
                remotedir = self.url + '/' + levdir + '/'
                localdir  = os.path.join(self.local_path, levdir)
            elif self.url[0:2].lower() == 'p:':
                ipdrive = True
                remotedir = self.url + '\\' + levdir + '\\'
                localdir  = os.path.join(self.local_path, levdir)
            else:
                # Tiles are stored locally
                localdir  = os.path.join(self.local_path, levdir)
                remotedir = localdir
                
            # Tiles
            # Array with x and y origin of available tiles
            # all_tiles_x0 = np.arange(x0, x0 + (nnx)*tile_size_x, tile_size_x)
            # all_tiles_y0 = np.arange(y0, y0 + (nny)*tile_size_y, tile_size_y)
            all_tiles_x0 = np.linspace(x0, x0 + (nnx - 1)*tile_size_x, num=nnx)
            all_tiles_y0 = np.linspace(y0, y0 + (nny - 1)*tile_size_y, num=nny)
            
            # Make sure that tiles are read east +180 deg lon.
            all_tiles_index_x = np.arange(0, nnx)
            all_tiles_index_y = np.arange(0, nny)
            
            if self.coord_ref_sys_kind == 'geographic' and nnx*tile_size_x>350.0:
                # Probably a global dataset
                all_tiles_x0 = np.concatenate((all_tiles_x0 - 360.0,
                                               all_tiles_x0,
                                               all_tiles_x0 + 360.0))
                all_tiles_index_x = np.concatenate((all_tiles_index_x,
                                                    all_tiles_index_x,
                                                    all_tiles_index_x))
            
            # Required tile indices
            if all_tiles_x0[0]>xl[1] or all_tiles_y0[0]>yl[1] or all_tiles_x0[-1] + tile_size_x<xl[0] or all_tiles_y0[-1] + tile_size_y<yl[0]:
                ok = False
                print('Tiles are outside of search range')
                return x, y, z
             
            ix1 = find_last(all_tiles_x0<=xl[0])
            if ix1 == None:
                ix1 = 0
            ix2 = find_last(all_tiles_x0<xl[1])
            iy1 = find_last(all_tiles_y0<=yl[0])
            if iy1 == None:
                iy1 = 0
            iy2 = find_last(all_tiles_y0<yl[1])
          
            # Total number of tiles to read in x and y direction
            nnnx = ix2 - ix1 + 1
            nnny = iy2 - iy1 + 1

            # Indices of tiles to be loaded
            tiles_index_x = all_tiles_index_x[ix1:ix2 + 1]
            tiles_index_y = all_tiles_index_y[iy1:iy2 + 1]
            # Origins of tiles to be loaded
            tiles_x0 = all_tiles_x0[ix1:ix2 + 1]
            tiles_y0 = all_tiles_y0[iy1:iy2 + 1]
            tiles_x1 = tiles_x0 + tile_size_x
            npixx = int(np.round((tiles_x1[-1] - tiles_x0[0])/dx))

            if not just_get_tiles:
                # Mesh of horizontal coordinates
                # x = np.arange(tiles_x0[0],
                #                tiles_x0[-1] + tile_size_x, dx)
                # y = np.arange(tiles_y0[0],
                #                tiles_y0[-1] + tile_size_y, dy)

                x = np.linspace(tiles_x0[0], tiles_x0[0] + (npixx - 1) * dx, num=npixx)
#                x = np.linspace(tiles_x0[0], tiles_x0[-1] + tile_size_x - dx, num=nnnx*nx)
                y = np.linspace(tiles_y0[0], tiles_y0[-1] + tile_size_y - dy, num=nnny*ny)

                # Allocate z
                z = np.empty((nnny*ny, npixx))
                z[:] = np.nan

            # Start indices for each tile in larger matrix
            istartx = []
            for i in range(nnnx):
                iii1 = find_first(abs(x - tiles_x0[i]) == min(abs(x - tiles_x0[i])))
                istartx.append(iii1)

            tilen = 0 # Tile number index (only used for waitbox)
            ntiles = nnnx*(iy2 - iy1 + 1) # Total number of tiles
            
#             if not quiet
#                 wb = awaitbar(0,'Getting tiles ...')
            
            # Now get the tiles
            for i in range(nnnx):
                
                itile = tiles_index_x[i]
                
                for j in range(nnny):

                    jtile = tiles_index_y[j]
    
                    tilen = tilen + 1

                    zzz = np.empty((ny,nx)) # make empty tile
                    zzz[:] = np.nan

                    # First check whether required file exists at at all

                    file_name = name + '.' + levdir + '.' + str(itile + 1).zfill(5) + '.' + str(jtile + 1).zfill(5) + '.nc'

                    tile_exists = False
                    if self.type == "netcdf_tiles_v2":
                        if iav[jtile, itile] == 1:
                            tile_exists = True
                            idirname = os.path.join(localdir, str(itile + 1).zfill(5))
                            full_file_name = os.path.join(idirname, file_name)
                            var_str = "value"
                    else:      
                        both_ok = (iav == itile) * (jav == jtile)
                        if both_ok.any():
                            tile_exists = True
                            full_file_name = os.path.join(localdir, file_name)
                            var_str = "depth"

                    if tile_exists:
    
                        if iopendap:
                            if self.use_cache:
                                # First check if file is available locally
                                idownload = False
                                if not os.path.exists(full_file_name):
                                    # File not available locally
                                    idownload = True
                                else:
                                    # Check if the file size seems right
                                    fsize = os.path.getsize(full_file_name)
                                    if fsize<1000:
                                        # Probably something wrong with
                                        # this file. Delete it and download
                                        # again.
                                        idownload = True
                                        os.remove(full_file_name);
    
                                if idownload:
                                    # Make localdir if it does not yet exist
                                    if not os.path.exists(localdir):
                                        os.mkdir(localdir)
                                    # Download file    
                                    try:
                                        print("Downloading tile " + file_name + " ...") 
                                        urllib.request.urlretrieve(remotedir + file_name,
                                                                       full_file_name)
                                    except:
                                        print('Could not download tile ...')
    
                                ncfile = full_file_name # name of local netcdf file
    
                            else:
    
                                # Don't use cache
                                ncfile = remotedir + file_name
                                
                        elif ipdrive:
                            if self.use_cache:
                                # First check if file is available locally
                                icopy = False
                                if not os.path.exists(full_file_name):
                                    # File not available locally
                                    icopy = True
    
                                if icopy:
                                    # Make localdir if it does not yet exist
                                    if not os.path.exists(localdir):
                                        os.mkdir(localdir)
                                    # Download file    
                                    try:
                                        copyfile(os.path.join(remotedir,file_name),
                                                                       full_file_name)
                                    except:
                                        print('Could not copy tile ...')
    
                                ncfile = full_file_name # name of local netcdf file
    
                            else:
    
                                # Don't use cache
                                ncfile = remotedir + file_name

                        else:
    
                            ncfile = full_file_name
                            
#                        print(ncfile)    
                        if not just_get_tiles:
                            
                            # Read the data in the tile
                            if os.path.exists(ncfile):
                                ds  = nc.Dataset(ncfile)
                                zzz = ds[var_str][:]
    #                            fill_value = nc_attget(ncfile, 'depth', 'fill_value')
                        
                            # Now stick the tile data in the large array

                            i1 = istartx[i]
                            i2 = istartx[i] + nx

                            j1 = j*ny
                            j2 = (j + 1)*ny

                            z[j1:j2,i1:i2] = zzz
                    
            
#             # # Close waitbar
#             # if ~quiet
#             #     if ~isempty(hh)
#             #         close(wb);
#             #     end
#             # end

            # Now crop the data to the requested limits
            if not just_get_tiles:
                
#                z(z<-15000)=NaN # Set values to NaN

                ix1 = find_last(x<=xl[0])
                if ix1 == None:
                    ix1 = 0
                ix2 = find_first(x>xl[1])
                if ix2 == None:
                    ix2 = len(x)
                iy1 = find_last(y<=yl[0])
                if iy1 == None:
                    iy1 = 0
                iy2 = find_first(y>yl[1])
                if iy2 == None:
                    iy2 = len(y)
                
                x = x[ix1:ix2]
                y = y[iy1:iy2]
                
                data_in_cell_centres = True

                if data_in_cell_centres:
                    x = x + 0.5*dx # This really should be added as the data are defined in the cell centres!!!               
                    y = y + 0.5*dy
                
                z = z[iy1:iy2, ix1:ix2]
                
                # Convert to metres
                if self.vertical_units == 'cm':
                    z = z*0.01
                elif self.vertical_units == 'ft':
                    z = z*0.3048

        return x, y, z
           
def find_first(a):
    if not a.any():
        i = None
    else:
        i = np.nonzero(a)[0][0]
    return i

def find_last(a):
    if not a.any():
        i = None
    else:
        i = np.nonzero(a)[0][-1]
    return i

def dict2yaml(file_name, dct, sort_keys=False):
    yaml_string = yaml.dump(dct, sort_keys=sort_keys)    
    file = open(file_name, "w")  
    file.write(yaml_string)
    file.close()

def yaml2dict(file_name):
    file = open(file_name,"r")
    dct = yaml.load(file, Loader=yaml.FullLoader)
    return dct
