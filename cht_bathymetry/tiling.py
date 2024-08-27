# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 13:48:55 2022

@author: ormondt
"""
import os
import rasterio
import xarray as xr
import numpy as np
import math
import wget
import yaml

from pyproj import CRS

import cht.misc.fileops as fo
import cht.misc.xmlkit as xml

class TiledDatasetLevel:
    def __init__(self, name, path, ilev):

        self.name = name
        self.path = os.path.join(path, "zl" + str(ilev + 1).zfill(2))
        self.available = []            

    def count_available_tiles(self):
        # Check available tiles
        folders = fo.list_folders(os.path.join(self.path,'*'), basename=True)
        self.available = np.zeros((self.ntilesy, self.ntilesx), dtype=np.int8)
        for ifolder in folders:
            i = int(ifolder) - 1
            files = fo.list_files(os.path.join(self.path, ifolder, "*.nc"))
            for file in files:
                j = int(file[-8:-3]) - 1
                self.available[j, i] = 1
        

class TiledDataset:
    def __init__(self):
        pass

    def get(self):
        pass
    
    def build(self,
              name,
              source,
              database_path,
              tfw=False,
              vrtfile=None,
              url=None,
              tile_size=256,
              make_base_tiles=True,
              make_high_level_tiles=True,
              first_high_level=2,
              minlev=0,
              keep_intermediate=False,
              delete_download=False,
              crs=None,
              skip_file=None,
              attrs={}):

        self.name      = name
        self.source        = source
        self.database_path = database_path

        self.tfw        = tfw
        self.tile_size  = tile_size
        self.url        = url
        self.delete_download = delete_download
        self.temp_path  = None
        self.path       = os.path.join(self.database_path, self.name)
        self.level      = []
        self.attributes = attrs

        first_high_level -= 1

        if crs:
            if isinstance(crs, str) or isinstance(crs, int):
                self.crs = CRS(crs)
            else:
                self.crs = crs
        else:
            self.crs = crs
            

        # Check what sort of file the source is
        if isinstance(source, list):
            # List with file names
            self.get_dimensions_from_list()
        else:
            # Must be a file name, get the extension           
            ext = os.path.splitext(source)[-1]
            if ext=="vrt":
                self.read_vrt()
            else:
                # File with list of source files
                self.source = []
                fid = open(source, 'r')
                lines = fid.readlines()
                self.source = []
                for line in lines:
                    self.source.append(line.strip())
                self.get_dimensions_from_list()

        if make_high_level_tiles:
            if skip_file:
                # Read file with tifs to skip
                if os.path.exists(skip_file):
                    with open(skip_file, 'r') as f:
                        skip_list = f.readlines()
                    f.close()
                    for skip in skip_list:
                        skip = skip.strip()
                        if skip in self.source_file_names:
                            index = self.source_file_names.index(skip)
                            self.source_file_names.pop(index)
                            self.no_data.pop(index)
                else:
                    print("Warning! Skip file does not exist.")
            

        # Determine the dimensions of each level        
        self.level_dimensions()

        if make_base_tiles:
            self.make_base_tiles()
            self.level[0].count_available_tiles()
            
        if make_high_level_tiles:
            self.make_high_level_tiles(first_level=first_high_level)
            
        for level in self.level:
            level.count_available_tiles()
            
        self.write_metadata_netcdf()
        self.write_metadata_yml()
            
            
    def level_dimensions(self):
        # Base tiles
        ilev = 0
        
        ntilesx = math.ceil(self.nx / self.tile_size)
        ntilesy = math.ceil(self.ny / self.tile_size)
        dx      = self.dx
        dy      = self.dy
        
        level = TiledDatasetLevel(self.name, self.path, ilev)
        level.ntilesx = ntilesx
        level.ntilesy = ntilesy
        level.dx      = dx
        level.dy      = dy
        level.ilev    = ilev
        
        self.level.append(level)

        while True:
            ilev += 1
            dx = dx*2
            dy = dy*2
            ntilesx = math.ceil(ntilesx / 2)
            ntilesy = math.ceil(ntilesy / 2)
            level = TiledDatasetLevel(self.name, self.path, ilev)
            level.ntilesx = ntilesx
            level.ntilesy = ntilesy
            level.dx      = dx
            level.dy      = dy
            self.level.append(level)
            if ntilesx <= 1 or ntilesy<=1:
                break
                
    def read_vrt(self):
        
        # Read vrt xml file
        xml_obj = xml.xml2obj(self.vrt_file)
        
        raster_vrt = rasterio.open(self.vrt_file)
        crs = raster_vrt.crs
        if not self.crs:
            epsg = crs.to_epsg()
            self.crs = CRS(epsg)
        
        res = [float(value) for value in xml_obj.GeoTransform[0].value.split(', ')]
        self.flipud = False

        if res[5]<0.0:
            self.flipud = True
        
        self.x0 = raster_vrt.bounds.left
        self.y0 = raster_vrt.bounds.bottom
        self.x1 = raster_vrt.bounds.right
        self.y1 = raster_vrt.bounds.top
        self.nx = int(xml_obj.rasterXSize)
        self.ny = int(xml_obj.rasterYSize)
        self.dx = (self.x1 - self.x0) / self.nx
        self.dy = (self.y1 - self.y0) / self.ny

        self.source_file_names = []
        self.no_data = []
        for src in xml_obj.VRTRasterBand[0].ComplexSource:
            if src.SourceFilename[0].value[0] == "." and src.SourceFilename[0].value[1] == "/":
                    src.SourceFilename[0].value = src.SourceFilename[0].value[2:]
            self.source_file_names.append(src.SourceFilename[0].value)            
            self.no_data.append(float(src.NODATA[0].value))

    def get_dimensions_from_list(self):
        
        self.source_file_names = []
        self.no_data           = []

        self.x0 = 1.0e9
        self.x1 = -1.0e9
        self.y0 = 1.0e9
        self.y1 = -1.0e9

        self.flipud = True
        
        for file in self.source:
            self.source_file_names.append(file)
            dataset = rasterio.open(file)
            self.crs = CRS(dataset.crs.to_epsg())
            self.no_data.append(dataset.nodata)
            left   = dataset.bounds[0]
            bottom = dataset.bounds[1]
            right  = dataset.bounds[2]
            top    = dataset.bounds[3]
            nx     = dataset.width
            ny     = dataset.height
            dx     = (right - left) / nx
            dy     = (top - bottom) / ny
            self.x0 = min(self.x0, left)
            self.x1 = max(self.x1, right)
            self.y0 = min(self.y0, bottom)
            self.y1 = max(self.y1, top)
            # if self.tfw:
            #     name = os.path.splitext(file)[-2]
            #     tfw_file = name + ".ext"
            #     tfw = np.loadtxt(tfw_file)
            #     dx = tfw[0]
            #     dy = tfw[3]
            #     x0 = tfw[4]
            #     y0 = tfw[5]
        self.nx = int(round((self.x1 - self.x0) / dx))
        self.ny = int(round((self.y1 - self.y0) / dy))
        self.dx = dx
        self.dy = dy
            

    def make_base_tiles(self, skip_list=None):
        
        print("Making base tiles Level 1...")

        ilev = 0
                
        level = self.level[0]

        fo.mkdir(level.path)
        
        if self.temp_path:
            fo.mkdir(self.temp_path)
                
        for isrc, source_file_name in enumerate(self.source_file_names):
            
            print("Processing " + source_file_name + " ...")
            no_data = self.no_data[isrc]

            # Download data
            if self.url:
                if not fo.exists(os.path.basename(source_file_name)):
                    print("Downloading ...")                
                    wget.download(self.url + source_file_name)
            
            source_file_name = os.path.basename(source_file_name)  

            if not os.path.exists(source_file_name):
                continue

            dataset = rasterio.open(source_file_name)
            z       = dataset.read(1)
            if self.flipud:
                z = np.flipud(z)
            x0d     = dataset.bounds.left
            y0d     = dataset.bounds.bottom
            x1d     = dataset.bounds.right
            y1d     = dataset.bounds.top
            dataset.close()
            
            ixd0 = round((x0d - self.x0)/level.dx) # first index of total dataset in x direction
            iyd0 = round((y0d - self.y0)/level.dy) # first index of total dataset in y direction
            ixd1 = round((x1d - self.x0)/level.dx) - 1 # last index of total data in x direction
            iyd1 = round((y1d - self.y0)/level.dy) - 1 # last index of total data in y direction

            ixb0 = int(ixd0 / self.tile_size)      # first tile in x direction
            iyb0 = int(iyd0 / self.tile_size)      # first tile in y direction
            ixb1 = int(ixd1 / self.tile_size)      # last tile in x direction
            iyb1 = int(iyd1 / self.tile_size)      # last tile in y direction
            
            zb    = np.empty((self.tile_size, self.tile_size), dtype=np.float32)
            zb[:] = np.nan
            
            iprogress = 0
            ntiles = (ixb1 + 1 - ixb0) * (iyb1 + 1 - iyb0)
            
            # Loop through tiles
            for ixb in range(ixb0, ixb1 + 1):
                
                istr = str(ixb + 1).zfill(5)
                ipath = os.path.join(level.path, istr)
                ipathexists = fo.exists(ipath)
                
                for iyb in range(iyb0, iyb1 + 1):
                    
                    iprogress += 1
                    
                    progress = 100 * iprogress / ntiles
                    
                    tile_file_name = self.name + "." + "zl" + str(ilev + 1).zfill(2) + "." \
                        + str(ixb + 1).zfill(5) + "." + str(iyb + 1).zfill(5) + ".nc"

                    print("Processing " + str(round(progress, 1)) + " - "  + tile_file_name)
                    
                    # Get indices of tile in complete dataset 
                    ixt0 = ixb * self.tile_size
                    ixt1 = (ixb + 1) * self.tile_size - 1
                    iyt0 = iyb * self.tile_size
                    iyt1 = (iyb + 1) * self.tile_size - 1
        
                    # Indices in total dataset for actual data in this tile (possibly just a subset)
                    ixt0r = max(ixt0, ixd0)            
                    iyt0r = max(iyt0, iyd0)            
                    ixt1r = min(ixt1, ixd1)            
                    iyt1r = min(iyt1, iyd1)            
                    
                    # Indices in tile with actual data
                    ixt0t = ixt0r - ixt0
                    iyt0t = iyt0r - iyt0
                    ixt1t = ixt1r - ixt0
                    iyt1t = iyt1r - iyt0

                    # Indices of data block (e.g. geotiff file) with actual data        
                    ixb0r = ixt0r - ixd0
                    iyb0r = iyt0r - iyd0
                    ixb1r = ixt1r - ixd0
                    iyb1r = iyt1r - iyd0
                    
                    zb[:] = np.nan
                    
                    zb[iyt0t:iyt1t + 1, ixt0t:ixt1t + 1] = z[iyb0r:iyb1r + 1, ixb0r:ixb1r + 1]
                    zb[np.where(zb==no_data)]=np.nan

                    zb[np.where(zb<-999999.0)]=np.nan
                    zb[np.where(zb>999999.0)]=np.nan
                    
                    # Check for non-nan values
                    if np.isnan(zb).all():
                        continue
        
                    # Check if this tile already exists
                    
                    if fo.exists(os.path.join(ipath, tile_file_name)):
                        isnan = np.where(np.isnan(zb))
                        # If so, read and merge
                        ds0 = xr.open_dataset(os.path.join(ipath, tile_file_name))
                        z0  = ds0["value"].values
                        ds0.close()
                        z0[np.where(zb==no_data)] = np.nan
                        zb[isnan] = z0[isnan]                        
                    
                    zb[np.where(zb==np.nan)] = no_data
                    
                    x = np.arange(self.x0 + 0.5*level.dx + ixt0 * level.dx, self.x0 + 0.5*level.dx + (ixt1 + 1) * level.dx - 1.0e-9, level.dx)
                    y = np.arange(self.y0 + 0.5*level.dy + iyt0 * level.dy, self.y0 + 0.5*level.dy + (iyt1 + 1) * level.dy - 1.0e-9, level.dy)
        
                    # Make xarray
                    dic = {}
                    dic["value"] = xr.DataArray(data=zb,
                                                dims=["y", "x"],
                                                coords={"y": y, "x": x},
                                                attrs={"_FillValue": np.float32(no_data)})
                    
                    ds = xr.Dataset(dic, attrs=self.attributes)
       
                    # Write to nc
                    if not ipathexists:
                        fo.mkdir(ipath)         
                        ipathexists = True
                    
                    fname_full = os.path.join(ipath, tile_file_name)

                    ds.to_netcdf(path=fname_full, mode="w")
                    

            if self.url and self.delete_download:
                os.remove(source_file_name)     
                
            # if self.temp_path:
            #     print("Moving files to remote folder ...")
            #     fo.move_file(os.path.join(self.temp_path, "*.nc"), level.path)

    def make_high_level_tiles(self, first_level=1):

        print("Making high level tiles ...")
        
        no_data = -99999.0

        for ilev, level in enumerate(self.level):
                        
            if ilev<first_level:
                # Skip level 0
                continue

            print("Level " + str(ilev + 1) + " ...")
            fo.mkdir(level.path)
            
            level0 = self.level[ilev - 1]

            # i=x
            # j=y
            
            # Get available tiles in higher level            
            level0.count_available_tiles()
            
            iimin = np.max(level0.available, axis=0)
            jjmin = np.max(level0.available, axis=1)
            
            iok = np.where(iimin==1)
            jok = np.where(jjmin==1)

            iimin = int(np.min(iok)/2)
            iimax = int(np.max(iok)/2)
            jjmin = int(np.min(jok)/2)
            jjmax = int(np.max(jok)/2)

            ni0 = np.size(iok)
            nj0 = np.size(jok)

            iprogress = 0
            ntiles = (iimax - iimin + 1) * (jjmax - jjmin + 1)

            for ib in range(iimin, iimax + 1):

                istr = str(ib + 1).zfill(5)
                ipath = os.path.join(level.path, istr)
                ipathexists = fo.exists(ipath)

                print(str(round(100*iprogress/ntiles, 1)) + "%\r")

                for jb in range(jjmin, jjmax + 1):
                    
                    iprogress += 1
                    
                    
                    okay = False
                    
                    zzz = np.empty((self.tile_size*2, self.tile_size*2), dtype=np.float32)
                    zzz[:] = np.nan

                    lstr = str(ilev).zfill(2)
                    
                    # Lower-left
                    iii = ib*2 + 0
                    jjj = jb*2 + 0
                    if iii<ni0 and jjj<nj0:
                        istr = str(iii + 1).zfill(5)
                        jstr = str(jjj + 1).zfill(5)
                        if level0.available[jjj, iii] == 1:
                            okay = True
                            fname = self.name + "." + "zl" + lstr + "." \
                                + istr + "." + jstr + ".nc"
                            fname = os.path.join(level0.path, istr, fname)    
                            ds0   = xr.open_dataset(fname)
                            zzz[0:self.tile_size, 0:self.tile_size] = ds0["value"].values


                    # Upper-left
                    iii = ib*2 + 0
                    jjj = jb*2 + 1
                    istr = str(iii + 1).zfill(5)
                    jstr = str(jjj + 1).zfill(5)
                    if iii<ni0 and jjj<nj0:
                        if level0.available[jjj, iii] == 1:
                            okay = True
                            fname = self.name + "." + "zl" + lstr + "." \
                                + istr + "." + jstr + ".nc"
                            fname = os.path.join(level0.path, istr, fname)    
                            ds0   = xr.open_dataset(fname)
                            zzz[self.tile_size:2*self.tile_size, 0:self.tile_size] = ds0["value"].values

                    # Lower-right
                    iii = ib*2 + 1
                    jjj = jb*2 + 0
                    istr = str(iii + 1).zfill(5)
                    jstr = str(jjj + 1).zfill(5)
                    if iii<ni0 and jjj<nj0:
                        if level0.available[jjj, iii] == 1:
                            okay = True
                            fname = self.name + "." + "zl" + lstr + "." \
                                + istr + "." + jstr + ".nc"
                            fname = os.path.join(level0.path, istr, fname)    
                            ds0   = xr.open_dataset(fname)
                            zzz[0:self.tile_size, self.tile_size:2*self.tile_size] = ds0["value"].values

                    # Upper-right
                    iii = ib*2 + 1
                    jjj = jb*2 + 1
                    istr = str(iii + 1).zfill(5)
                    jstr = str(jjj + 1).zfill(5)
                    if iii<ni0 and jjj<nj0:
                        if level0.available[jjj, iii] == 1:
                            okay = True
                            fname = self.name + "." + "zl" + lstr + "." \
                                + istr + "." + jstr + ".nc"
                            fname = os.path.join(level0.path, istr, fname)    
                            ds0   = xr.open_dataset(fname)
                            zzz[self.tile_size:2*self.tile_size, self.tile_size:2*self.tile_size] = ds0["value"].values

                    if okay:                    

                        zz = np.empty((4, self.tile_size, self.tile_size), dtype=np.float32)
                        zz[:] = np.nan
                        
                        zz[0, :, :] = zzz[0:2*self.tile_size:2, 0:2*self.tile_size:2]
                        zz[1, :, :] = zzz[1:2*self.tile_size:2, 0:2*self.tile_size:2]
                        zz[2, :, :] = zzz[0:2*self.tile_size:2, 1:2*self.tile_size:2]
                        zz[3, :, :] = zzz[1:2*self.tile_size:2, 1:2*self.tile_size:2]
                         
                        zb = np.nanmean(zz, axis=0)
    
                        # Check for non-nan values
                        if np.isnan(zb).all():
                            continue
                                
                        zb[np.where(zb==np.nan)] = no_data
                        
                        x = np.arange(self.x0 + 0.5*level.dx + ib * self.tile_size * level.dx, \
                                      self.x0 + 0.5*level.dx + (ib + 1) * self.tile_size * level.dx - 1.0e-9, level.dx)
                        y = np.arange(self.y0 + 0.5*level.dx + jb * self.tile_size * level.dy, \
                                      self.y0 + 0.5*level.dx + (jb + 1) * self.tile_size * level.dy - 1.0e-9, level.dy)

                        # Make xarray
                        dic = {}
                        dic["value"] = xr.DataArray(data=zb,
                                                    dims=["y", "x"],
                                                    coords={"y": y, "x": x},
                                                    attrs={"_FillValue": np.float32(no_data)})
                        
                        ds = xr.Dataset(dic, attrs=self.attributes)
            
                        # Write to nc
                        lstr = str(ilev + 1).zfill(2)
                        istr = str(ib + 1).zfill(5)
                        jstr = str(jb + 1).zfill(5)
                        ipath = os.path.join(level.path, istr)
                        fname = self.name + "." + "zl" + lstr + "." \
                            + istr + "." + jstr + ".nc"
                        tile_file_name = os.path.join(level.path, istr, fname)    
                        if not ipathexists:
                            fo.mkdir(ipath)
                            ipathexists = True
                        
                        ds.to_netcdf(path=tile_file_name, mode="w")

    def write_metadata_netcdf(self):            
        
        nlevs = len(self.level)
        nc_grid_size_x = np.empty(nlevs, dtype=np.float64)
        nc_grid_size_y = np.empty(nlevs, dtype=np.float64)
        nc_x0          = np.empty(nlevs, dtype=np.float64)
        nc_y0          = np.empty(nlevs, dtype=np.float64)
        nc_nx          = np.empty(nlevs, dtype=np.int32)                      
        nc_ny          = np.empty(nlevs, dtype=np.int32)                      
        nc_ntilesx     = np.empty(nlevs, dtype=np.int32)                      
        nc_ntilesy     = np.empty(nlevs, dtype=np.int32)                      

        for ilev, level in enumerate(self.level):

            nc_grid_size_x[ilev] = level.dx
            nc_grid_size_y[ilev] = level.dy
            nc_x0[ilev]          = self.x0                      
            nc_y0[ilev]          = self.y0        
            nc_nx[ilev]          = self.tile_size                      
            nc_ny[ilev]          = self.tile_size          
            nc_ntilesx[ilev]     = level.ntilesx 
            nc_ntilesy[ilev]     = level.ntilesy     
           
        # Make metadata nc file            
        dic = {}
        if self.crs.is_projected:
            kind_string = "projected"
        else:
            kind_string = "geographic"
            
        dic["crs"] = xr.DataArray(
                                  attrs={"coord_ref_sys_name": self.crs.name,
                                         "coord_ref_sys_kind": kind_string,
                                         "coord_ref_sys_code": self.crs.to_epsg()})
        
        dic["grid_size_x"] = xr.DataArray(data=nc_grid_size_x,
                                          dims=["zoom_levels"]
                                         )
        
        dic["grid_size_y"] = xr.DataArray(data=nc_grid_size_y,
                                          dims=["zoom_levels"]
                                         )
        
        dic["x0"] = xr.DataArray(data=nc_x0,
                                          dims=["zoom_levels"]
                                         )
        
        dic["y0"] = xr.DataArray(data=nc_y0,
                                          dims=["zoom_levels"]
                                         )
        
        dic["tile_size_x"] = xr.DataArray(data=nc_nx,
                                          dims=["zoom_levels"]
                                         )
        
        dic["tile_size_y"] = xr.DataArray(data=nc_ny,
                                          dims=["zoom_levels"]
                                         )
        
        dic["nr_tiles_x"] = xr.DataArray(data=nc_ntilesx,
                                          dims=["zoom_levels"]
                                         )
        
        dic["nr_tiles_y"] = xr.DataArray(data=nc_ntilesy,
                                          dims=["zoom_levels"]
                                         )
        for ilev, level in enumerate(self.level):

            dataname = "available_zl" + str(ilev + 1).zfill(2)
            dimname_i = "available_x_zl" + str(ilev + 1)        
            dimname_j = "available_y_zl" + str(ilev + 1)        
            dic[dataname] = xr.DataArray(data=level.available,
                                         dims=[dimname_j, dimname_i])
        
        ds = xr.Dataset(dic, attrs=self.attributes)
        
        metadata_file_name = os.path.join(self.path, self.name + ".nc")
        
        # Write to nc
        ds.to_netcdf(path=metadata_file_name, mode="w")

    def write_metadata_yml(self):   

        attrs = self.attributes
        
        if self.crs.is_projected:
            kind_string = "projected"
        else:
            kind_string = "geographic"
        
        attrs["coord_ref_sys_name"]       = self.crs.name
        attrs["coord_ref_sys_kind"]       = kind_string
        attrs["coord_ref_sys_code"]       = self.crs.to_epsg()

        
        yaml_file_name = os.path.join(self.path, self.name + ".yml")

        with open(yaml_file_name, 'w') as file:
            yaml.dump(attrs, file, sort_keys=False)
