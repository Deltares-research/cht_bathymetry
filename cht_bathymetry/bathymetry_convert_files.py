import os
import xarray as xr
import numpy as np
from scipy.interpolate import RegularGridInterpolator
import xml.etree.ElementTree as ET
# from cht_sfincs import SFINCS
# from cht_tiling import make_index_tiles, make_topobathy_tiles_v2, make_topobathy_overlay
# from cht_tiling.utils import get_zoom_level_for_resolution

import os
from cht_utils import fileops as fo
import cht_utils.xmlkit as xml
import yaml
import toml

def dict2yaml(file_name, dct, sort_keys=False):
    yaml_string = yaml.dump(dct, sort_keys=sort_keys)    
    file = open(file_name, "w")  
    file.write(yaml_string)
    file.close()

def yaml2dict(file_name):
    file = open(file_name,"r")
    dct = yaml.load(file, Loader=yaml.FullLoader)
    return dct

path = "c:\\work\\delftdashboard\\data\\bathymetry"

xml_file = os.path.join(path, "bathymetry.xml")
xml_root = ET.parse(xml_file).getroot()

dct = {}
dct["dataset"] = []

for xml_dataset in xml_root.findall('dataset'):


    dsdct = {}

    # find name
    for prop in xml_dataset:
        tag = prop.tag.lower()
        txt = prop.text
        if tag == "name":
            name = txt
            print(name)
            dsdct["name"] = name
            break

    ncfile = os.path.join(path, name, name + ".nc")
    if os.path.exists(ncfile) == False:
        continue    


    dct["dataset"].append(dsdct)        
    # dct[name] = {}        

    dsdct = {}

    # Set attributes
    for prop in xml_dataset:
        tag = prop.tag.lower()
        txt = prop.text
        if txt is not None:
            if tag == "name":
                continue
            if tag == "url":
                if txt[0] != "h":
                    continue
            if tag == "edit":
                continue
            if tag == "version":
                continue
            if tag == "usecache":
                continue
            if tag == "type":
                tag = "format"
                if txt.lower() == "netcdftiles":
                    txt = "netcdf_tiles_v1"

            dsdct[tag] = txt
#            print(tag + " : " + txt)

    # check for existing yaml file (netcdfv2)
    yfile = os.path.join(path, name, name + ".yml")
    if os.path.exists(yfile):
        tmpdct = yaml2dict(yfile)
        for key in tmpdct:
            dsdct[key] = tmpdct[key]

    ds = xr.open_dataset(ncfile)

    attrs = ds["crs"].attrs
    for attr in attrs:
        if attr == "difference_with_msl":
            attrs[attr] = float(attrs[attr])
        dsdct[attr] = attrs[attr]

    attrs = ds.attrs
    dsdct["description"] = {}
    for attr in attrs:
        if attr == "Conventions":
            continue
        if attr == "CF:featureType":
            continue
        dsdct["description"][attr.lower()] = attrs[attr]
    xxx=1

    ds.close()

#    dict2yaml(name + ".yml", dsdct, sort_keys=False)
    with open(os.path.join(path, name, name + ".tml"), "w") as f:
        new_toml_string = toml.dump(dsdct, f)            

#dict2yaml("bathymetry.yml", dct, sort_keys=False)
xxx=1

with open(os.path.join(path, "bathymetry.tml"), "w") as f:
    new_toml_string = toml.dump(dct, f)            
#        setattr(dataset, prop.tag.lower(), prop.text)
#         # Set attributes