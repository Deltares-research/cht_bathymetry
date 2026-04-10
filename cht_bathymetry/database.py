"""
Bathymetry database: top-level registry and interpolation engine.

Provides :class:`BathymetryDatabase`, which loads dataset metadata from a
TOML catalogue file, optionally synchronises with an S3 bucket, and exposes
methods for interpolating depth onto arbitrary grids or point sets.

Helper utilities (:func:`dict2yaml`, :func:`yaml2dict`, :func:`inpolygon`)
are included for YAML I/O and polygon masking.
"""

import os
from typing import Optional

import boto3
import geopandas as gpd
import numpy as np
import rioxarray
import toml
import yaml
from botocore import UNSIGNED
from botocore.client import Config
from cht_utils.interpolation import interp2
from matplotlib import path
from pyproj import Transformer

from .cog import BathymetryDatasetCOG
from .dataarray import BathymetryDatasetDataArray
from .netcdf_tiles_v1 import BathymetryDatasetNetCDFTilesV1
from .netcdf_tiles_v2 import BathymetryDatasetNetCDFTilesV2
from .tiled_web_map import BathymetryDatasetTiledWebMap


class BathymetryDatabase:
    """
    Registry of bathymetry datasets with grid/point interpolation support.

    Reads a ``bathymetry.tml`` catalogue from *path*, instantiates the
    appropriate dataset class for each entry, and provides methods to
    interpolate depth data onto user-supplied grids or point clouds.

    Parameters
    ----------
    path : str or None, optional
        Local directory that contains ``bathymetry.tml`` and subdirectories
        for each dataset.
    s3_bucket : str or None, optional
        S3 bucket name used for remote dataset synchronisation.
    s3_key : str or None, optional
        S3 key prefix under which datasets are stored.
    s3_region : str or None, optional
        AWS region of the S3 bucket.
    check_online : bool, optional
        When ``True``, synchronise the local catalogue against the S3 bucket
        on construction.  Default is ``False``.
    """

    def __init__(
        self,
        path: Optional[str] = None,
        s3_bucket: Optional[str] = None,
        s3_key: Optional[str] = None,
        s3_region: Optional[str] = None,
        check_online: bool = False,
    ) -> None:
        self.dataset = []
        self.s3_client = None
        self.s3_bucket = s3_bucket
        self.s3_key = s3_key
        self.s3_region = s3_region
        self.path = path
        self.dataset = []
        self.read()
        if check_online:
            self.check_online_database()
        self.initialized = True

    def read(self) -> None:
        """
        Load metadata for all datasets listed in ``bathymetry.tml``.

        Creates the local database directory when it does not exist yet.
        Skips any dataset whose metadata file cannot be found and prints a
        warning.
        """

        if self.path is None:
            print("Path to bathymetry database not set !")
            return

        # Check if the path exists. If not, create it.
        if not os.path.exists(self.path):
            os.makedirs(self.path)

        # Read in database
        tml_file = os.path.join(self.path, "bathymetry.tml")

        if not os.path.exists(tml_file):
            print(f"Warning! Bathymetry database file not found: {tml_file}")
            return

        datasets = toml.load(tml_file)

        for d in datasets["dataset"]:
            name = d["name"]

            if "path" in d:
                path = d["path"]
            else:
                path = os.path.join(self.path, name)

            # Read the meta data for this dataset
            fname = os.path.join(path, f"{name}.tml")

            if os.path.exists(fname):
                metadata = toml.load(fname)
                dataset_format = metadata["format"]
            elif os.path.exists(os.path.join(path, "metadata.tml")):
                # More likely to exist
                metadata = toml.load(os.path.join(path, "metadata.tml"))
                dataset_format = metadata["format"]
            else:
                print(
                    f"Could not find metadata file for dataset {name} ! Skipping dataset."
                )
                continue

            if dataset_format == "netcdf_tiles_v1":
                dataset = BathymetryDatasetNetCDFTilesV1(name, path)
            elif dataset_format == "netcdf_tiles_v2":
                dataset = BathymetryDatasetNetCDFTilesV2(name, path)
            elif dataset_format == "tiled_web_map":
                dataset = BathymetryDatasetTiledWebMap(name, path)
            elif dataset_format == "cog":
                dataset = BathymetryDatasetCOG(name, path)

            dataset.database = self

            self.dataset.append(dataset)

    def load_dataset(self, name: str) -> None:
        """
        Load a single dataset by name and add (or replace) it in the registry.

        Parameters
        ----------
        name : str
            Dataset name; a subdirectory with this name is expected under
            ``self.path``.
        """
        path = os.path.join(self.path, name)
        metadata = toml.load(os.path.join(path, "metadata.tml"))
        dataset_format = metadata["format"]
        if dataset_format == "netcdf_tiles_v1":
            dataset = BathymetryDatasetNetCDFTilesV1(name, path)
        elif dataset_format == "netcdf_tiles_v2":
            dataset = BathymetryDatasetNetCDFTilesV2(name, path)
        elif dataset_format == "tiled_web_map":
            dataset = BathymetryDatasetTiledWebMap(name, path)
        elif dataset_format == "cog":
            dataset = BathymetryDatasetCOG(name, path)
        dataset.database = self
        # Check if dataset already exists in database
        for d in self.dataset:
            if d.name == name:
                # Replace existing dataset
                d = dataset
                return
        self.dataset.append(dataset)

    def add_dataset_from_file(self, name: str, filename: str) -> None:
        """
        Add a dataset to the database by opening a GeoTIFF file.

        Parameters
        ----------
        name : str
            Identifier to assign to the new dataset.
        filename : str
            Path to a GeoTIFF (or any rasterio-readable raster) file.
        """
        da = rioxarray.open_rasterio(filename)
        self.add_dataarray_dataset(da, name)

    def add_dataarray_dataset(self, da, name: str) -> None:
        """
        Wrap an :class:`xarray.DataArray` and register it as a dataset.

        Parameters
        ----------
        da : xr.DataArray
            Source DataArray with a ``rio`` CRS accessor.
        name : str
            Identifier for the dataset.
        """
        dataset = BathymetryDatasetDataArray(da, name)
        # Get the CRS from the DataArray
        dataset.crs = da.rio.crs
        self.dataset.append(dataset)

    def add_datasets_from_files(self, flist: list[str]) -> None:
        """
        Add multiple datasets from a list of GeoTIFF file paths.

        Parameters
        ----------
        flist : list[str]
            Paths to raster files.  Each file is opened with
            :func:`rioxarray.open_rasterio` and registered under its file
            path as the dataset name.
        """
        for f in flist:
            print(f)
            # Add dataset to database
            # Read geotiff as xarray DataArray
            da = rioxarray.open_rasterio(f)
            self.add_dataarray_dataset(da, f)

    def write(self) -> None:
        """
        Persist the database catalogue to disk.

        Currently a no-op placeholder.
        """
        pass

    def check_online_database(self) -> None:
        """
        Synchronise the local catalogue with the S3 bucket.

        Downloads ``bathymetry.tml`` from S3 and compares it with the local
        catalogue.  For each dataset present in S3 but not locally, the
        metadata (and, where applicable, supporting files such as
        ``available_tiles.nc`` or ``index.html``) are downloaded and the
        local catalogue is updated.
        """
        if self.s3_client is None:
            self.s3_client = boto3.client(
                "s3", config=Config(signature_version=UNSIGNED)
            )
        if self.s3_bucket is None:
            return
        # First download a copy of bathymetry.tml and call it bathymetry_s3.tml
        key = f"{self.s3_key}/bathymetry.tml"
        filename = os.path.join(self.path, "bathymetry_s3.tml")
        print("Updating bathymetry database ...")
        try:
            self.s3_client.download_file(
                Bucket=self.s3_bucket,
                Key=key,
                Filename=filename,
            )
        except Exception:
            print(
                f"Failed to download {key} from {self.s3_bucket}. Database will not be updated."
            )
            return

        # Read bathymetry_s3.tml
        short_name_list, long_name_list, source_name_list = self.dataset_names()
        datasets_s3 = toml.load(filename)
        datasets_added = False
        added_names = []
        # Loop through s3 datasets, and check whether they exist in the local database.
        # If so, check if the metadata also exists. If not, make local folder and download the metadata.
        # Additionally, check if available_tiles.nc in s3 and not in local database, download it.

        for d in datasets_s3["dataset"]:
            # Get list of existing datasets
            dataset_name = d["name"]

            if dataset_name not in short_name_list:
                # Dataset not in local database
                print(f"Adding bathymetry dataset {dataset_name} to local database ...")

                # Create folder and download metadata
                path = os.path.join(self.path, dataset_name)
                os.makedirs(path, exist_ok=True)

                key = f"{self.s3_key}/{dataset_name}/metadata.tml"
                filename = os.path.join(path, "metadata.tml")

                # Download metadata
                try:
                    self.s3_client.download_file(
                        Bucket=self.s3_bucket,
                        Key=key,
                        Filename=filename,
                    )
                except Exception as e:
                    print(e)
                    print(f"Failed to download {key}. Skipping dataset.")
                    continue

                # Read metadata.tml (can't use BathymetryDataset.read_metadata(),
                # because it is not a dataset yet)
                metadata = toml.load(filename)
                format = metadata["format"]

                dataset_s3_bucket = metadata["s3_bucket"]
                dataset_s3_key = metadata["s3_key"]

                if format == "tiled_web_map":
                    # For tiled web maps, this is a good idea,
                    # because not much data needs to be downloaded
                    if "available_tiles" in metadata:
                        if metadata["available_tiles"]:
                            # Download available_tiles.nc
                            key = f"{self.s3_key}/{dataset_name}/available_tiles.nc"
                            filename = os.path.join(path, "available_tiles.nc")
                            try:
                                self.s3_client.download_file(
                                    Bucket=self.s3_bucket,
                                    Key=key,
                                    Filename=filename,
                                )
                            except Exception:
                                print(f"Failed to download {key}. Skipping dataset.")
                                continue
                    key = f"{self.s3_key}/{dataset_name}/index.html"
                    filename = os.path.join(path, "index.html")
                    try:
                        self.s3_client.download_file(
                            Bucket=self.s3_bucket,
                            Key=key,
                            Filename=filename,
                        )
                    except Exception:
                        print("index.html could not be downloaded")

                elif format == "cog":
                    # Do not download now, but wait until the data is needed
                    # We'll do it in the COG class
                    pass

                # Necessary data has been downloaded
                datasets_added = True
                added_names.append(dataset_name)

        # Write new local bathymetry.tml
        if datasets_added:
            d = {}
            d["dataset"] = []
            for name in short_name_list:
                d["dataset"].append({"name": name})
            for name in added_names:
                d["dataset"].append({"name": name})
            # Now write the new bathymetry.tml
            with open(os.path.join(self.path, "bathymetry.tml"), "w") as tml:
                toml.dump(d, tml)
            # Read the database again
            self.dataset = []
            self.read()
        else:
            print("No new datasets were added to the local database.")

    def get_bathymetry_on_points(
        self,
        xz: np.ndarray,
        yz: np.ndarray,
        dxmin: float,
        crs,
        bathymetry_list: list,
        method: str = "linear",
    ) -> np.ndarray:
        """
        Interpolate bathymetry onto a cloud of points.

        Parameters
        ----------
        xz : np.ndarray
            1-D array of x coordinates of the target points.
        yz : np.ndarray
            1-D array of y coordinates of the target points.
        dxmin : float
            Minimum spacing between points, used as the resolution hint.
        crs : pyproj.CRS
            CRS of the target coordinates.
        bathymetry_list : list[dict]
            Ordered list of dataset descriptors (each a dict with at least
            ``"name"`` and optionally ``"zmin"``, ``"zmax"``, ``"polygon"``).
        method : str, optional
            Interpolation method passed to :func:`cht_utils.interpolation.interp2`.
            Default is ``"linear"``.

        Returns
        -------
        zz : np.ndarray
            1-D depth array aligned with *xz*/*yz* (NaN where no data).
        """
        zz = self.get_bathymetry_on_grid(
            xz, yz, crs, bathymetry_list, method=method, coords="points", dxmin=dxmin
        )
        return zz

    def get_bathymetry_on_grid(
        self,
        xz: np.ndarray,
        yz: np.ndarray,
        crs,
        bathymetry_list: list,
        method: str = "linear",
        coords: str = "grid",
        dxmin: float = 1.0e6,
        waitbox=None,
        buffer: float = 0.2,
    ) -> np.ndarray:
        """
        Interpolate bathymetry onto a regular grid or arbitrary point set.

        Iterates through *bathymetry_list* in order, filling NaN cells from
        each successive dataset.  Applies optional depth clipping
        (``zmin``/``zmax``) and polygon masking per dataset entry.

        Parameters
        ----------
        xz : np.ndarray
            X coordinates of the target grid or points.  For a 1-D regular
            grid with ``coords="grid"`` this is the x-axis vector; for a 2-D
            curvilinear grid or point set it is a full coordinate array.
        yz : np.ndarray
            Y coordinates (same shape convention as *xz*).
        crs : pyproj.CRS
            CRS of the target coordinates.
        bathymetry_list : list[dict]
            Ordered list of dataset descriptors.
        method : str, optional
            Interpolation method.  Default is ``"linear"``.
        coords : str, optional
            ``"grid"`` (default) or ``"points"``.
        dxmin : float, optional
            Resolution hint used when ``coords="points"``.  Default is
            ``1.0e6``.
        waitbox : optional
            Reserved for a progress-dialog object; currently unused.
        buffer : float, optional
            Fractional buffer added around the bounding box when fetching
            source data.  Default is ``0.2``.

        Returns
        -------
        zz : np.ndarray
            Depth array matching the shape of the target grid/points (NaN
            where no data from any dataset).
        """

        if xz.ndim == 2:
            # xy and yz are a grid
            zz = np.full(xz.shape, np.nan)
            dx = np.sqrt((xz[0, 1] - xz[0, 0]) ** 2 + (yz[0, 1] - yz[0, 0]) ** 2)
            dy = np.sqrt((xz[1, 0] - xz[0, 0]) ** 2 + (yz[1, 0] - yz[0, 0]) ** 2)
        else:
            if coords == "grid":
                zz = np.full((len(yz), len(xz)), np.nan)
                dx = xz[1] - xz[0]
                dy = yz[1] - yz[0]
                xz, yz = np.meshgrid(xz, yz)
            else:
                zz = np.full(xz.shape, np.nan)

        # Determine resolution to get bathy data
        if coords == "grid":
            # Resolution follow from grid
            if crs.is_geographic:
                dx = min(
                    111111.0 * dx * np.cos(np.pi * np.max(np.abs(yz)) / 180.0),
                    111111.0 * dy,
                )
            else:
                dx = min(dx, dy)
        else:
            dx = dxmin

        # Loop through bathymetry datasets
        for ibathy, bathymetry in enumerate(bathymetry_list):
            if "zmin" not in bathymetry:
                bathymetry["zmin"] = -1.0e9
            if "zmax" not in bathymetry:
                bathymetry["zmax"] = 1.0e9

            dataset_name = bathymetry["name"]  # name
            zmin = bathymetry["zmin"]
            zmax = bathymetry["zmax"]

            dataset = self.get_dataset(dataset_name)

            transformer = Transformer.from_crs(crs, dataset.crs, always_xy=True)

            if np.isnan(zz).any():
                xzb, yzb = transformer.transform(xz, yz)
                xmin = np.nanmin(np.nanmin(xzb))
                xmax = np.nanmax(np.nanmax(xzb))
                ymin = np.nanmin(np.nanmin(yzb))
                ymax = np.nanmax(np.nanmax(yzb))
                ddx = buffer * (xmax - xmin)
                ddy = buffer * (ymax - ymin)
                xl = [xmin - ddx, xmax + ddx]
                yl = [ymin - ddy, ymax + ddy]

                # Get DEM data
                xb, yb, zb = dataset.get_data(xl, yl, max_cell_size=dx, waitbox=waitbox)

                # If zb equal np.nan, then there is no data
                if not np.isnan(zb).all():
                    zb[np.where(zb < zmin)] = np.nan
                    zb[np.where(zb > zmax)] = np.nan
                    zz1 = interp2(
                        xb, yb, zb, xzb, yzb, method=method
                    )  # bathymetry from this dataset
                    # Check if a polygon is given
                    if "polygon_file" in bathymetry:
                        # Read the polygon file
                        bathymetry["polygon"] = gpd.read_file(
                            bathymetry["polygon_file"]
                        ).to_crs(dataset.crs)
                    if "polygon" in bathymetry:
                        # Mask out values outside the polygon
                        if isinstance(bathymetry["polygon"], gpd.GeoDataFrame):
                            # Loop through polygons in gdf
                            inpols = np.full(xzb.shape, False)
                            for ip, polygon in bathymetry["polygon"].iterrows():
                                inpol = inpolygon(xzb, yzb, polygon["geometry"])
                                inpols = np.logical_or(inpols, inpol)

                            zz1[~inpols] = np.nan

                    isn = np.where(np.isnan(zz))
                    zz[isn] = zz1[isn]

        return zz

    def get_dataset(self, name: str):
        """
        Retrieve a dataset by name.

        Parameters
        ----------
        name : str
            Dataset identifier.

        Returns
        -------
        dataset : BathymetryDataset or None
            The matching dataset, or ``None`` if not found.
        """
        for dataset in self.dataset:
            if dataset.name == name:
                return dataset
        return None

    def get_lon_lat_range(self, name: str) -> tuple[list[float] | None, list[float] | None]:
        """
        Return the geographic lon/lat range for a named dataset.

        Parameters
        ----------
        name : str
            Dataset identifier.

        Returns
        -------
        lon_range : list[float] or None
        lat_range : list[float] or None
        """
        dataset = self.get_dataset(name)
        if dataset is None:
            print(f"Dataset {name} not found in database.")
            return None, None
        lon_range, lat_range = dataset.get_lon_lat_range()
        return lon_range, lat_range

    def dataset_names(
        self, source: Optional[str] = None
    ) -> tuple[list[str], list[str], list[str]]:
        """
        Return lists of dataset names, optionally filtered by source.

        Parameters
        ----------
        source : str or None, optional
            When given, only datasets whose ``source`` attribute matches are
            included.

        Returns
        -------
        short_name_list : list[str]
            Short (identifier) names.
        long_name_list : list[str]
            Human-readable long names.
        source_name_list : list[str]
            Source labels.
        """
        short_name_list = []
        long_name_list = []
        source_name_list = []
        for dataset in self.dataset:
            ok = False
            if source:
                if dataset.source == source:
                    ok = True
            else:
                ok = True
            if ok:
                short_name_list.append(dataset.name)
                long_name_list.append(dataset.long_name)
                source_name_list.append(dataset.source)
        return short_name_list, long_name_list, source_name_list

    def sources(self) -> tuple[list[str], list]:
        """
        Return all unique sources and their associated datasets.

        Returns
        -------
        source_names : list[str]
            Unique source labels.
        sources : list[BathymetrySource]
            Corresponding :class:`BathymetrySource` objects, each holding the
            datasets that belong to that source.
        """

        sources = []
        source_names = []

        for dataset in self.dataset:
            source = dataset.source
            if source in source_names:
                # Existing source
                for src in sources:
                    if src.name == source:
                        src.dataset.append(dataset)
            else:
                # New source
                src = BathymetrySource(source)
                src.dataset.append(dataset)
                sources.append(src)
                source_names.append(source)

        return source_names, sources


class BathymetrySource:
    """
    Simple container grouping datasets under a common source label.

    Parameters
    ----------
    name : str
        Source identifier (e.g. ``"GEBCO"``, ``"USGS"``).
    """

    def __init__(self, name: str) -> None:
        self.name = name
        self.dataset = []


def dict2yaml(file_name: str, dct: dict, sort_keys: bool = False) -> None:
    """
    Serialise a dictionary to a YAML file.

    Parameters
    ----------
    file_name : str
        Destination file path.
    dct : dict
        Data to serialise.
    sort_keys : bool, optional
        Whether to sort dictionary keys in the output.  Default is ``False``.
    """
    yaml_string = yaml.dump(dct, sort_keys=sort_keys)
    with open(file_name, "w") as file:
        file.write(yaml_string)


def yaml2dict(file_name: str) -> dict:
    """
    Load a YAML file into a dictionary.

    Parameters
    ----------
    file_name : str
        Path to the YAML file.

    Returns
    -------
    dict
        Parsed contents of the YAML file.
    """
    with open(file_name, "r") as file:
        dct = yaml.load(file, Loader=yaml.FullLoader)
    return dct


def inpolygon(xq: np.ndarray, yq: np.ndarray, p) -> np.ndarray:
    """
    Test which query points fall inside a Shapely polygon.

    Parameters
    ----------
    xq : np.ndarray
        X coordinates of the query points (any shape).
    yq : np.ndarray
        Y coordinates of the query points (same shape as *xq*).
    p : shapely.geometry.Polygon
        Polygon to test against.

    Returns
    -------
    np.ndarray
        Boolean array with the same shape as *xq*; ``True`` where the
        corresponding point lies inside *p*.
    """
    shape = xq.shape
    xq = xq.reshape(-1)
    yq = yq.reshape(-1)
    q = [(xq[i], yq[i]) for i in range(xq.shape[0])]
    p = path.Path([(crds[0], crds[1]) for i, crds in enumerate(p.exterior.coords)])
    return p.contains_points(q).reshape(shape)
