"""
Tiled web map bathymetry dataset.

Wraps a :class:`cht_tiling.TiledWebMap` instance as a bathymetry dataset so
that tile-based elevation services (e.g. Mapbox Terrain-RGB) can be queried
through the standard :class:`~cht_bathymetry.dataset.BathymetryDataset`
interface.  Data are always returned in EPSG:3857.
"""

from cht_tiling import TiledWebMap
from pyproj import CRS

from .dataset import BathymetryDataset


class BathymetryDatasetTiledWebMap(BathymetryDataset):
    """
    Bathymetry dataset backed by a tiled web map (XYZ/TMS tile service).

    Coordinates are in EPSG:3857 (Web Mercator).
    """

    def __init__(self, name: str, path: str) -> None:
        """
        Initialise a tiled web map dataset.

        Parameters
        ----------
        name : str
            Short identifier for the dataset; used for metadata lookup and
            tile file naming.
        path : str
            Local directory containing the dataset metadata and cached tiles.
        """
        super().__init__()

        self.name = name
        self.path = path
        self.local_path = path
        self.read_metadata()
        self.data = TiledWebMap(self.local_path, name, parameter="elevation")
        self.crs = CRS(3857)

    def get_data(
        self,
        xl: list[float],
        yl: list[float],
        max_cell_size: float = 1000.0,
        waitbox: None = None,
    ) -> tuple:
        """
        Read elevation data from the tiled web map service.

        Delegates entirely to :meth:`cht_tiling.TiledWebMap.get_data`.
        Coordinates and depth values are returned in EPSG:3857.

        Parameters
        ----------
        xl : list[float]
            ``[x_min, x_max]`` bounds in EPSG:3857 (metres).
        yl : list[float]
            ``[y_min, y_max]`` bounds in EPSG:3857 (metres).
        max_cell_size : float, optional
            Maximum pixel size in metres used to select the zoom level.
            Default is ``1000.0``.
        waitbox : None, optional
            Reserved for a progress-dialog object; currently unused.

        Returns
        -------
        x : np.ndarray
            1-D array of x coordinates (EPSG:3857, metres).
        y : np.ndarray
            1-D array of y coordinates (EPSG:3857, metres).
        z : np.ndarray
            2-D elevation array (NaN where no data).
        """
        x, y, z = self.data.get_data(xl, yl, max_cell_size, waitbox=waitbox)

        return x, y, z
