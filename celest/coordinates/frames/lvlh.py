

from celest.coordinates.frames.base_positions import Position3d
from celest.coordinates.ground_location import GroundLocation
from celest.file_save import _save_data_as_txt
from celest.units.core import Unit
from celest.units.quantity import Quantity
from celest import units as u
import numpy as np


class LVLH(Position3d):
    """Coordinates in the level-horizontal-level-vertical or Hill frame.

    The local-vertical local-horizontal frame (also known as the Hill
    frame) is a body frame where the z-axis is algigned with the negative
    of the geocentric position vector, the y-axis is aligned with the
    negative orbit normal, and the x-axis completes the right handed triad.

    Parameters
    ----------
    julian : np.ndarray
        1-D array containing time in the J2000 epoch.
    x, y, z : np.ndarray
        1-D array containing the coordinate data.
    unit : Unit
        The unit of the spatial data.

    Attributes
    ----------
    x, y, z : Quantity
        Coordinate data.
    location : GroundLocation
        Origin of the AzEl frame.
    time : Quantity
        Times associated with coordinate dimensions.

    Methods
    -------
    save_text_file(file_name)
        Save data as a pretty text file.

    See Also
    --------
    Attitude : Satellite attitude.
    AzEl : Azimuth elevation coordinates.
    GCRS : Geocentric Celestial Reference System.
    ITRS : International Terrestrial Reference System.
    WGS84 : Geographical coordinates.

    Notes
    -----
    The LVLH frame definition was taken from NASA's technical memorandum
    on coordinate frames for the space shuttle program. [NASA1974]_

    References
    ----------
    .. [NASA1974] Coordinate Systems for the Space Shuttle Program, Lyndon
       B. Johnson Space Center, Houston, Texas 77058, Oct 1974, no. NASA
       TM X-58153.

    Examples
    --------
    Let `julian`, `x`, `y`, and `z` be `np.ndarray` instances. A `LVLH`
    coordinate object can be initialized:

    >>> lvlh = LVLH(julian, x, y, z, u.km)

    Save `LVLH` data to a text file:

    >>> lvlh.save_text_file("lvlh_data")
    """

    def __init__(self, julian: np.ndarray, x: np.ndarray, y: np.ndarray,
                 z: np.ndarray, unit: Unit) -> None:

        super().__init__(x, unit, y, unit, z, unit, julian, u.jd2000)

    def save_text_file(self, file_name: str) -> None:
        header = "Satellite LVLH Coordinate Data"
        data = [
            ["Time", self.time],
            ["X", self.x],
            ["Y", self.y],
            ["Z", self.z]
        ]
        _save_data_as_txt(file_name, header, data=data)
