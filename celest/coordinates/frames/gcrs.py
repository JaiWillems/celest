

from celest.coordinates.frames.base_positions import Position3d
from celest.file_save import _save_data_as_txt
from celest.units.core import Unit
from celest.units.quantity import Quantity
from celest import units as u
import numpy as np


class GCRS(Position3d):
    """Coordinates in the Geocentric Celestial Reference System.

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
    time : Quantity
        Times associated with coordinate dimensions.

    Methods
    -------
    convert_to(frame, **kwargs)
        Convert current frame into a new reference frame.
    save_text_file(file_name)
        Save data as a pretty text file.

    See Also
    --------
    AzEl : Azimuth elevation coordinates.
    ITRS : International Terrestrial Reference System.
    LVLH : Local vertical local horizontal coordinates.
    WGS84 : Geographical coordinates.

    Examples
    --------
    Let `julian`, `x`, `y`, and `z` be `np.ndarray` instances. A `GCRS`
    coordinate object can be initialized:

    >>> gcrs = GCRS(julian, x, y, z, u.km)

    Convert to the `ITRS` frame:

    >>> itrs = gcrs.convert_to(ITRS)

    Save `GCRS` data to a text file:

    >>> gcrs.save_text_file("gcrs_data")
    """

    def __init__(self, julian: np.ndarray, x: np.ndarray, y: np.ndarray, z:
                 np.ndarray, unit: Unit) -> None:

        super().__init__(x, unit, y, unit, z, unit, julian, u.jd2000)

    def save_text_file(self, file_name: str) -> None:
        header = "GCRS Coordinate Data"
        data = [
            ["Time", self.time],
            ["X", self.x],
            ["Y", self.y],
            ["Z", self.z]
        ]
        _save_data_as_txt(file_name, header, data=data)
