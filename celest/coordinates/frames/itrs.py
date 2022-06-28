

from celest.coordinates.base_positions import Position3d
from celest.file_save import _save_data_as_txt
from celest.units.quantity import Quantity
from celest import units as u
import numpy as np


# TODO: Move preprocessing into external functions to remove repeated code.
class ITRS(Position3d):
    """Coordinates in the Internation Terrestrial Reference System.

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
    GCRS : Geocentric Celestial Reference System.
    LVLH : Local vertical local horizontal coordinates.
    WGS84 : Geographical coordinates.

    Examples
    --------
    Let `julian`, `x`, `y`, and `z` be `np.ndarray` instances. An `ITRS`
    coordinate object can be initialized:

    >>> itrs = ITRS(julian, x, y, z, u.km)

    Convert to the `GCRS` frame:

    >>> gcrs = itrs.convert_to(GCRS)

    Save `ITRS` data to a text file:

    >>> itrs.save_text_file("gcrs_data")
    """

    def __init__(self, julian, x, y, z, unit):

        if julian.ndim != 1 or x.ndim != 1 or y.ndim != 1 or z.ndim != 1:
            raise ValueError("Input arrays should be one dimensional.")
        if julian.size != x.size != y.size != z.size:
            raise ValueError("Input arrays should have the same length.")

        x_quantity = Quantity(x, unit)
        y_quantity = Quantity(y, unit)
        z_quantity = Quantity(z, unit)
        julian_quantity = Quantity(julian, u.jd2000)

        super().__init__(x_quantity, y_quantity, z_quantity, julian_quantity)

    def save_text_file(self, file_name):
        header = "ITRS Coordinate Data"
        data = [
            ["Time", self.time],
            ["X", self.x],
            ["Y", self.y],
            ["Z", self.z]
        ]
        _save_data_as_txt(file_name, header, data=data)
