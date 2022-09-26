

from celest.coordinates.frames.base_positions import Position3d
from celest.file_save import TextFileWriter
from celest.units.core import Unit
from celest.units.quantity import Quantity
from celest import units as u
import numpy as np


class GCRS(Position3d):
    """GCRS(julian, x, y, z, unit)

    Coordinates in the Geocentric Celestial Reference System.

    Parameters
    ----------
    julian : np.ndarray
        1-D array containing time in the JD2000 epoch.
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
    save_text_file(file_name)
        Save data as a pretty text file.

    See Also
    --------
    Attitude : Satellite attitude.
    AzEl : Azimuth-elevation coordinates.
    ITRS : International Terrestrial Reference System.
    LVLH : Local vertical local horizontal coordinates.
    WGS84 : Geographical coordinates.

    Examples
    --------
    Let `julian`, `x`, `y`, and `z` be `np.ndarray` instances. A `GCRS`
    coordinate object can be initialized:

    >>> gcrs = GCRS(julian, x, y, z, u.km)

    Save `GCRS` data to a text file:

    >>> gcrs.save_text_file("gcrs_data")
    """

    def __init__(self, julian: np.ndarray, x: np.ndarray, y: np.ndarray, z:
                 np.ndarray, unit: Unit) -> None:
        """Coordinates in the Geocentric Celestial Reference System.

        Parameters
        ----------
        julian : np.ndarray
            1-D array containing time in the JD2000 epoch.
        x, y, z : np.ndarray
            1-D array containing the coordinate data.
        unit : Unit
            The unit of the spatial data.
        """

        super().__init__(x, unit, y, unit, z, unit, julian, u.jd2000)

    @property
    def x(self) -> Quantity:
        return self._get_x()

    @property
    def y(self) -> Quantity:
        return self._get_y()

    @property
    def z(self) -> Quantity:
        return self._get_z()

    def save_text_file(self, file_name: str) -> None:
        """Save data as a pretty text file.

        Parameters
        ----------
        file_name : str
            Name of the text file for the saved data.
        """

        header = "GCRS Coordinate Data"
        data = [
            ["Time", self.time],
            ["X", self.x],
            ["Y", self.y],
            ["Z", self.z]
        ]
        writer = TextFileWriter(file_name, header)
        writer.add_layer(data=data)
        writer.save()
