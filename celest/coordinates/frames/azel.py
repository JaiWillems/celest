

from celest.coordinates.frames.base_positions import Position2d
from celest.coordinates.ground_location import GroundLocation
from celest.file_save import _save_data_as_txt
from celest.units.core import Unit
from celest.units.quantity import Quantity
from celest import units as u
import numpy as np


class AzEl(Position2d):
    """Coordinates in the horizontal system.

    Parameters
    ----------
    julian : np.ndarray
        1-D array containing time in the J2000 epoch.
    azimuth, elevation : np.ndarray
        1-D array containing the coordinate data.
    unit : Unit
        The unit of the spatial data.
    location : GroundLocation
        Origin of the azel coordinate frame.

    Attributes
    ----------
    azimuth, elevation : Quantity
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
    GCRS : Geocentric Celestial Reference System.
    ITRS : International Terrestrial Reference System.
    LVLH : Local vertical local horizontal coordinates.
    WGS84 : Geographical coordinates.

    Examples
    --------
    Let `julian`, `azimuth`, and `elevation` be `np.ndarray` instances and
    `location` be a `GroundLocation` instance. An `AzEl`coordinate object can
    be initialized:

    >>> azel = AzEl(julian, azimuth, elevation, u.deg, location)

    Save `AzEl` data to a text file:

    >>> azel.save_text_file("azel_data")
    """

    def __init__(self, julian: np.ndarray, azimuth: np.ndarray, elevation:
                 np.ndarray, unit: Unit, location: GroundLocation) -> None:

        super().__init__(azimuth, unit, elevation, unit, julian, u.jd2000)

        self._location = location

    @property
    def azimuth(self) -> Quantity:
        return self.x

    @property
    def elevation(self) -> Quantity:
        return self.y

    @property
    def location(self) -> GroundLocation:
        return self._location

    def save_text_file(self, file_name: str) -> None:
        header = "Azimuth-Elevation Coordinate Data"
        parameters = [
            ["Location", str(self._location)]
        ]
        data = [
            ["Time", self.time],
            ["Azimuth", self.x],
            ["Elevation", self.y]
        ]
        _save_data_as_txt(file_name, header, parameters, data)
