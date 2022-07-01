

from celest.coordinates.frames.base_positions import Position3d
from celest.coordinates.ground_location import GroundLocation
from celest.file_save import _save_data_as_txt
from celest.units.core import Unit
from celest.units.quantity import Quantity
from celest import units as u
import numpy as np


class LVLH(Position3d):
    """Coordinates in the level-horizontal-level-vertical or Hill frame.

    Parameters
    ----------
    julian : np.ndarray
        1-D array containing time in the J2000 epoch.
    roll, pitch, yaw : np.ndarray
        1-D array containing the coordinate data.
    unit : Unit
        The unit of the spatial data.
    location : GroundLocation
        Target for satellite orientations.

    Attributes
    ----------
    roll, pitch, yaw : Quantity
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
    AzEl : Azimuth elevation coordinates.
    GCRS : Geocentric Celestial Reference System.
    ITRS : International Terrestrial Reference System.
    LVLH : Local vertical local horizontal coordinates.
    WGS84 : Geographical coordinates.

    Examples
    --------
    Let `julian`, `roll`, `pitch`, and `yaw` be `np.ndarray` instances and
    `location` be a `GroundLocation` instance. A `LVLH`coordinate object can be
    initialized:

    >>> lvlh = LVLH(julian, roll, pitch, yaw, u.deg, location)

    Save `LVLH` data to a text file:

    >>> lvlh.save_text_file("lvlh_data")
    """

    def __init__(self, julian: np.ndarray, roll: np.ndarray, pitch: np.ndarray,
                 yaw: np.ndarray, unit: Unit, location: GroundLocation) -> None:

        super().__init__(roll, unit, pitch, unit, yaw, unit, julian, u.jd2000)

        self._location = location

    @property
    def roll(self) -> Quantity:
        return self.x

    @property
    def pitch(self) -> Quantity:
        return self.y

    @property
    def yaw(self) -> Quantity:
        return self.z

    @property
    def location(self) -> GroundLocation:
        return self._location

    def save_text_file(self, file_name: str) -> None:
        header = "Satellite Attitude Data"
        parameters = [
            ["Location", str(self._location)]
        ]
        data = [
            ["Time", self.time],
            ["Roll", self.x],
            ["Pitch", self.y],
            ["Yaw", self.z]
        ]
        _save_data_as_txt(file_name, header, parameters, data)
