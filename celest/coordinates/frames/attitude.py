

from celest.coordinates.frames.base_positions import Position3d
from celest.coordinates.ground_location import GroundLocation
from celest.file_save import TextFileWriter
from celest.units.core import Unit
from celest.units.quantity import Quantity
from celest import units as u
import numpy as np


class Attitude(Position3d):
    """Attitude(julian, roll, pitch, yaw, unit, location)

    Satellite attitude.

    The attitude of a satellite is defined by the roll, pitch, and yaw angles
    that take the satellite from the lvlh frame to a ground-target-pointing
    orientation.

    Parameters
    ----------
    julian : np.ndarray
        1-D array containing time in the JD2000 epoch.
    roll, pitch, yaw : np.ndarray
        1-D array containing the coordinate data.
    unit : Unit
        The unit of the angular data.
    location : GroundLocation
        The ground location associated with the attitude data.

    Attributes
    ----------
    roll, pitch, yaw : Quantity
        Coordinate data.
    location : GroundLocation
        The ground location associated with the attitude data.
    time : Quantity
        Times associated with coordinate dimensions.

    Methods
    -------
    save_text_file(file_name)
        Save data as a pretty text file.

    See Also
    --------
    AzEl : Azimuth-elevation coordinates.
    GCRS : Geocentric Celestial Reference System.
    ITRS : International Terrestrial Reference System.
    LVLH : Local vertical local horizontal coordinates.
    WGS84 : Geographical coordinates.

    Examples
    --------
    Let `julian`, `roll`, `pitch`, and `yaw` be `np.ndarray` instances and
    `location` be a `GroundLocation` instance. An `Attitude` object can be
    initialized:

    >>> attitude = Attitude(julian, roll, picth, yaw, u.deg, location)

    Save `Attitude` data to a text file:

    >>> attitude.save_text_file("attitude_data")
    """

    def __init__(self, julian: np.ndarray, roll: np.ndarray, pitch: np.ndarray,
                 yaw: np.ndarray, unit: Unit, location: GroundLocation) -> None:
        """Satellite attitude.

        The attitude of a satellite is defined by the roll, pitch, and yaw
        angles that take the satellite from the lvlh frame to a
        ground-target-pointing orientation.

        Parameters
        ----------
        julian : np.ndarray
            1-D array containing time in the JD2000 epoch.
        roll, pitch, yaw : np.ndarray
            1-D array containing the coordinate data.
        unit : Unit
            The unit of the angular data.
        location : GroundLocation
            The ground location associated with the attitude data.
        """

        super().__init__(roll, unit, pitch, unit, yaw, unit, julian, u.jd2000)

        self._location = location

    @property
    def roll(self) -> Quantity:
        return self._get_x()

    @property
    def pitch(self) -> Quantity:
        return self._get_y()

    @property
    def yaw(self) -> Quantity:
        return self._get_z()

    @property
    def location(self) -> GroundLocation:
        return self._location

    def save_text_file(self, file_name: str) -> None:
        """Save data as a pretty text file.

        Parameters
        ----------
        file_name : str
            Name of the text file for the saved data.
        """

        header = "Attitude Coordinate Data"
        parameters = [
            ["Location", str(self._location)]
        ]
        data = [
            ["Time", self.time],
            ["Roll", self.roll],
            ["Pitch", self.pitch],
            ["Yaw", self.yaw]
        ]
        writer = TextFileWriter(file_name, header)
        writer.add_layer(parameters=parameters, data=data)
        writer.save()
