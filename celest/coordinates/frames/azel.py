

from celest.coordinates.frames.base_positions import Position2d
from celest.coordinates.ground_location import GroundLocation
from celest.file_save import TextFileWriter
from celest.units.core import Unit
from celest.units.quantity import Quantity
from celest import units as u
import numpy as np


class AzEl(Position2d):
    """AzEl(julian, azimuth, elevation, unit, location)

    Coordinates in the horizontal system.

    The horizontal system has its origin located at a ground location and has
    measures of azimith and elevation. Azimuth is defined as the angle of the
    satellite in the horizontal plane clockwise from north. Elevation is defined
    as the angle of the satellite above the horzontal plane.

    The horizontal frame is sometimes refered to as the azimuth-elevation
    (az-el) or altitude-azimuth (alt-az) frame.

    Parameters
    ----------
    julian : np.ndarray
        1-D array containing time in the JD2000 epoch.
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
    Attitude : Satellite attitude.
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
        """Coordinates in the horizontal system.

        The horizontal system has its origin located at a ground location and
        has measures of azimith and elevation. Azimuth is defined as the angle
        of the satellite in the horizontal plane clockwise from north. Elevation
        is defined as the angle of the satellite above the horzontal plane.
    
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
        """

        super().__init__(azimuth, unit, elevation, unit, julian, u.jd2000)

        self._location = location

    @property
    def azimuth(self) -> Quantity:
        return self._get_x()

    @property
    def elevation(self) -> Quantity:
        return self._get_y()

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

        header = "Azimuth-Elevation Coordinate Data"
        parameters = [
            ["Location", str(self._location)]
        ]
        data = [
            ["Time", self.time],
            ["Azimuth", self.azimuth],
            ["Elevation", self.elevation]
        ]
        writer = TextFileWriter(file_name, header)
        writer.add_layer(parameters=parameters, data=data)
        writer.save()
