

from celest.coordinates.frames.base_positions import Position3d
from celest.file_save import _save_data_as_txt
from celest.units.core import Unit
from celest.units.quantity import Quantity
from celest import units as u
import numpy as np


class WGS84(Position3d):
    """Coordinates in the WGS84 Earth ellipsoid model.

    Parameters
    ----------
    julian : np.ndarray
        1-D array containing time in the J2000 epoch.
    latitude, longitude, height : np.ndarray
        1-D array containing the coordinate data.
    angular_unit : Unit
        The unit of the latitude and longitude data.
    length_unit : Unit
        The unit of the height data.

    Attributes
    ----------
    latitude, longitude, height : Quantity
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
    ITRS : International Terrestrial Reference System.
    LVLH : Local vertical local horizontal coordinates.

    Examples
    --------
    Let `julian`, `latitude`, `longitude`, and `height` be `np.ndarray`
    instances. A `WGS84` coordinate object can be initialized:

    >>> wgs84 = WGS84(julian, latitude, longitude, height, u.deg, u.km)

    Save `WGS84` data to a text file:

    >>> wgs84.save_text_file("wgs84_data")
    """

    def __init__(self, julian: np.ndarray, latitude: np.ndarray, longitude:
                 np.ndarray, height: np.ndarray, angular_unit: Unit,
                 length_unit: Unit) -> None:

        super().__init__(latitude, angular_unit, longitude, angular_unit,
                         height, length_unit, julian, u.jd2000)

    @property
    def latitude(self) -> Quantity:
        return self.x

    @property
    def longitude(self) -> Quantity:
        return self.y

    @property
    def height(self) -> Quantity:
        return self.z

    def save_text_file(self, file_name: str) -> None:
        header = "WGS84 Coordinate Data"
        data = [
            ["Time", self.time],
            ["Latitude", self.x],
            ["Longitude", self.y],
            ["Height", self.z]
        ]
        _save_data_as_txt(file_name, header, data=data)
