

from celest.coordinates.frames.base_positions import Position3d
from celest.file_save import TextFileWriter
from celest.units.core import Unit
from celest.units.quantity import Quantity
from celest import units as u
import numpy as np


class WGS84(Position3d):
    """WGS84(julian, latitude, longitude, height, angular_unit, length_unit)

    Coordinates in the WGS84 Earth ellipsoid model.

    The World Geodetic System 84 (WGS84) is an Earth ellipsoid model with its
    origin located at the Earth's center of mass. The WGS84 meridian of zero
    longitude is located at the IERS reference meridian with the parallel of
    zero latitude located at the WGS84 reference meridian plane. WGS84 is the
    standard model used by the Global Positioning System (GPS).

    Parameters
    ----------
    julian : np.ndarray
        1-D array containing time in the JD2000 epoch.
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
    save_text_file(file_name)
        Save data as a pretty text file.

    See Also
    --------
    Attitude : Satellite attitude.
    AzEl : Azimuth-elevation coordinates.
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
        """Coordinates in the WGS84 Earth ellipsoid model.
        
        The World Geodetic System 84 (WGS84) is an Earth ellipsoid model with
        its origin located at the Earth's center of mass. The WGS84 meridian of
        zero longitude is located at the IERS reference meridian with the
        parallel of zero latitude located at the WGS84 reference meridian plane.
        WGS84 is the standard model used by the Global Positioning System (GPS).

        Parameters
        ----------
        julian : np.ndarray
            1-D array containing time in the JD2000 epoch.
        latitude, longitude, height : np.ndarray
            1-D array containing the coordinate data.
        angular_unit : Unit
            The unit of the latitude and longitude data.
        length_unit : Unit
            The unit of the height data.
        """

        super().__init__(latitude, angular_unit, longitude, angular_unit,
                         height, length_unit, julian, u.jd2000)

    @property
    def latitude(self) -> Quantity:
        return self._get_x()

    @property
    def longitude(self) -> Quantity:
        return self._get_y()

    @property
    def height(self) -> Quantity:
        return self._get_z()

    def save_text_file(self, file_name: str) -> None:
        """Save data as a pretty text file.

        Parameters
        ----------
        file_name : str
            Name of the text file for the saved data.
        """

        header = "WGS84 Coordinate Data"
        data = [
            ["Time", self.time],
            ["Latitude", self.latitude],
            ["Longitude", self.longitude],
            ["Height", self.height]
        ]
        writer = TextFileWriter(file_name, header)
        writer.add_layer(data=data)
        writer.save()
