

from celest.coordinates.base_positions import Position3d
from celest.file_save import _save_data_as_txt
from celest.units.quantity import Quantity
from celest import units as u


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

    def __init__(self, julian, latitude, longitude, height, angular_unit,
                 length_unit):

        if julian.ndim != 1 or latitude.ndim != 1 or longitude.ndim != 1 or height.ndim != 1:
            raise ValueError("Input arrays should be one dimensional.")
        if julian.size != latitude.size != longitude.size != height.size:
            raise ValueError("Input arrays should have the same length.")

        x_quantity = Quantity(latitude, angular_unit)
        y_quantity = Quantity(longitude, angular_unit)
        z_quantity = Quantity(height, length_unit)
        julian_quantity = Quantity(julian, u.jd2000)

        super().__init__(x_quantity, y_quantity, z_quantity, julian_quantity)

    @property
    def latitude(self):
        return self.x

    @property
    def longitude(self):
        return self.y

    @property
    def height(self):
        return self.z

    def save_text_file(self, file_name):
        header = "WGS84 Coordinate Data"
        data = [
            ["Time", self.time],
            ["Latitude", self.x],
            ["Longitude", self.y],
            ["Height", self.z]
        ]
        _save_data_as_txt(file_name, header, data=data)
