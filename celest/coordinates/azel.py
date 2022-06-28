

from celest.coordinates.base_positions import Position2d
from celest.file_save import _save_data_as_txt
from celest.units.quantity import Quantity
from celest import units as u


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
    `location` be a `GroundPosition` instance. An `AzEl`coordinate object
    can be initialized:

    >>> azel = AzEl(julian, azimuth, elevation, u.deg, location)

    Save `AzEl` data to a text file:

    >>> azel.save_text_file("azel_data")
    """

    def __init__(self, julian, azimuth, elevation, unit, location):

        if julian.ndim != 1 or azimuth.ndim != 1 or elevation.ndim != 1:
            raise ValueError("Input arrays should be one dimensional.")
        if julian.size != azimuth.size != elevation.size:
            raise ValueError("Input arrays should have the same length.")

        x_quantity = Quantity(azimuth, unit)
        y_quantity = Quantity(elevation, unit)
        julian_quantity = Quantity(julian, u.jd2000)

        super().__init__(x_quantity, y_quantity, julian_quantity)

        self._location = location

    @property
    def azimuth(self):
        return self.x

    @property
    def elevation(self):
        return self.y

    @property
    def location(self):
        return self._location

    def save_text_file(self, file_name):
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
