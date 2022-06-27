

from celest.units.quantity import Quantity
from celest.units.core import Unit
from celest import units as u
from celest.constants import WGS84_MINOR_AXIS_KM, WGS84_MAJOR_AXIS_KM
from celest.angle_strings import _ISO6709_representation
from math import cos, sin, sqrt


class GroundLocation:
    """GroundLocation(latitude, longitude, height, angular_unit, length_unit)

    Specify Earth bound location.

    Parameters
    ----------
    latitude : float
    longitude : float
    height : float, optional
        Height of the location above the standard ellipsoid.
    angular_unit : Unit
        Unit of the `latitude` and `longitude` measures.
    length_unit : Unit
        Unit of the height measure.

    Attributes
    ----------
    latitude : Quantity
    longitude : Quantity
    height : Quantity
    radius : Quantity
        Earth radius at ground level at the specified position.

    """

    def __init__(self, latitude: float, longitude: float, height: float,
                 angular_unit: Unit, length_unit: Unit) -> None:

        self._latitude = Quantity(latitude, angular_unit)
        self._longitude = Quantity(longitude, angular_unit)
        self._height = Quantity(height, length_unit)

    def __str__(self) -> str:
        return _ISO6709_representation(self._latitude, self.longitude,
                                       self._height)

    def __repr__(self) -> str:
        return "GroundLocation(" + str(self._latitude.data) + ", " +\
            str(self._longitude.data) + ", " + str(self._height.data) +\
            ", " + repr(self._latitude.unit) + ", " + repr(self._height.unit)

    @property
    def latitude(self) -> Quantity:
        return self._latitude

    @property
    def longitude(self) -> Quantity:
        return self._longitude

    @property
    def height(self) -> Quantity:
        return self._height

    @property
    def radius(self) -> Quantity:
        """Calculate geocentric radius using WGS84.

        Returns
        -------
        Quantity
            Earth's geocentric radius at the given location.

        Notes
        -----
        The WGS84 Earth ellipsoid model is used as discussed in "Earth Radius
        by Latitude (WGS 84)" by Timur. [1]_

        References
        ----------
        .. [1] Timur. Earth Radius by Latitude (WGS 84). 2018. url:
           https://planetcalc.com/7721/.
        """
        cos_latitude = cos(self._latitude.to(u.rad).data)
        sin_latitude = sin(self._latitude.to(u.rad).data)

        numerator = (WGS84_MAJOR_AXIS_KM ** 2 * cos_latitude) ** 2 + \
                    (WGS84_MINOR_AXIS_KM ** 2 * sin_latitude) ** 2
        denominator = (WGS84_MAJOR_AXIS_KM * cos_latitude) ** 2 + \
                      (WGS84_MINOR_AXIS_KM * sin_latitude) ** 2
        radius_km = sqrt(numerator / denominator) + self._height.to(u.km).data

        return Quantity(radius_km, u.km)

    @property
    def itrs_x(self):
        m = self._meridional_radius_of_curvature()
        itrs_x = (m + self._height.to(u.km).data) * \
            cos(self._latitude.to(u.rad).data) * \
            cos(self._longitude.to(u.rad).data)

        return Quantity(itrs_x, u.km)

    def _ellipse_eccentricity(self):
        return sqrt(1 - WGS84_MINOR_AXIS_KM ** 2 / WGS84_MAJOR_AXIS_KM ** 2)

    def _meridional_radius_of_curvature(self):
        return WGS84_MAJOR_AXIS_KM / sqrt(1 - self._ellipse_eccentricity() **
                                          2 * sin(self._latitude.to(u.rad).data) ** 2)

    @property
    def itrs_y(self):
        m = self._meridional_radius_of_curvature()
        itrs_y = (m + self._height.to(u.km).data) * \
            cos(self._latitude.to(u.rad).data) * \
            sin(self._longitude.to(u.rad).data)

        return Quantity(itrs_y, u.km)

    @property
    def itrs_z(self):
        e = self._ellipse_eccentricity()
        m = self._meridional_radius_of_curvature()
        itrs_z = (m * (1 - e ** 2) + self._height.to(u.km).data) * \
            sin(self._latitude.to(u.rad).data)

        return Quantity(itrs_z, u.km)
