

from celest.angle_strings import _ISO6709_representation
from celest.coordinates.astronomical_quantities import radius_of_curvature_in_prime_vertical
from celest.constants import (
    WGS84_MINOR_AXIS_KM,
    WGS84_MAJOR_AXIS_KM,
    WGS83_FIRST_ECCENTRICITY
)
from celest.units.quantity import Quantity
from celest.units.core import Unit
from celest import units as u
from math import cos, sin, sqrt


class GroundLocation:
    """GroundLocation(latitude, longitude, height, angular_unit, length_unit)

    Specify an Earth bound location.

    The `GroundLocation` class models the Earth with the WGS84 reference
    ellipsoid to provide accurate representations of ground locations.

    Parameters
    ----------
    latitude : float
        The location's latitude.
    longitude : float
        The location's longitude.
    height : float, optional
        The location's height above the WGS84 ellipsoid.
    angular_unit : Unit
        Unit of the `latitude` and `longitude` measures.
    length_unit : Unit
        Unit of the `height` measure.

    Attributes
    ----------
    latitude : Quantity
    longitude : Quantity
    height : Quantity
    radius : Quantity
        Earth radius at ground level at the specified position.
    itrs_x : Quantity
        Location's X coordinate in the ITRS frame.
    itrs_y : Quantity
        Location's Y coordinate in the ITRS frame.
    itrs_z : Quantity
        Location's Z coordinate in the ITRS frame.
    """

    def __init__(self, latitude: float, longitude: float, height: float,
                 angular_unit: Unit, length_unit: Unit) -> None:
        """Specify an Earth bound location.

        The `GroundLocation` class models the Earth with the WGS84 reference
        ellipsoid to provide accurate representations of ground locations.

        Parameters
        ----------
        latitude : float
            The location's latitude.
        longitude : float
            The location's longitude.
        height : float, optional
            The location's height above the WGS84 ellipsoid.
        angular_unit : Unit
            Unit of the `latitude` and `longitude` measures.
        length_unit : Unit
            Unit of the `height` measure.
        """

        if latitude < -90 or latitude > 90:
            raise ValueError("Latitude must be between -90 and 90 degrees.")
        if longitude < -180 or longitude > 180:
            raise ValueError("Longitude must be between -180 and 180 degrees.")

        self._latitude = Quantity(latitude, angular_unit)
        self._longitude = Quantity(longitude, angular_unit)
        self._height = Quantity(height, length_unit)

    def __str__(self) -> str:
        return _ISO6709_representation(
            self._latitude,
            self._longitude,
            self._height
        )

    def __repr__(self) -> str:
        return "GroundLocation(" + ", ".join([
            str(self._latitude.to(self._latitude.unit)),
            str(self._longitude.to(self._longitude.unit)),
            str(self._height.to(self._height.unit)),
            repr(self._latitude.unit),
            repr(self._height.unit)
        ]) + ")"

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
        """Return the location's geocentric radius using WGS84.

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

        cos_latitude = cos(self._latitude.to(u.rad))
        sin_latitude = sin(self._latitude.to(u.rad))

        numerator = (WGS84_MAJOR_AXIS_KM ** 2 * cos_latitude) ** 2 + \
                    (WGS84_MINOR_AXIS_KM ** 2 * sin_latitude) ** 2
        denominator = (WGS84_MAJOR_AXIS_KM * cos_latitude) ** 2 + \
                      (WGS84_MINOR_AXIS_KM * sin_latitude) ** 2
        radius_km = sqrt(numerator / denominator) + self._height.to(u.km)

        return Quantity(radius_km, u.km)

    @property
    def itrs_x(self) -> Quantity:
        """
        Location's x coordinate in the itrs frame.
        """

        m = radius_of_curvature_in_prime_vertical(self._latitude)
        return (m + self._height) * cos(self._latitude.to(u.rad)) * \
            cos(self._longitude.to(u.rad))

    @property
    def itrs_y(self) -> Quantity:
        """
        Location's y coordinate in the itrs frame.
        """

        m = radius_of_curvature_in_prime_vertical(self.latitude)
        return (m + self._height) * cos(self._latitude.to(u.rad)) * \
            sin(self._longitude.to(u.rad))

    @property
    def itrs_z(self) -> Quantity:
        """
        Location's z coordinate in the itrs frame.
        """

        m = radius_of_curvature_in_prime_vertical(self.latitude)
        return (m * (1 - WGS83_FIRST_ECCENTRICITY ** 2) + self._height) * \
            sin(self._latitude.to(u.rad))
