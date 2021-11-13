"""Localize ground position information.

This module contains the `GoundPosition` class which stores ground location
data and is an instrumental part of encounter planning by managing the ground
location and satellite encounter relationship.
"""


from celest.core.decorators import set_module
import numpy as np


@set_module('celest.encounter')
class GroundPosition(object):
    """Localize Earth bound location information.

    The `GoundPosition` class stores ground location information and is an
    instrumental part of encounter planning functionality.

    Parameters
    ----------
    latitude : float
        Locations geodetic lattitude in degrees and decimals.
    longitude : float
        Locations geodetic longitude in degrees and decimals.

    Attributes
    ----------
    lat : float
        Geodetic latitude for the location given in degrees and decimals.
    lon : float
        Geodetic longitude for the location given in degrees and decimals.
    radius : float
        Earth's radius at the given coordinates.
    """

    def __init__(self, latitude: float, longitude: float) -> None:
        """Initialize attributes."""

        self.lat = latitude
        self.lon = longitude
        self.radius = self._radius(latitude)

    def __str__(self) -> str:
        """Define information string."""

        output = [
            f"coor={(self.lat, self.lon)}",
            f"radius={round(self.radius, 5)}"
        ]

        return ", ".join(output)

    def _radius(self, latitude: float) -> float:
        """Calculate the Earth's geocentric radius using WGS84.

        Parameters
        ----------
        latitude : float
            Latitude of the ground location in degrees and decimals.

        Returns
        -------
        float
            Earth's geocentric radius in kilometres and decimals.

        Notes
        -----
        By using an Earth ellipsoid with the WGS84 parameters of
        :math:`a=6378.137` and :math:`b=6356.7523142`, the geocentric radius
        can be calculated using the following formulation:

        .. math:: r = \sqrt{\frac{(a^2\cos(\beta))^2 + (b^2\sin(\beta))^2}{(a\cos(\beta))^2 + (b\sin(\beta))^2}}

        where :math:`\beta` is the observer's latitude.[1]_

        References
        ----------
        .. [1] Timur. Earth Radius by Latitude (WGS 84). 2018. url:
           https://planetcalc.com/7721/.
        """

        # Get lattidue parameter.
        phi = np.radians(latitude)

        # Define WGS84 Parameters.
        semi_major = 6378.137
        semi_minor = 6356.752314245

        c_phi, s_phi = np.cos(phi), np.sin(phi)

        num = (semi_major ** 2 * c_phi) ** 2 + (semi_minor ** 2 * s_phi) ** 2
        denom = (semi_major * c_phi) ** 2 + (semi_minor * s_phi) ** 2
        radius = np.sqrt(num / denom)

        return radius
