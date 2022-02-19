

import numpy as np


class GroundPosition(object):
    """GroundPosition(latitude, longitude)

    Localize Earth surface location information.

    Parameters
    ----------
    latitude : float
        Lattitude of the location in decimal degrees.
    longitude : float
        Longitude of the location in decimal degrees.

    Attributes
    ----------
    lat, lon : float
        Latitude and longitude of the location in decimal degrees.
    radius : float
        Earth radius at (`latitude`, `longitude`).
    """

    def __init__(self, latitude: float, longitude: float) -> None:

        self.lat = latitude
        self.lon = longitude
        self.radius = self._radius(latitude)

    def __str__(self) -> str:

        return f"coor={(self.lat, self.lon)}, radius={round(self.radius, 5)}"

    def _radius(self, latitude: float) -> float:
        """Calculate geocentric radius using WGS84.

        Parameters
        ----------
        latitude : float
            Location's latitude in decimal degrees.

        Returns
        -------
        float
            Earth's geocentric radius in decimal kilometers.

        Notes
        -----
        The WGS84 Earth ellipsoid model is used as discussed in "Earth Radius
        by Latitude (WGS 84)" by Timur. [1]_

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
