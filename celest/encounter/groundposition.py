

import numpy as np


WGS84_MAJOR_AXIS = 6378.137
WGS84_MINOR_AXIS = 6356.752314245


class GroundPosition:
    """GroundPosition(latitude, longitude, height=0)

    Localize Earth surface location information.

    Parameters
    ----------
    latitude : float
        Lattitude of the location in decimal degrees.
    longitude : float
        Longitude of the location in decimal degrees.
    height : float, optional
        Height of the location above the standard ellipsoid in decimal
        kilometers.

    Attributes
    ----------
    latitude, longitude : float
        Latitude and longitude of the location in decimal degrees.
    radius : float
        Earth radius at (`latitude`, `longitude`).
    """

    def __init__(self, latitude: float, longitude: float, height: float=0) -> None:

        self.latitude = latitude
        self.longitude = longitude
        self.height = height
        self.radius = self._radius(latitude) + height

    def __str__(self) -> str:

        return f"coor={(self.latitude, self.longitude)}, radius={round(self.radius, 5)}"

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

        cos_latitude = np.cos(np.radians(latitude))
        sin_latitude = np.sin(np.radians(latitude))

        numerator = (WGS84_MAJOR_AXIS ** 2 * cos_latitude) ** 2 + \
            (WGS84_MINOR_AXIS ** 2 * sin_latitude) ** 2
        denominator = (WGS84_MAJOR_AXIS * cos_latitude) ** 2 + \
            (WGS84_MINOR_AXIS * sin_latitude) ** 2

        return np.sqrt(numerator / denominator)
