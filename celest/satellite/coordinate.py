"""
"""


from datetime import time
from celest.core.decorators import set_module
from celest.satellite import Time, Interpolation
from celest.encounter import GroundPosition
import numpy as np


@set_module('celest.satellite')
class Coordinate(object):
    """
    """

    def __init__(self, basePos: np.array, type: str, timeData: Time, factor: int=0) -> None:
        """Initialize attributes."""

        self.time = timeData

        self._GEO = None
        self._ERA = None
        self._ICRS = None
        self._GCRS = None
        self._ITRS = None
        self._CIRS = None
        self._altitude = None
        self._equatorial = None
        self._ecliptic = None
        self._galactic = None
        self._super_galactic = None

        self.interp = Interpolation()

        self.length = None
        self._set_base_position(basePos, type, factor)
    
    def _set_base_position(self, basePos: np.array, type: str, factor: int) -> None:
        pass
    
    def GEO(self, **kwargs) -> np.array:
        pass
    
    def ERA(self, **kwargs) -> np.array:
        pass
    
    def ICRS(self, **kwargs) -> np.array:
        pass
    
    def GCRS(self, **kwargs) -> np.array:
        pass
    
    def ITRS(self, **kwargs) -> np.array:
        pass
    
    def CIRS(self, **kwargs) -> np.array:
        pass
    
    def horizontal(self, goundPos: GroundPosition, **kwargs) -> np.array:
        pass
    
    def off_nadir(self, groundPos: GroundPosition, **kwargs) -> np.array:
        pass

    def _WGS84_radius(self, lattitude: np.array) -> np.array:
        """Calculate the Earth's radius using WGS84.

        Parameters
        ----------
        latitude : np.array
            Array of shape (n,) representing the lattitude of a ground
            location.

        Returns
        -------
        np.array
            Earth's radius at each row in `lattitude` using WGS84.

        Notes
        -----
        The Earth can be modeled as an ellipsoid given by the following
        equation:

        .. math:: r = \sqrt{\\frac{(6378.14)^2(6356.75)^2}{(6378.14)^2\sin{\phi}^2+(6356.75)^2\cos{\phi}^2}}

        where :math:`\phi` is the observers lattitude.
        """

        # Get lattidue parameter.
        phi = np.radians(lattitude)

        # Define WGS84 Parameters.
        semiMajor = 6378.137**2
        semiMinor = 6356.752314245**2

        numerator = semiMajor * semiMinor
        denominator = semiMajor * np.sin(phi)**2 + semiMinor * np.cos(phi)**2
        radius = np.sqrt(numerator / denominator)

        return radius
    
    def altitude(self, **kwargs) -> np.array:
        """Get the altitude above Earth's surface.

        Parameters
        ----------
        kwargs : dict, optional
            Optional keyword arguments passed into the `Interpolation()`
            callable.

        Returns
        -------
        np.array
            Array of shape (n,) containing time varying altitudes.
        
        Notes
        -----
        This method implements WGS84 to calculate the Earth's radius before
        computing the positions altitude.
        """

        if type(self._altitude) != type(None):
            if kwargs:
                return self.interp(self._altitude)
            else:
                return self._altitude

        if type(self._ITRS) == type(None):
            self.ITRS()

        ITRSdata = self._ITRS
        x = ITRSdata[:, 0]
        y = ITRSdata[:, 1]
        z = ITRSdata[:, 2]
        arg = np.sqrt(x**2 + y**2) / z
        lattitude = np.arctan(arg)

        earthRadius = self._WGS84_radius(lattitude)
        satRadius = np.linalg.norm(ITRSdata, axis=1)

        altitude = satRadius - earthRadius
        self._altitude = altitude

        if kwargs:
            altitude = self.interp(altitude, **kwargs)

        return altitude
    
    def distance(self, groundPos: GroundPosition, **kwargs) -> np.array:
        """Get distance to a ground location.

        Parameters
        ----------
        groundPos : GroundPosition
            `GroundPosition` object instantiated as per its documentation.
        kwargs : dict, optional
            Optional keyword arguments passed into the `Interpolation()`
            callable.

        Returns
        -------
        np.array
            Array of shape (n,) containing time varying distances between to a
            ground location.
        """

        if type(self._ITRS) == type(None):
            self.ITRS()

        ITRSdata = self._ITRS
        groundGEO = np.repeat(np.array(list(groundPos.coor)), self.length, axis=0)
        position = Coordinate(basePos=groundGEO, type="GEO", timeData=self.time)
        groundITRS = np.repeat(position.ITRS, self.length)

        # Find LOS vector norm.
        LOSvec = ITRSdata - groundITRS
        distances = np.linalg.norm(LOSvec, axis=1)

        if kwargs:
            distances = self.interp(distances, **kwargs)

        return distances
    
    def equatorial(self, **kwargs) -> np.array:
        pass
    
    def ecliptic(self, **kwargs) -> np.array:
        pass
    
    def galactic(self, **kwargs) -> np.array:
        pass
    
    def super_galactic(self, **kwargs) -> np.array:
        pass
    
    def sexagesimal(self, angles: np.array) -> np.array:
        """Convert decimal angles into sexagesimal angles.

        Parameters
        ----------
        angles : np.array
            Array of shape (n,) containing angles in decimal degrees.

        Returns
        -------
        np.array
            Array of shape (n,) containing sexagesimal angles as strings.
        """

        length = angles.shape[0]
        sexagesimalAng = np.empty((length,), dtype="<U32")

        for i in range(length):

            num = angles[i]
            sign = "+" if num >= 0 else "-"
            degrees = str(int(abs(num) - num % 1)).zfill(2)
            minutes = str(int(60 * (num % 1) - (60 * (num % 1)) % 1)).zfill(2)
            seconds = "{:.3f}".format(60 * ((60 * (num % 1)) % 1)).zfill(5)
            degree_symbol = u"\u00B0"
            minute_symbol = u"\u2032"
            second_symbol = u"\u2033"
            sexagesimalAng[i] = f"{sign}{degrees}{degree_symbol}{minutes}{minute_symbol}{seconds}{second_symbol}"

        return sexagesimalAng
