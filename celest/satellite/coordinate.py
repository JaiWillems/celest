"""Localize position information representations.

The  coordinate module provieds the `Coordinate` class to allow for a simple
user interface into converting a position type into various representations.
The class is also used for position inputs and outputs for various other
`Celest` functionality.
"""


from celest.core.decorators import set_module
from celest.encounter import GroundPosition
from celest.satellite import Interpolation, Time
from typing import Any, Dict, Literal
import numpy as np


@set_module('celest.satellite')
class Coordinate(object):
    """Localize position information representations.

    The `Coordinate` class provides a simple user interface into converting a
    position type into various representations. The class is also used for
    position inputs and outputs for various other `Celest` functionality.

    Parameters
    ----------
    basePos : np.array
        Base position to initialize the `Coodinate` class.
    type : {"GEO", "ECEF", "ECI"}
        Specifying the inputed position type.
    timeData : Time
        Times associated with the position data. The length of the `timeData`
        parameter must match the length of the `basePos` parameter.
    factor : int, optional
        Interpolation factor used to interpolate the inputted position data.

    Attributes
    ----------
    times : Time
        Times corresponding to the position data.
    interp : Interpolation
        Callable attribute for general interpolation.
    length : int
        Length of the position arrays.

    Methods
    -------
    GEO(**kwargs)
        Return the geographical position representation.
    ERA(**kwargs)
        Return the Earth Rotation Angles.
    ECI(**kwargs)
        Return ECI position data.
    ECEF(**kwargs)
        Return ECEF position data.
    horizontal(groundPos, **kwargs)
        Return horizontal position data.
    off_nadir(groundPos, **kwargs)
        Return the off-nadir angle to a ground location.
    altitude(**kwargs)
        Return the altitude above Earth's surface.
    distance(groundPos, **kwargs)
        Return the distance to a ground location.
    sexagesimal(angles)
        Convert decimal angles into sexagesimal angles.

    Notes
    -----
    Units of the `Coordinate` class are represented in the metric system.
    Specific units will be detailed in method documentation strings.
    """

    def __init__(self, basePos: np.array, type: str, timeData: Time, factor: int=0) -> None:
        """Initialize attributes."""

        self.time = timeData

        self._GEO = None
        self._ERA = None
        self._ECI = None
        self._ECEF = None
        self._altitude = None

        self.interp = Interpolation()

        self.length = None
        self._set_base_position(basePos, type, factor)
    
    def _set_base_position(self, basePos: np.array, type: Literal["ECEF", "ECI","GEO"], factor: int) -> None:
        """Initialize the base position.

        Parameters
        ----------
        basePos : np.array
            Array containing the inputted position data.
        type : {"ECEF", "ECI", "GEO"}
            String identifier of the type of inputted position data.
        factor : int
            Interpolation factor used to interpolate inputed data.
        """

        if factor > 0:
            basePos = self.interp(basePos, factor=factor)

        if type == "GEO":
            if basePos.shape[1] == 2:
                radius = self._WGS84_radius(basePos[:, 0])
                basePos = np.concatenate((basePos, radius), axis=1)
            self._GEO = basePos
        elif type == "ECI":
            self._ECI = basePos
        elif type == "ECEF":
            self._ECEF = basePos
        
        self.length = basePos.shape[0]
    
    def _GEO_to_ECEF(self, geoPos: np.array) -> np.array:
        """Convert from geographical to ECEF coordinates.

        Parameters
        ----------
        geoPos : np.array
            Array of shape (n, 3) containing Geodetic coordinates for a location
            with columns of lattitude, longitude, radius given in decimal
            degrees and km.
        radius : np.array, optional
            Radius of the position from the Earth's surface.
        
        Returns
        -------
        np.array
            Array of shape (n, 3) containing the XYZ ECEF position data.
        """
        radius = geoPos[:, 2]
        
        # Convert geographical to spherical.
        theta = geoPos[:, 1]
        negLon = np.where(theta < 0)
        posLon = np.where(theta > 0)
        theta[negLon] = np.radians(theta[negLon] + 360)
        theta[posLon] = np.radians(theta[posLon])
        phi = np.radians(90 - geoPos[:, 0])

        # Convert spherical to cartesian.
        x = radius * np.cos(theta) * np.sin(phi).reshape(-1, 1)
        y = radius * np.sin(theta) * np.sin(phi).reshape(-1, 1)
        z = radius * np.cos(phi).reshape(-1, 1)
        ECEFpos = np.concatenate((x, y, z), axis=1)

        return ECEFpos

    def _ECEF_to_GEO(self, ECEFpos: np.array) -> np.array:
        """Convert from ECEF to geographical coordinates.

        Parameters
        ----------
        ECEFpos : np.array
            Array of shape (n, 3) with rows containing XYZ ECEF position data.
        
        Returns
        -------
        np.array
            Array of shape (n, 3) with columns of lattitude, longitude, radius
            data in decimal degrees and km.
        """
        # Cartesian coordinates.
        x = ECEFpos[:, 0]
        y = ECEFpos[:, 1]
        z = ECEFpos[:, 2]

        # Convert to spherical coordinates.
        radius = np.sqrt(x**2 + y**2 + z**2)
        theta = np.degrees(np.arccos(np.divide(z, radius)))
        phi = np.degrees(np.arctan(np.divide(y, x)))
        phi[np.where(x < 0)[0]] = phi[np.where(x < 0)[0]] - 180

        # Convert to geographical coordinates.
        lat = 90 - theta
        lon = phi
        lon[np.where(lon < 180)[0]] = phi[np.where(lon < 180)[0]] + 360
        lon[np.where(lon > 180)[0]] = phi[np.where(lon > 180)[0]] - 360

        # Fomulate output array.
        lat = lat.reshape(-1, 1)
        lon = lon.reshape(-1, 1)
        radius = radius.reshape(-1, 1)
        geo = np.concatenate((lat, lon, radius), axis=1)

        return geo

    def GEO(self, **kwargs: Dict[str, Any]) -> np.array:
        """Return the geographical position representation.

        Parameters
        ----------
        kwargs : dict, optional
            Optional keyword arguments passed into the `Interpolation()`
            callable.
        
        Returns
        -------
        np.array
            Array of shape (n, 3) with rows of lattitude, longitude, radius
            position data.
        """
        if type(self._GEO) == type(None):
            if kwargs:
                return self.interp(self._GEO, **kwargs)
            return self._GEO
        if type(self._ECEF) == type(None):
            self.ECEF()

        GEOdata = self._ECEF_to_GEO(ECEFpos=self._ECEF)
        self._GEO = GEOdata

        if kwargs:
            GEOdata = self.interp(GEOdata, **kwargs)

        return GEOdata
    
    def ERA(self, **kwargs: Dict[str, Any]) -> np.array:
        """Return the Earth Rotation Angles.

        Parameters
        ----------
        kwargs : dict, optional
            Optional keyword arguments passed into the `Interpolation()`
            callable.

        Returns
        -------
        np.array
            Array of shape (n,) containing radian earth rotation angles.

        Notes
        -----
        The method implements the earth rotation angle formula:

        .. math:: \gamma = 360.99(\Delta T) + 280.46

        Examples
        --------
        >>> finch = Satellite()
        >>> UTCTimeData = np.array(["2020-06-01 12:00:00.0340", ...,
        ...                         "2020-06-01 12:01:00.0340"])
        >>> ERAangles = finch.ERA(timeData=UTCtimeData)
        """

        angArr = np.zeros((self.length,))
        Julian = self.time.julian

        # Multiply time elapsed since J2000 by ERA and add J2000 orientation.
        dJulian = Julian - 2451545
        angArr = (360.9856123035484 * dJulian + 280.46) % 360

        angArr = np.radians(angArr)
        self.ERAdata = angArr

        if kwargs:
            angArr = self.interp(angArr, **kwargs)

        return angArr
    
    def _ECI_and_ECEF(self, posData: np.array, type: str) -> np.array:
        """Convert between ECI and ECEF positions.

        Parameters
        ----------
        posData : np.array
            Array of shape (n, 3) representing the input data as XYZ cartesian
            data.
        type : {"ECI", "ECEF"}
            The type of inputted data.
        
        Returns
        -------
        np.array
            Array of shape (n, 3) representing the output data as XYZ cartesian
            data.
        """
        if type == "ECI":
            theta = -self._ERA
        else:
            theta = self._ERA
        
        # Construct rotational matrix.
        outVec = np.zeros((self.length, 3))
        A11 = np.cos(theta)
        A12 = -np.sin(theta)
        A21 = np.sin(theta)
        A22 = np.cos(theta)

        # Rotate position data around z-axis by ERA.
        outVec[:, 0] = np.add(np.multiply(A11, self.posData[:, 0]),
                              np.multiply(A12, self.posData[:, 1]))
        outVec[:, 1] = np.add(np.multiply(A21, self.posData[:, 0]),
                              np.multiply(A22, self.posData[:, 1]))
        outVec[:, 2] = self.posData[:, 2]

        return outVec
    
    def ECI(self, **kwargs: Dict[str, Any]) -> np.array:
        """Return ECI position data.

        Parameters
        ----------
        kwargs : dict, optional
            Optional keyword arguments passed into the `Interpolation()`
            callable.

        Returns
        -------
        np.array
            Array of shape (n,3) with columns of X, Y, Z ECI position data.

        See Also
        --------
        ECEF : Return ECEF position data.
        """

        if type(self._ERA) == type(None):
            self.ERA()
        if type(self._ECEF) == type(None):
            self.ECEF()

        ECIdata = self._ECI_and_ECEF(posData=self._ECEF, type="ECEF")
        self._ECI = ECIdata

        if kwargs:
            ECIdata = self.interp(ECIdata, **kwargs)

        return ECIdata
    
    def ECEF(self, **kwargs: Dict[str, Any]) -> np.array:
        """Return ECEF position data.

        Parameters
        ----------
        kwargs : dict, optional
            Optional keyword arguments passed into the `Interpolation()`
            callable.

        Returns
        -------
        np.array
            Array of shape (n,3) with columns of X, Y, Z ECEF position data.

        See Also
        --------
        ECI : Return ECI position data.
        """

        if type(self._ECI) == type(None) and type(self._GEO) == type(None):
            if type == "ECI":
                self.ECI()

            elif type == "GEO":
                self.GEO()

        if type(self._GEO) != type(None):
            ECEFdata = self._GEO_to_ECEF(posData=self._GEO)

        elif type(self._ECI) != type(None):
            if type(self._ERA) == type(None):
                self.ERA()
            ECEFdata = self._ECI_and_ECEF(posData=self.ECIdata, type="ECI")

        if kwargs:
            ECEFdata = self.interp(ECEFdata, **kwargs)

        return self.ECEFdata
    
    def _get_ang(self, vecOne: np.array, vecTwo: np.array) -> float:
        """Calculate degree angle bewteen two vectors.
        
        Parameters
        ----------
        vecOne, vecTwo : np.array
            Arrays of shape (n,3) with rows of ECEF data.
        
        Returns
        -------
        float
            Degree angle between the two arrays.
        """
        # Use simple linalg formula.
        dividend = np.einsum("ij, ij->i", vecOne, vecTwo)
        divisor = np.multiply(np.linalg.norm(
            vecOne, axis=1), np.linalg.norm(vecTwo, axis=1))
        arg = np.divide(dividend, divisor)
        ang = np.degrees(np.arccos(arg))

        return ang
    
    def horizontal(self, groundPos: GroundPosition, **kwargs: Dict[str, Any]) -> np.array:
        """Return horizontal position data.

        Parameters
        ----------
        groundPos : GroundPosition object
            GroundPosition object instantiated as per its documentation.
        kwargs : dict, optional
            Optional keyword arguments passed into the `Interpolation()`
            callable.

        Returns
        -------
        tuple
            The altitude/azimuth data is returned in a tuple where both the
            altitude and azimuth are arrays of shape (n,).
        """

        if type(self._ECEF) == type(None):
            self.ECEF()

        # Convert observer position into spherical then cartesian.
        obsCoor = groundPos.coor
        radius = groundPos.radius
        xyzObs = np.full((self.length, 3), self._geo_to_ECEF(obsCoor, radius))

        # Determine line of sight vector then altitude.
        ECEFvec = self.ECEFdata
        xyzLOS = np.subtract(ECEFvec, xyzObs)
        Alt = 90 - self._get_ang(xyzLOS, xyzObs)

        # Find surface tangent vector passing through z-axis.
        kHat = np.full((self.length, 3), np.array([0, 0, 1]))
        beta = np.pi/2 - np.radians(90 - obsCoor[0])
        tangentVec = np.subtract((kHat.T * radius/np.sin(beta)).T, xyzObs)

        # Find LOS projection on tangent plane.
        coeff = np.einsum("ij, ij->i", xyzLOS, xyzObs)/radius**2
        normProj = (xyzObs.T * coeff).T
        projLOS = np.subtract(xyzLOS, normProj)

        # Determing azimuth.
        vecOne = np.cross(tangentVec, xyzObs)
        normOne = 1/np.linalg.norm(vecOne, axis=1).reshape((self.length, 1))
        vecOneUnit = normOne*vecOne
        vecTwo = (vecOneUnit.T * np.einsum("ij, ij->i", projLOS, vecOneUnit)).T
        normTwo = 1/np.linalg.norm(vecTwo, axis=1).reshape((self.length, 1))
        vecTwoUnit = normTwo*vecTwo

        posEq = np.isclose(vecOneUnit, vecTwoUnit, atol=0.01)
        negEq = np.isclose(vecOneUnit, -vecTwoUnit, atol=0.01)
        posInd = np.where(np.all(posEq, axis=1))[0]
        negInd = np.where(np.all(negEq, axis=1))[0]

        Az = np.zeros((self.length,))
        Az[posInd] = self._get_ang(tangentVec[posInd], projLOS[posInd])
        Az[negInd] = 360 - self._get_ang(tangentVec[negInd], projLOS[negInd])

        if kwargs:
            Alt = self.interp(Alt, **kwargs)
            Az = self.interp(Az, **kwargs)

        return Alt, Az
    
    def off_nadir(self, groundPos: GroundPosition, **kwargs: Dict[str, Any]) -> np.array:
        """Return the off-nadir angle to a ground location.

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
            Array of shape (n,) with off-nadir angle data in decimal degrees.
        """

        if type(self._ECEF) == type(None):
            self.ECEF()

        obsCoor = groundPos.coor
        radius = groundPos.radius

        geoPos = self._geo_to_ECEF(np.array(list(obsCoor).append(radius)))[0]

        xyzObs = np.repeat(geoPos, self.length, 0)

        ECEFvec = self._ECEF
        xyzLOS = np.subtract(ECEFvec, xyzObs)

        ang = self._get_ang(xyzLOS, ECEFvec)

        if kwargs:
            ang = self.interp(ang, **kwargs)

        return ang

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
    
    def altitude(self, **kwargs: Dict[str, Any]) -> np.array:
        """Return the altitude above Earth's surface.

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

        if type(self._ECEF) == type(None):
            self.ECEF()

        ECEFdata = self._ECEF
        x = ECEFdata[:, 0]
        y = ECEFdata[:, 1]
        z = ECEFdata[:, 2]
        arg = np.sqrt(x**2 + y**2) / z
        lattitude = np.arctan(arg)

        earthRadius = self._WGS84_radius(lattitude)
        satRadius = np.linalg.norm(ECEFdata, axis=1)

        altitude = satRadius - earthRadius
        self._altitude = altitude

        if kwargs:
            altitude = self.interp(altitude, **kwargs)

        return altitude
    
    def distance(self, groundPos: GroundPosition, **kwargs: Dict[str, Any]) -> np.array:
        """Return the distance to a ground location.

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

        if type(self._ECEF) == type(None):
            self.ECEF()

        ECEFdata = self._ECEF
        groundGEO = np.repeat(np.array(list(groundPos.coor)), self.length, axis=0)
        position = Coordinate(basePos=groundGEO, type="GEO", timeData=self.time)
        groundECEF = np.repeat(position.ECEF, self.length)

        # Find LOS vector norm.
        LOSvec = ECEFdata - groundECEF
        distances = np.linalg.norm(LOSvec, axis=1)

        if kwargs:
            distances = self.interp(distances, **kwargs)

        return distances
    
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
