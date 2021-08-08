"""Localize position information representations.

The  coordinate module provieds the `Coordinate` class to allow for a simple
user interface into converting a position type into various representations.
The class is also used for position inputs and outputs for various other
`Celest` functionality.
"""


from celest.core.decorators import set_module
from celest.encounter import GroundPosition
from celest.satellite import Interpolation
from typing import Any, Dict, Literal, Tuple, TYPE_CHECKING
import numpy as np

if TYPE_CHECKING:
    from celest.satellite import Time

@set_module('celest.satellite')
class Coordinate(Interpolation):
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
    equatorial(groundPos, **kwargs)
        Return Equatorial coordinates.
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
            basePos = self._interp(basePos, factor=factor)

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
        neg_lon = np.where(theta < 0)
        pos_lon = np.where(theta > 0)
        theta[neg_lon] = np.radians(theta[neg_lon] + 360)
        theta[pos_lon] = np.radians(theta[pos_lon])
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
        kwargs : Dict, optional
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
                return self._interp(self._GEO, **kwargs)
            return self._GEO
        if type(self._ECEF) == type(None):
            self.ECEF()

        GEO_data = self._ECEF_to_GEO(ECEFpos=self._ECEF)
        self._GEO = GEO_data

        if kwargs:
            GEO_data = self._interp(GEO_data, **kwargs)

        return GEO_data
    
    def ERA(self, **kwargs: Dict[str, Any]) -> np.array:
        """Return the Earth Rotation Angles.

        Parameters
        ----------
        kwargs : Dict, optional
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

        ang = np.zeros((self.length,))
        Julian = self.time.julian

        # Multiply time elapsed since J2000 by ERA and add J2000 orientation.
        dJulian = Julian - 2451545
        ang = (360.9856123035484 * dJulian + 280.46) % 360

        ang = np.radians(ang)
        self._ERA = ang

        if kwargs:
            ang = self._interp(ang, **kwargs)

        return ang
    
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

        if type(self._ERA) == type(None):
            self.ERA()

        if type == "ECI":
            theta = -self._ERA
        else:
            theta = self._ERA
        
        # Construct rotational matrix.
        A11 = np.cos(theta)
        A12 = -np.sin(theta)
        A21 = np.sin(theta)
        A22 = np.cos(theta)

        # Rotate position data around z-axis by ERA.
        output = np.zeros((self.length, 3))
        output[:, 0] = A11 * posData[:, 0] + A12 * posData[:, 1]
        output[:, 1] = A21 * posData[:, 0] + A22 * posData[:, 1]
        output[:, 2] = posData[:, 2]

        return output
    
    def ECI(self, **kwargs: Dict[str, Any]) -> np.array:
        """Return ECI position data.

        Parameters
        ----------
        kwargs : Dict, optional
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

        if type(self._ECEF) == type(None):
            self.ECEF()

        ECI_data = self._ECI_and_ECEF(posData=self._ECEF, type="ECEF")
        self._ECI = ECI_data

        if kwargs:
            ECI_data = self._interp(ECI_data, **kwargs)

        return ECI_data
    
    def ECEF(self, **kwargs: Dict[str, Any]) -> np.array:
        """Return ECEF position data.

        Parameters
        ----------
        kwargs : Dict, optional
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
            ECEF_data = self._GEO_to_ECEF(posData=self._GEO)

        elif type(self._ECI) != type(None):
            ECEF_data = self._ECI_and_ECEF(posData=self._ECI, type="ECI")

        if kwargs:
            ECEF_data = self._interp(ECEF_data, **kwargs)

        return self.ECEF_data
    
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
        dividend = np.sum(vecOne, vecTwo, axis=1)
        divisor = np.linalg.norm(vecOne, axis=1) * np.linalg.norm(vecTwo, axis=1)
        arg = np.divide(dividend, divisor)
        ang = np.degrees(np.arccos(arg))

        return ang
    
    def horizontal(self, groundPos: GroundPosition, **kwargs: Dict[str, Any]) -> np.array:
        """Return horizontal position data.

        Parameters
        ----------
        groundPos : GroundPosition object
            GroundPosition object instantiated as per its documentation.
        kwargs : Dict, optional
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

        # Convert observer position into cartesian coordinates.
        obs_coor, radius = groundPos.coor, groundPos.radius
        GEO_data = np.array([obs_coor[0], obs_coor[1], radius]).reshape((1, -1))
        obs = np.repeat(self._GEO_to_ECEF(GEO_data), self.length, axis=0)

        # Determine line of sight vector then altitude.
        LOS = self._ECEF - obs
        Alt = 90 - self._get_ang(LOS, obs)

        # Find surface tangent vector passing through z-axis.
        k_hat = np.repeat(np.array([[0, 0, 1]]), self.length, axis=0)
        beta = np.radians(obs_coor[0])
        tangent = (k_hat.T * radius / np.sin(beta)).T - obs

        # Find LOS projection on tangent plane.
        coeff = np.sum(LOS * obs, axis=1) / radius**2
        norm_proj = (obs.T * coeff).T
        proj_LOS = LOS - norm_proj

        # Determing azimuth.
        reference = np.cross(tangent, obs)
        neg_ind = np.where(np.sum(proj_LOS * reference, axis=1) < 0)[0]
        Az = self._get_ang(tangent, proj_LOS)
        Az[neg_ind] = 360 - self._get_ang(tangent[neg_ind], proj_LOS[neg_ind])

        if kwargs:
            Alt = self._interp(Alt, **kwargs)
            Az = self._interp(Az, **kwargs)

        return Alt, Az
    
    def off_nadir(self, groundPos: GroundPosition, **kwargs: Dict[str, Any]) -> np.array:
        """Return the off-nadir angle to a ground location.

        Parameters
        ----------
        groundPos : GroundPosition
            `GroundPosition` object instantiated as per its documentation.
        kwargs : Dict, optional
            Optional keyword arguments passed into the `Interpolation()`
            callable.

        Returns
        -------
        np.array
            Array of shape (n,) with off-nadir angle data in decimal degrees.
        """

        if type(self._ECEF) == type(None):
            self.ECEF()

        obs_coor = groundPos.coor
        radius = groundPos.radius

        GEO_data = self._GEO_to_ECEF(np.array([obs_coor[0], obs_coor[1], radius]))
        obs = np.repeat(GEO_data, self.length, 0)

        ECEF_data = self._ECEF
        LOS = np.subtract(ECEF_data, obs)

        ang = self._get_ang(LOS, ECEF_data)

        if kwargs:
            ang = self._interp(ang, **kwargs)

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
        semi_major = 6378.137**2
        semi_minor = 6356.752314245**2

        numerator = semi_major * semi_minor
        denominator = semi_major * np.sin(phi)**2 + semi_minor * np.cos(phi)**2
        radius = np.sqrt(numerator / denominator)

        return radius
    
    def altitude(self, **kwargs: Dict[str, Any]) -> np.array:
        """Return the altitude above Earth's surface.

        Parameters
        ----------
        kwargs : Dict, optional
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
                return self._interp(self._altitude)
            else:
                return self._altitude

        if type(self._ECEF) == type(None):
            self.ECEF()

        ECEF_data = self._ECEF
        x = ECEF_data[:, 0]
        y = ECEF_data[:, 1]
        z = ECEF_data[:, 2]
        arg = np.sqrt(x**2 + y**2) / z
        lattitude = np.arctan(arg)

        earth_radius = self._WGS84_radius(lattitude)
        sat_radius = np.linalg.norm(ECEF_data, axis=1)

        altitude = sat_radius - earth_radius
        self._altitude = altitude

        if kwargs:
            altitude = self._interp(altitude, **kwargs)

        return altitude
    
    def distance(self, groundPos: GroundPosition, **kwargs: Dict[str, Any]) -> np.array:
        """Return the distance to a ground location.

        Parameters
        ----------
        groundPos : GroundPosition
            `GroundPosition` object instantiated as per its documentation.
        kwargs : Dict, optional
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

        ECEF_data = self._ECEF
        gnd_GEO = np.repeat(np.array(list(groundPos.coor)), self.length, axis=0)
        position = Coordinate(basePos=gnd_GEO, type="GEO", timeData=self.time)
        gnd_ECEF = np.repeat(position.ECEF, self.length)

        # Find LOS vector norm.
        LOS = ECEF_data - gnd_ECEF
        distances = np.linalg.norm(LOS, axis=1)

        if kwargs:
            distances = self._interp(distances, **kwargs)

        return distances
    
    def equatorial(self, groundPos: GroundPosition, **kwargs: Dict[str, Any]) -> Tuple:
        """Return Equatorial coordinates.

        Parameters
        ----------
        groundPos : GroundPosition
            `GroundPosition` object instantiated as per its documentation.
        kwargs : Dict, optional
            Optional keyword arguments passed into the `Interpolation()`
            callable.
        
        Returns
        -------
        tuple
            Return a tuple of Numpy arrays containing the right ascension and
            declination of the coordinate object.
        """

        latitude, longitude = groundPos.coor[0], groundPos.coor[1]
        alt, az = self.horizontal(groundPos)

        latitude, alt, az = np.radians(latitude), np.radians(alt), np.radians(az + 180)

        tanH = np.sin(az) / (np.cos(az) * np.sin(latitude) + np.tan(alt) * np.cos(latitude))
        sinDelta = np.sin(latitude) * np.sin(alt) - np.cos(latitude) * np.cos(alt) * np.cos(az)

        delta = np.arcsin(sinDelta)
        alpha = self.time.LMST(longitude) - np.arctan(tanH)

        if kwargs:
            delta = self._interp(delta, **kwargs)
            alpha = self._interp(alpha, **kwargs)

        return alpha, delta
    
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
        sexagesimal_ang = np.empty((length,), dtype="<U32")

        for i in range(length):

            num = angles[i]
            sign = "+" if num >= 0 else "-"
            degrees = str(int(abs(num) - num % 1)).zfill(2)
            minutes = str(int(60 * (num % 1) - (60 * (num % 1)) % 1)).zfill(2)
            seconds = "{:.3f}".format(60 * ((60 * (num % 1)) % 1)).zfill(5)
            degree_symbol = u"\u00B0"
            minute_symbol = u"\u2032"
            second_symbol = u"\u2033"
            sexagesimal_ang[i] = f"{sign}{degrees}{degree_symbol}{minutes}{minute_symbol}{seconds}{second_symbol}"

        return sexagesimal_ang
