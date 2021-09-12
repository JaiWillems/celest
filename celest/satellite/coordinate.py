"""Localize position information representations.

This module provides the `Coordinate` class to allow for a simple interface
for converting a position type into various representations. The class is also
used for position inputs and outputs for other `Celest` functionality.
"""


from celest.core.decorators import set_module
from celest.satellite import Interpolation
from typing import Any, Dict, Literal, Tuple
import numpy as np


@set_module('celest.satellite')
class Coordinate(Interpolation):
    """Localize position information representations.

    The `Coordinate` class provides a simple user interface for converting a
    position type into various representations. The class is also used for
    position inputs and outputs for other `Celest` functionality.

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
        Interpolation factor used to interpolate the input position data.

    Attributes
    ----------
    times : Time
        Times corresponding to the position data.
    length : int
        Length of the input position and time arrays.

    Methods
    -------
    GEO(iso=False, **kwargs)
        Return geographical position data.
    ERA(**kwargs)
        Return the Earth rotation angles in radians and decimals.
    ECI(**kwargs)
        Return cartesian ECI position data.
    ECEF(**kwargs)
        Return cartesian ECEF position data.
    horizontal(groundPos, **kwargs)
        Return horizontal position data in degrees and decimals.
    off_nadir(groundPos, **kwargs)
        Return the off-nadir angle to a ground location.
    altitude(**kwargs)
        Return the altitude above Earth's surface in kilometres.
    distance(groundPos, **kwargs)
        Return the distance to a ground location.
    sexagesimal(angles)
        Convert decimal angles into sexagesimal angles.
    """

    def __init__(self, basePos: np.array, type: str, timeData: Any, factor:
                 int=0) -> None:
        """Initialize attributes."""

        self.times = timeData

        self._GEO = None
        self._ECI = None
        self._ECEF = None

        self.length = None
        self._set_base_position(basePos, type, factor)

    def _set_base_position(self, basePos: np.array, type:
                           Literal["ECEF", "ECI", "GEO"], factor: int) -> None:
        """Initialize base position.

        This method takes an input position to initialize the object's base
        position.

        Parameters
        ----------
        basePos : np.array
            Array of shape (n, 2) or (n, 3) containing the inputted position
            data.
        type : {"ECEF", "ECI", "GEO"}
            String defining the type of input position data as either
            Earth-centered-Earth-fixed (ECEF), Earth-centered-inertial (ECI),
            or geographical (GEO) data.
        factor : int
            Interpolation factor used to interpolate input data of length `n`
            to length `factor*n`.

        Notes
        -----
        The input data must be of shape (n, 3) if `type="ECEF"` or `type="ECI"`
        where the columns are XYZ cartesian data. The data can be of the shape
        (n, 2) or (n, 3) if `type="GEO"` where the columns are geodetic
        latitude, terrestrial longitude, and geodetic altitude. When
        geographical data is entered of shape (n, 2), the height data is
        assumed to be zero.
        """

        if factor > 0:
            basePos = self._interp(basePos, factor=factor)

        self.length = basePos.shape[0]

        if type == "GEO":
            if basePos.shape[1] == 2:
                height = np.zeros((self.length, 1))
                basePos = np.concatenate((basePos, height), axis=1)

            self._GEO = basePos

        elif type == "ECI":
            self._ECI = basePos

        elif type == "ECEF":
            self._ECEF = basePos

    def _GEO_to_ECEF(self, posData: np.array) -> np.array:
        """Convert geographical to ECEF coordinates.

        Parameters
        ----------
        posData : np.array
            Array of shape (n, 3) containing geographical coordinates of a
            position with columns of geodetic latitude, terrestrial longitude,
            and geodetic altitude given in decimal degrees and kilometres.

        Returns
        -------
        np.array
            Array of shape (n, 3) containing the XYZ ECEF position data.

        See Also
        --------
        _ECEF_to_GEO : Convert ECEF to geographical coordinates.

        Notes
        -----
        This method uses an ellipsoid based model of the Earth to convert a
        geographical position to ECEF cartesian coordinates using the methods
        described in "Coordinate Systems in Geodesy" by E. J. Krakiwsky and
        D.E. Wells as presented by Christopher Lum.[1]_[2]_

        References
        ----------
        .. [1] E. J. Krakiwsky and D. E. Wells. Coordinate Systems in Geodesy.
           Jan. 1998.
        .. [2] Christopher Lum. Geodetic Coordinates: Computing Latitude and
           Longitude.June 2020.url:https://www.youtube.com/watch?v=4BJ-GpYbZlU.
        """

        a = 6378.137
        b = 6356.752314245

        lat, lon = np.radians(posData[:, 0]), np.radians(posData[:, 1])

        e = np.sqrt(1 - b ** 2 / a ** 2)
        N = a / np.sqrt(1 - e ** 2 * np.sin(lat) ** 2)
        h = posData[:, 2]

        x = ((N + h) * np.cos(lat) * np.cos(lon)).reshape((-1, 1))
        y = ((N + h) * np.cos(lat) * np.sin(lon)).reshape((-1, 1))
        z = ((N * (1 - e ** 2) + h) * np.sin(lat)).reshape((-1, 1))
        ECEF = np.concatenate((x, y, z), axis=1)

        return ECEF

    def _ECEF_to_GEO(self, posData: np.array) -> np.array:
        """Convert ECEF to geographical coordinates.

        Parameters
        ----------
        posData : np.array
            Array of shape (n, 3) with rows containing XYZ ECEF position data.

        Returns
        -------
        np.array
            Array of shape (n, 3) with columns of geodetic latitude,
            terrestrial longitude, and geodetic altitude data in degrees and
            kilometres.

        See Also
        --------
        _GEO_to_ECEF : Convert geographical to ECEF coordinates.

        Notes
        -----
        Let :math:`x`, :math:`y`, and :math:`z` be the ECEF vector components.
        We can then calculate the latitude, :math:`\phi`, using the following
        equation:

        .. math:: \phi = 90^\circ - \cos^{-1}\left(\frac{z}{R_\bigoplus}\right)

        where :math:`R_\bigoplus` is the geocentric radius of the Earth. We can
        also calculate the longitude, :math:`\lambda` using the following:

        .. math::\lambda = \tan^{-1}_2\left(y, x\right)
        """

        # Cartesian coordinates.
        x = posData[:, 0]
        y = posData[:, 1]
        z = posData[:, 2]

        # Convert to spherical coordinates.
        radius = np.sqrt(x ** 2 + y ** 2 + z ** 2)
        theta = np.degrees(np.arccos(z / radius))
        phi = np.degrees(np.arctan2(y, x))

        # Get Geodetic altitude.
        alt = self.altitude()

        # Fomulate output array.
        lat = (90 - theta).reshape(-1, 1)
        lon = phi.reshape(-1, 1)
        alt = alt.reshape(-1, 1)
        geo = np.concatenate((lat, lon, alt), axis=1)

        return geo

    def _ISO6709_representation(self, geoPos: np.array) -> np.array:
        """Format geographical data in an international standard.

        This method formats the input geographical point location data in the
        ISO6709 standard format.

        Parameters
        ----------
        geoPos : np.array
            Array of shape (n, 3) with columns of geodetic latitude,
            terrestrial longitude, and geodetic altitude data.

        Returns
        -------
        np.array
            Array of shape (n,) containing position strings in accordance to
            ISO6709 standard.
        """

        deg = "\u00B0"
        min = "\u2032"
        sec = "\u2033"

        out_arr = np.empty((geoPos.shape[0],), dtype="<U37")

        for i in range(geoPos.shape[0]):

            lat, lon, h = geoPos[i, 0], geoPos[i, 1], geoPos[i, 2]

            # Format latitude string.
            degree = int(abs(lat))
            minute = int(60 * np.round(abs(lat) - degree, 2))
            second = np.round(abs(60 * (60 * (abs(lat) - degree) - minute)), 2)
            direction = "N" if lat >= 0 else "S"

            degree = "%.2s" % str(degree).zfill(2)
            minute = "%.2s" % str(minute).zfill(2)
            second = str("%.2f" % second).zfill(5)

            lat_str = f"{degree}{deg}{minute}{min}{second}{sec}{direction} "

            # Format longitude string.
            degree = int(abs(lon))
            minute = int(60 * np.round(abs(lon) - degree, 2))
            second = np.round(abs(60 * (60 * (abs(lon) - degree) - minute)), 2)
            direction = "E" if lon >= 0 else "W"

            degree = "%.3s" % str(degree).zfill(2)
            minute = "%.2s" % str(minute).zfill(2)
            second = str("%.2f" % second).zfill(5)

            lon_str = f"{degree}{deg}{minute}{min}{second}{sec}{direction} "

            # Format height string.
            h_str = f"{h}km"

            out_arr[i] = lat_str + lon_str + h_str

        return out_arr

    def GEO(self, iso=False, **kwargs: Dict[str, Any]) -> np.array:
        """Return geographical position data.

        Parameters
        ----------
        iso : bool, optional
            Formats the output as ISO6709 sexagesimal position strings if true.
        kwargs : Dict, optional
            Optional keyword arguments passed into the
            `Interpolation()._interp()` method.

        Returns
        -------
        np.array
            Array of shape (n, 3) with columns of geodetic latitude,
            terrestrial longitude, and geodetic altitude data. If
            `iso=True`, the output is an array of shape (n,) containing
            standard representation position strings.

        See Also
        --------
        ECEF : Return ECEF position data.
        ECI : Return ECI position data.

        Examples
        --------
        >>> timeData = Time(julian=np.array([2454545]))
        >>> basePos = np.array([[6343.82, -2640.87, -11.26]])
        >>> coor = Coordinate(basePos=basePos, type="ECEF", timeData=timeData)
        >>> coor.GEO()
        np.array([[-9.38870528e-02, -2.26014826e+01, 5.04126976e+02]])

        We can generate user-friendly position strings by setting `iso=True`.

        >>> coor.GEO(iso=True)
        np.array(['00°05′37.99″S 22°36′05.34″W 504.12697'])
        """

        if self._GEO is None:
            if self._ECEF is not None:
                self._GEO = self._ECEF_to_GEO(posData=self._ECEF)
            else:
                self._ECEF = self._ECI_and_ECEF(posData=self._ECI, type="ECI")
                self._GEO = self._ECEF_to_GEO(posData=self._ECEF)

        GEO = self._interp(self._GEO, **kwargs) if kwargs else self._GEO

        if iso:
            GEO = self._ISO6709_representation(geoPos=GEO)

        return GEO

    def ERA(self, **kwargs: Dict[str, Any]) -> np.array:
        """Return the Earth rotation angles in radians and decimals.

        Parameters
        ----------
        kwargs : Dict, optional
            Optional keyword arguments passed into the
            `Interpolation()._interp()` method.

        Returns
        -------
        np.array
            Array of shape (n,) containing Earth rotation angles in radians and
            decimals.

        Notes
        -----
        The Earth rotation angle denotes the diurnal rotation component of the
        coordinate misalignment between the ECI and ECEF frames. The angle can
        be calculated from the following formulation:

        .. math:: \gamma^\circ = 360.9856123035484\Delta T + 280.46

        where :math:`\Delta T=JD-2451545` is the elapsed days since the J2000
        epoch where :math:`JD` is the Julian day.[1]_

        References
        ----------
        .. [1] Don Koks. Changing Coordinates in the Context of Orbital
           Mechanics. Cyber and Electronic Warfare Division, Defence Science,
           and Technology Group, Jan.2017, p. 12 - 13.

        Examples
        --------
        >>> timeData = Time(julian=np.array([2454545]))
        >>> basePos = np.array([[6343.82, -2640.87, -11.26]])
        >>> posData = Coordinate(basePos=basePos, type="ECEF", timeData=timeData)
        >>> posData.ERA()
        np.array([6.2360075])
        """

        ang = np.zeros((self.length,))
        jul_data = self.times.julian()

        # Multiply time elapsed since J2000 by ERA and add J2000 orientation.
        dJulian = jul_data - 2451545
        ang = (360.9856123035484 * dJulian + 280.46) % 360

        ang = np.radians(ang)

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
            The type of input data.

        Returns
        -------
        np.array
            Array of shape (n, 3) representing the output data as XYZ cartesian
            data.

        Notes
        -----
        This method applies a simplification to converting between ECI and ECEF
        coordinate by disregarding precession, nutation, and polar motion
        effects (which will be incorporated in future releases). By disregarding
        such factors, we can assume the alignment of the ECI and ECEF axes and
        simplify the conversion to a rotation about the z-axis.

        We can convert between the ECI and ECEF frames using the following
        equations:

        .. math:: \vec{V}_{ECI} = E_3^\gamma\;\vec{V}_{ECEF}

        .. math:: \vec{V}_{ECEF} = E_3^{-\gamma}\;\vec{V}_{ECI}

        where :math:`\gamma` is the Earth rotation angle, :math:`E_3^\theta` is
        the rotation matrix that rotates a vector around the z-axis by an angle
        :math:`\theta`, and :math:`\vec{V}` is a position vector in either the
        ECI or ECEF frames.[1]_

        References
        ----------
        .. [1] Don Koks. Changing Coordinates in the Context of Orbital
           Mechanics. Cyber and Electronic Warfare Division, Defence Science,
           and Technology Group, Jan.2017, p. 21.
        """

        theta = -self.ERA() if type == "ECI" else self.ERA()

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
        """Return cartesian ECI position data.

        Parameters
        ----------
        kwargs : Dict, optional
            Optional keyword arguments passed into the
            `Interpolation()._interp()` method.

        Returns
        -------
        np.array
            Array of shape (n,3) with columns of XYZ ECI position data.

        See Also
        --------
        ECEF : Return cartesian ECEF position data.
        GEO : Return geographical position data.

        Examples
        --------
        >>> timeData = Time(julian=np.array([2454545]))
        >>> basePos = np.array([[6343.82, -2640.87, -11.26]])
        >>> coor = Coordinate(basePos=basePos, type="ECEF", timeData=timeData)
        >>> coor.ECI()
        np.array([[6212.21719598, -2937.10811161, -11.26]])
        """

        if self._ECI is None:
            if self._ECEF is None:
                self._ECEF = self._GEO_to_ECEF(self._GEO)
            self._ECI = self._ECI_and_ECEF(self._ECEF, type="ECEF")

        ECI = self._interp(self.ECI, **kwargs) if kwargs else self._ECI

        return ECI

    def ECEF(self, **kwargs: Dict[str, Any]) -> np.array:
        """Return cartesian ECEF position data.

        Parameters
        ----------
        kwargs : Dict, optional
            Optional keyword arguments passed into the
            `Interpolation()._interp()` method.

        Returns
        -------
        np.array
            Array of shape (n,3) with columns of XYZ ECEF position data.

        See Also
        --------
        ECI : Return cartesian ECI position data.
        GEO : Return geographical position data.

        Examples
        --------
        >>> timeData = Time(julian=np.array([2454545]))
        >>> basePos = np.array([[6343.82, -2640.87, -11.26]])
        >>> coor = Coordinate(basePos=basePos, type="ECI", timeData=timeData)
        >>> coor.ECEF()
        np.array([[6461.30569276, -2338.75507354, -11.26]])
        """

        if self._ECEF is None:
            if self._ECI is not None:
                self._ECEF = self._ECI_and_ECEF(posData=self._ECI, type="ECI")
            else:
                self._ECEF = self._GEO_to_ECEF(posData=self._GEO)

        ECEF = self._interp(self._ECEF, **kwargs) if kwargs else self._ECEF

        return ECEF

    def _get_ang(self, vecOne: np.array, vecTwo: np.array) -> float:
        """Calculate degree angle bewteen two vectors.

        Parameters
        ----------
        vecOne, vecTwo : np.array
            Arrays of shape (n,3) with rows of XYZ cartesian data.

        Returns
        -------
        float
            Degree angle between the two arrays.
        """

        num = np.sum(vecOne * vecTwo, axis=1)
        denom = np.linalg.norm(vecOne, axis=1) * np.linalg.norm(vecTwo, axis=1)
        ang = np.degrees(np.arccos(num / denom))

        return ang

    def horizontal(self, groundPos: Any, **kwargs: Dict[str, Any]) -> Tuple:
        """Return horizontal position data in degrees and decimals.

        In this method, the azimuth is calculated as increasing degrees
        clockwise from North and ranges from 0 to 360 degrees. The altitude is
        calculated as increasing degrees above the observer's local horizon and
        ranges from 0 to 90 degrees.

        Parameters
        ----------
        groundPos : GroundPosition
            Ground location defining the center of the horizontal system.
        kwargs : Dict, optional
            Optional keyword arguments passed into the
            `Interpolation()._interp()` method.

        Returns
        -------
        tuple
            Tuple of length two containing the altitude and azimuth data as
            decimals and degrees in NumPy arrays of shape (n,).

        Examples
        --------
        >>> timeData = Time(julian=np.array([2454545]))
        >>> basePos = np.array([[6343.82, -2640.87, -11.26]])
        >>> groundPos = GroundPosition(name="Saskatoon", coor=obsCoor)
        >>> coor = Coordinate(basePos=basePos, type="ECEF", timeData=timeData)
        >>> coor.horizontal(groundPos=groundPos)
        (array([-40.8786098]), array([94.73615482]))
        """

        if self._ECEF is None:
            self.ECEF()

        # Convert observer position into cartesian coordinates.
        coor, radius = groundPos.coor, groundPos.radius
        GEO_data = self._GEO_to_ECEF(np.array([[coor[0], coor[1], 0]]))
        obs = np.repeat(GEO_data, self.length, 0)

        # Determine line of sight vector then altitude.
        LOS = self._ECEF - obs
        Alt = 90 - self._get_ang(LOS, obs)

        # Find surface tangent vector passing through z-axis.
        k_hat = np.repeat(np.array([[0, 0, 1]]), self.length, axis=0)
        beta = np.radians(coor[0])
        tangent = (k_hat.T * radius / np.sin(beta)).T - obs

        # Find LOS projection on tangent plane.
        norm_proj = (obs.T * np.sum(LOS * obs, axis=1) / radius ** 2).T
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

    def off_nadir(self, groundPos: Any, **kwargs: Dict[str, Any]) -> np.array:
        """Return the off-nadir angle to a ground location.

        This method calculates the off-nadir angle, in degrees and decimals, to
        the input ground position at each satellite position.

        Parameters
        ----------
        groundPos : GroundPosition
            Calculate off-nadir angles for the specified ground location.
        kwargs : Dict, optional
            Optional keyword arguments passed into the
            `Interpolation()._interp()` method.

        Returns
        -------
        np.array
            Array of shape (n,) containing off-nadir angle data in degrees and
            decimals.

        Notes
        -----
        The off-nadir angle is the acute angle measured in increasing degrees
        from the satellites nadir to the line joining the satellite to the
        ground location of interest.
        """

        if self._ECEF is None:
            self.ECEF()

        coor = groundPos.coor
        GEO_data = self._GEO_to_ECEF(np.array([[coor[0], coor[1], 0]]))
        obs = np.repeat(GEO_data, self.length, 0)

        LOS = np.subtract(self._ECEF, obs)
        ang = self._get_ang(LOS, self._ECEF)

        if kwargs:
            ang = self._interp(ang, **kwargs)

        return ang

    def _WGS84_radius(self, lattitude: np.array) -> np.array:
        """Calculate the Earth's geocentric radius using WGS84.

        Parameters
        ----------
        latitude : np.array
            Array of shape (n,) representing the geocentric latitude of a
            ground location in degrees and decimals.

        Returns
        -------
        np.array
            Array of shape (n,) containing the Earth's geocentric radius in
            kilometres.

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
        phi = np.radians(lattitude)

        # Define WGS84 Parameters.
        a = 6378.137
        b = 6356.752314245

        c_phi, s_phi = np.cos(phi), np.sin(phi)

        num = (a ** 2 * c_phi) ** 2 + (b ** 2 * s_phi) ** 2
        denom = (a * c_phi) ** 2 + (b * s_phi) ** 2
        radius = np.sqrt(num / denom)

        return radius

    def altitude(self, **kwargs: Dict[str, Any]) -> np.array:
        """Return the geodetic altitude above Earth's surface in kilometres.

        This method uses the WGS84 reference ellipsoid to calculate the
        geodetic altitude above the Earth's surface.

        Parameters
        ----------
        kwargs : Dict, optional
            Optional keyword arguments passed into the
            `Interpolation()._interp()` method.

        Returns
        -------
        np.array
            Array of shape (n,) containing satellite altitudes in kilometres.

        Notes
        -----
        This method uses an ellipsoid based model of the Earth to calculate the
        ellipsoid height in an iterative manner described in "Coordinate
        Systems in Geodesy" by E. J. Krakiwsky and D.E. Wells.[1]_

        References
        ----------
        .. [1] E. J. Krakiwsky and D. E. Wells. Coordinate Systems in Geodesy.
           Jan. 1998, pp. 31–33.

        Examples
        --------
        >>> timeData = Time(julian=np.array([2454545]))
        >>> basePos = np.array([[6343.82, -2640.87, -11.26]])
        >>> coor = Coordinate(basePos=basePos, type="ECEF", timeData=timeData)
        >>> coor.altitude()
        np.array([504.1269764])
        """

        if self._ECEF is None:
            self.ECEF()

        ECEF_data = self._ECEF
        x = ECEF_data[:, 0]
        y = ECEF_data[:, 1]
        z = ECEF_data[:, 2]

        a = 6378.137
        b = 6356.752314245

        epsilon = 10e-10

        e = np.sqrt(1 - b ** 2 / a ** 2)
        p = np.sqrt(x ** 2 + y ** 2)

        N = a
        h = np.sqrt(x ** 2 + y ** 2 + z ** 2) - np.sqrt(a * b)
        phi = np.arctan((z / p) * (1 - (e ** 2 * N) / (N + h)) ** -1)

        def new_vals(a, b, e, p, N, h, phi, z):

            N = a / np.sqrt(np.cos(phi) ** 2 + b ** 2 / a ** 2 * np.sin(phi) ** 2)
            h = p / np.cos(phi) - N
            phi = np.arctan((z / p) * (1 - (e ** 2 * N) / (N + h)) ** -1)

            return N, h, phi

        Np, hp, phip = new_vals(a, b, e, p, N, h, phi, z)

        while (np.mean(hp - h) > a * epsilon) and (np.mean(phip - phi) > epsilon):

            N, h, phi = Np, hp, phip
            Np, hp, phip = new_vals(a, b, e, p, N, h, phi, z)

        if kwargs:
            h = self._interp(h, **kwargs)

        return h

    def distance(self, groundPos: Any, **kwargs: Dict[str, Any]) -> np.array:
        """Return the distance to a ground location.

        Parameters
        ----------
        groundPos : GroundPosition
            Calculate distances to the specified ground location.
        kwargs : Dict, optional
            Optional keyword arguments passed into the
            `Interpolation()._interp()` method.

        Returns
        -------
        np.array
            Array of shape (n,) containing distances in kilometres between the
            satellite and ground location for each satellite position.
        """

        if self._ECEF is None:
            self.ECEF()

        coor = groundPos.coor
        gnd_ECEF = self._GEO_to_ECEF(np.array([[coor[0], coor[1], 0]]))
        gnd_ECEF = np.repeat(gnd_ECEF, self.length, 0)

        # Find LOS vector norm.
        LOS = self._ECEF - gnd_ECEF
        distances = np.linalg.norm(LOS, axis=1)

        if kwargs:
            distances = self._interp(distances, **kwargs)

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

        Examples
        --------
        >>> angles = np.array([43.6532, -79.3832, -33.2833, 149.1000])
        >>> coor = Coordinate(...)
        >>> coor.sexagesimal(angles=angles)
        np.array(['+43°39′11.52″',
                  '-79°22′59.52″',
                  '-33°16′59.88″',
                  '+149°06′00.00″'])
        """

        deg = u"\u00B0"
        min = u"\u2032"
        sec = u"\u2033"

        length = angles.shape[0]
        out_arr = np.empty((length,), dtype="<U32")

        for i in range(length):

            ang = angles[i]

            sign = "+" if ang >= 0 else "-"
            degree = int(abs(ang))
            minute = int(60 * np.round(abs(ang) - degree, 2))
            second = np.round(abs(60 * (60 * (abs(ang) - degree) - minute)), 2)

            degree = "%.3s" % str(degree).zfill(2)
            minute = "%.2s" % str(minute).zfill(2)
            second = str("%.2f" % second).zfill(5)

            out_arr[i] = f"{sign}{degree}{deg}{minute}{min}{second}{sec}"

        return out_arr
