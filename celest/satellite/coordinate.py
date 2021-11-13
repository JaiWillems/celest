"""Localize position information representations.

This module provides the `Coordinate` class to allow for a simple interface
for converting a position type into various representations. The class is also
used for position inputs and outputs for other `Celest` functionality.
"""


from celest.core.decorators import set_module
from celest.satellite._angle_representations import _ISO6709_representation
from typing import Any, Literal, Tuple
import numpy as np


@set_module('celest.satellite')
class Coordinate(object):
    """Localize position information representations.

    The `Coordinate` class provides a simple user interface for converting a
    position type into various representations. The class is also used for
    position inputs and outputs for other `Celest` functionality.

    Parameters
    ----------
    position : np.ndarray
        Base position to initialize the `Coodinate` class.
    frame : {"gcrs", "geo", "itrs"}
        Specifies the input position frame.
    time : Time
        Times associated with the position data in the J2000 epoch. The length
        of the `time` parameter must match the length of the `position`
        parameter.

    Attributes
    ----------
    time : Time
        Times corresponding to the position data.
    length : int
        Length of the input position and time arrays.

    Methods
    -------
    geo(iso=False)
        Return geographical position data.
    era()
        Return the Earth rotation angles in radians and decimals.
    gcrs()
        Return cartesian gcrs position data.
    itrs()
        Return cartesian itrs position data.
    horizontal(location)
        Return horizontal position data in degrees and decimals.
    off_nadir(location)
        Return the off-nadir angle to a ground location.
    altitude()
        Return the altitude above Earth's surface in kilometres.
    distance(location)
        Return the distance to a ground location.
    """

    def __init__(self, position: np.ndarray, frame: Literal["gcrs", "geo", "itrs"], time: Any) -> None:
        """Initialize attributes."""
        
        if frame not in ["gcrs", "geo", "itrs"]:
            raise ValueError(f"{frame} is not a valid frame.")
        
        if position.shape[0] != len(time):
            raise ValueError(f"position and time data lengths are mismatched being {position.shape[0]} and {len(time)}")

        self.time = time

        self._GEO = None
        self._GCRS = None
        self._ITRS = None

        self.length = None
        self._set_base_position(position, frame)
    
    def __len__(self) -> int:
        """Return length of position and time data."""

        return self.length

    def _set_base_position(self, position: np.ndarray, frame:
                           Literal["gcrs", "geo", "itrs"]) -> None:
        """Initialize base position.

        This method takes an input position to initialize the object's base
        position.

        Parameters
        ----------
        position : np.ndarray
            Array of shape (n, 2) or (n, 3) containing the inputted position
            data.
        frame : {"gcrs", "geo", "itrs"}
            String defining the type of input position data as the Geocentric
            Celestial Reference System (gcrs), International Terrestrial
            Reference System (itrs), or geographical (geo) data.

        Notes
        -----
        The input data must be of shape (n, 3) if `type="itrs"` or `type="gcrs"`
        where the columns are XYZ cartesian data. The data can be of the shape
        (n, 2) or (n, 3) if `type="geo"` where the columns are geodetic
        latitude, terrestrial longitude, and geodetic altitude. When
        geographical data is entered of shape (n, 2), the height data is
        assumed to be zero.
        """
        
        basePos = position
        self.length = basePos.shape[0]

        if frame == "geo":
            if basePos.shape[1] == 2:
                height = np.zeros((self.length, 1))
                basePos = np.concatenate((basePos, height), axis=1)

            self._GEO = basePos

        elif frame == "gcrs":
            self._GCRS = basePos

        elif frame == "itrs":
            self._ITRS = basePos

    def _geo_to_itrs(self, position: np.ndarray) -> np.ndarray:
        """Convert geographical to itrs coordinates.

        Parameters
        ----------
        position : np.ndarray
            Array of shape (n, 3) containing geographical coordinates of a
            position with columns of geodetic latitude, terrestrial longitude,
            and geodetic altitude given in decimal degrees and kilometres.

        Returns
        -------
        np.ndarray
            Array of shape (n, 3) containing the XYZ itrs position data.

        See Also
        --------
        _itrs_to_geo : Convert itrs to geographical coordinates.

        Notes
        -----
        This method uses an ellipsoid based model of the Earth to convert a
        geographical position to itrs cartesian coordinates using the methods
        described in "Coordinate Systems in Geodesy" by E. J. Krakiwsky and
        D.E. Wells as presented by Christopher Lum. [KW98a]_ [Lum20]_

        References
        ----------
        .. [KW98a] E. J. Krakiwsky and D. E. Wells. Coordinate Systems in Geodesy.
           Jan. 1998.
        .. [Lum20] Christopher Lum. Geodetic Coordinates: Computing Latitude and
           Longitude.June 2020.url:https://www.youtube.com/watch?v=4BJ-GpYbZlU.
        """

        a = 6378.137
        b = 6356.752314245

        lat, lon = np.radians(position[:, 0]), np.radians(position[:, 1])

        e = np.sqrt(1 - b ** 2 / a ** 2)
        N = a / np.sqrt(1 - e ** 2 * np.sin(lat) ** 2)
        h = position[:, 2]

        x = ((N + h) * np.cos(lat) * np.cos(lon)).reshape((-1, 1))
        y = ((N + h) * np.cos(lat) * np.sin(lon)).reshape((-1, 1))
        z = ((N * (1 - e ** 2) + h) * np.sin(lat)).reshape((-1, 1))
        itrs = np.concatenate((x, y, z), axis=1)

        return itrs

    def _itrs_to_geo(self, position: np.ndarray) -> np.ndarray:
        """Convert itrs to geographical coordinates.

        Parameters
        ----------
        position : np.ndarray
            Array of shape (n, 3) with rows containing XYZ itrs position data.

        Returns
        -------
        np.ndarray
            Array of shape (n, 3) with columns of geodetic latitude,
            terrestrial longitude, and geodetic altitude data in degrees and
            kilometres.

        See Also
        --------
        _GEO_to_itrs : Convert geographical to itrs coordinates.

        Notes
        -----
        Let :math:`x`, :math:`y`, and :math:`z` be the itrs vector components.
        We can then calculate the latitude, :math:`\phi`, using the following
        equation:

        .. math:: \phi = 90^\circ - \cos^{-1}\left(\frac{z}{R_\bigoplus}\right)

        where :math:`R_\bigoplus` is the geocentric radius of the Earth. We can
        also calculate the longitude, :math:`\lambda` using the following:

        .. math::\lambda = \tan^{-1}_2\left(y, x\right)
        """

        # Cartesian coordinates.
        x = position[:, 0]
        y = position[:, 1]
        z = position[:, 2]

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

    def geo(self, iso: bool=False) -> np.ndarray:
        """Return geographical position data.

        Parameters
        ----------
        iso : bool, optional
            Formats the output as ISO6709 sexagesimal position strings if true.

        Returns
        -------
        np.ndarray
            Array of shape (n, 3) with columns of geodetic latitude,
            terrestrial longitude, and geodetic altitude data. If
            `iso=True`, the output is an array of shape (n,) containing
            standard representation position strings.

        See Also
        --------
        itrs : Return itrs position data.
        gcrs : Return gcrs position data.

        Examples
        --------
        >>> time = Time(julian=np.array([2454545]))
        >>> position = np.array([[6343.82, -2640.87, -11.26]])
        >>> coor = Coordinate(position=position, frame="itrs", time=time)
        >>> coor.geo()
        np.array([[-9.38870528e-02, -2.26014826e+01, 5.04126976e+02]])

        We can generate user-friendly position strings by setting `iso=True`.

        >>> coor.geo(iso=True)
        np.array(['00°05′37.99″S 22°36′05.34″W 504.12697'])
        """

        if self._GEO is None:
            if self._ITRS is not None:
                self._GEO = self._itrs_to_geo(position=self._ITRS)
            else:
                self._ITRS = self._gcrs_and_itrs(position=self._GCRS, frame="gcrs")
                self._GEO = self._itrs_to_geo(position=self._ITRS)

        geo = _ISO6709_representation(position=self._GEO) if iso else self._GEO

        return geo

    def era(self) -> np.ndarray:
        """Return the Earth rotation angles in degrees and decimals.

        Returns
        -------
        np.ndarray
            Array of shape (n,) containing Earth rotation angles in degrees and
            decimals.

        Notes
        -----
        The Earth rotation angle denotes the diurnal rotation component of the
        coordinate misalignment between the ECI and ECEF frames. The angle can
        be calculated from the following formulation:

        .. math:: \gamma^\circ = 360.9856123035484\Delta T + 280.46

        where :math:`\Delta T=JD-2451545` is the elapsed days since the J2000
        epoch where :math:`JD` is the Julian day. [Kok17a]_

        References
        ----------
        .. [Kok17a] Don Koks. Changing Coordinates in the Context of Orbital
           Mechanics. Cyber and Electronic Warfare Division, Defence Science,
           and Technology Group, Jan.2017, p. 12 - 13.

        Examples
        --------
        >>> time = Time(julian=np.array([2454545]))
        >>> position = np.array([[6343.82, -2640.87, -11.26]])
        >>> coor = Coordinate(position=position, type="itrs", time=time)
        >>> coor.era()
        np.array([6.2360075])
        """

        ang = np.zeros((self.length,))
        jul_data = self.time.julian()

        # Multiply time elapsed since J2000 by Earth rotation rate and add
        # J2000 orientation.
        dJulian = jul_data - 2451545
        ang = (360.9856123035484 * dJulian + 280.46) % 360

        return ang

    def _gcrs_and_itrs(self, position: np.ndarray, frame: Literal["itrs", "gcrs"]) -> np.ndarray:
        """Convert between gcrs and itrs positions.

        Parameters
        ----------
        position : np.ndarray
            Array of shape (n, 3) representing the input data as XYZ cartesian
            data.
        frame : {"itrs", "gcrs"}
            The type of input data.

        Returns
        -------
        np.ndarray
            Array of shape (n, 3) representing the output data as XYZ cartesian
            data.

        Notes
        -----
        This method applies a simplification to converting between gcrs and
        itrs coordinate by disregarding precession, nutation, and polar motion
        effects (which will be incorporated in future releases). By
        disregarding such factors, we can assume the alignment of the gcrs and
        itrs axes and simplify the conversion to a rotation about the z-axis.

        We can convert between the gcrs and itrs frames using the following
        equations:

        .. math:: \vec{V}_{gcrs} = E_3^\gamma\;\vec{V}_{itrs}

        .. math:: \vec{V}_{itrs} = E_3^{-\gamma}\;\vec{V}_{gcrs}

        where :math:`\gamma` is the Earth rotation angle, :math:`E_3^\theta` is
        the rotation matrix that rotates a vector around the z-axis by an angle
        :math:`\theta`, and :math:`\vec{V}` is a position vector in either the
        ECI or ECEF frames. [Kok17b]_

        References
        ----------
        .. [Kok17b] Don Koks. Changing Coordinates in the Context of Orbital
           Mechanics. Cyber and Electronic Warfare Division, Defence Science,
           and Technology Group, Jan.2017, p. 21.
        """

        theta = -self.era() if frame == "gcrs" else self.era()
        theta = np.radians(theta)

        # Construct rotational matrix.
        A11 = np.cos(theta)
        A12 = -np.sin(theta)
        A21 = np.sin(theta)
        A22 = np.cos(theta)

        # Rotate position data around z-axis by ERA.
        output = np.zeros((self.length, 3))
        output[:, 0] = A11 * position[:, 0] + A12 * position[:, 1]
        output[:, 1] = A21 * position[:, 0] + A22 * position[:, 1]
        output[:, 2] = position[:, 2]

        return output

    def gcrs(self) -> np.ndarray:
        """Return cartesian gcrs position data.

        Returns
        -------
        np.ndarray
            Array of shape (n,3) with columns of XYZ gcrs position data.

        See Also
        --------
        itrs : Return cartesian itrs position data.
        geo : Return geographical position data.

        Examples
        --------
        >>> time = Time(julian=np.array([2454545]))
        >>> position = np.array([[6343.82, -2640.87, -11.26]])
        >>> coor = Coordinate(position=position, type="itrs", time=time)
        >>> coor.gcrs()
        np.array([[6212.21719598, -2937.10811161, -11.26]])
        """

        if self._GCRS is None:
            if self._ITRS is None:
                self._ITRS = self._geo_to_itrs(self._GEO)
            self._GCRS = self._gcrs_and_itrs(self._ITRS, frame="itrs")

        return self._GCRS

    def itrs(self) -> np.ndarray:
        """Return cartesian itrs position data.

        Returns
        -------
        np.ndarray
            Array of shape (n,3) with columns of XYZ itrs position data.

        See Also
        --------
        gcrs : Return cartesian gcrs position data.
        geo : Return geographical position data.

        Examples
        --------
        >>> time = Time(julian=np.array([2454545]))
        >>> position = np.array([[6343.82, -2640.87, -11.26]])
        >>> coor = Coordinate(position=position, type="gcrs", time=time)
        >>> coor.itrs()
        np.array([[6461.30569276, -2338.75507354, -11.26]])
        """

        if self._ITRS is None:
            if self._GCRS is not None:
                self._ITRS = self._gcrs_and_itrs(position=self._GCRS, frame="gcrs")
            else:
                self._ITRS = self._geo_to_itrs(position=self._GEO)

        return self._ITRS

    def _get_ang(self, u: np.ndarray, v: np.ndarray) -> np.ndarray:
        """Calculate degree angle bewteen two vectors.

        Parameters
        ----------
        u, v : np.ndarray
            Arrays of shape (n,3) with rows of XYZ cartesian data.

        Returns
        -------
        np.ndarray
            Array of shape (n,) containing the degree angle between the two
            arrays.
        """

        num = np.sum(u * v, axis=1)
        denom = np.linalg.norm(u, axis=1) * np.linalg.norm(v, axis=1)
        ang = np.degrees(np.arccos(num / denom))

        return ang

    def horizontal(self, location: Any) -> Tuple:
        """Return horizontal position data in degrees and decimals.

        In this method, the azimuth is calculated as increasing degrees
        clockwise from North and ranges from 0 to 360 degrees. The altitude is
        calculated as increasing degrees above the observer's local horizon and
        ranges from 0 to 90 degrees.

        Parameters
        ----------
        location : GroundPosition
            Ground location defining the center of the horizontal system.

        Returns
        -------
        tuple
            Tuple of length two containing the altitude and azimuth data as
            decimals and degrees in NumPy arrays of shape (n,).

        Examples
        --------
        >>> time = Time(julian=np.array([2454545]))
        >>> position = np.array([[6343.82, -2640.87, -11.26]])
        >>> location = GroundPosition(52.1579, -106.6702)
        >>> coor = Coordinate(position=position, frame="ecef", time=time)
        >>> coor.horizontal(location=location)
        (array([-40.8786098]), array([94.73615482]))
        """

        if self._ITRS is None:
            self.itrs()

        # Convert observer position into cartesian coordinates.
        lat, lon, radius = location.lat, location.lon, location.radius
        GEO_data = self._geo_to_itrs(np.array([[lat, lon, 0]]))
        obs = np.repeat(GEO_data, self.length, 0)

        # Determine line of sight vector then altitude.
        LOS = self._ITRS - obs
        Alt = 90 - self._get_ang(LOS, obs)

        # Find surface tangent vector passing through z-axis.
        k_hat = np.repeat(np.array([[0, 0, 1]]), self.length, axis=0)
        beta = np.radians(lat)
        tangent = (k_hat.T * radius / np.sin(beta)).T - obs

        # Find LOS projection on tangent plane.
        norm_proj = (obs.T * np.sum(LOS * obs, axis=1) / radius ** 2).T
        proj_LOS = LOS - norm_proj

        # Determing azimuth.
        reference = np.cross(tangent, obs)
        neg_ind = np.where(np.sum(proj_LOS * reference, axis=1) < 0)[0]
        Az = self._get_ang(tangent, proj_LOS)
        Az[neg_ind] = 360 - self._get_ang(tangent[neg_ind], proj_LOS[neg_ind])

        return Alt, Az

    def off_nadir(self, location: Any) -> np.ndarray:
        """Return the off-nadir angle to a ground location.

        This method calculates the off-nadir angle, in degrees and decimals, to
        the input ground position at each satellite position.

        Parameters
        ----------
        location : GroundPosition
            Calculate off-nadir angles for the specified ground location.

        Returns
        -------
        np.ndarray
            Array of shape (n,) containing off-nadir angle data in degrees and
            decimals.

        Notes
        -----
        The off-nadir angle is the acute angle measured in increasing degrees
        from the satellites nadir to the line joining the satellite to the
        ground location of interest.
        """

        if self._ITRS is None:
            self.itrs()

        lat, lon = location.lat, location.lon
        geo_data = self._geo_to_itrs(np.array([[lat, lon, 0]]))
        obs = np.repeat(geo_data, self.length, 0)

        LOS = np.subtract(self._ITRS, obs)
        ang = self._get_ang(LOS, self._ITRS)

        return ang

    def _WGS84_radius(self, lattitude: np.ndarray) -> np.ndarray:
        """Calculate the Earth's geocentric radius using WGS84.

        Parameters
        ----------
        latitude : np.ndarray
            Array of shape (n,) representing the geocentric latitude of a
            ground location in degrees and decimals.

        Returns
        -------
        np.ndarray
            Array of shape (n,) containing the Earth's geocentric radius in
            kilometres.

        Notes
        -----
        By using an Earth ellipsoid with the WGS84 parameters of
        :math:`a=6378.137` and :math:`b=6356.7523142`, the geocentric radius
        can be calculated using the following formulation:

        .. math:: r = \sqrt{\frac{(a^2\cos(\beta))^2 + (b^2\sin(\beta))^2}{(a\cos(\beta))^2 + (b\sin(\beta))^2}}

        where :math:`\beta` is the observer's latitude. [Tim18]_

        References
        ----------
        .. [Tim18] Timur. Earth Radius by Latitude (WGS 84). 2018. url:
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

    def altitude(self) -> np.ndarray:
        """Return the geodetic altitude above Earth's surface in kilometres.

        This method uses the WGS84 reference ellipsoid to calculate the
        geodetic altitude above the Earth's surface.

        Returns
        -------
        np.ndarray
            Array of shape (n,) containing satellite altitudes in kilometres.

        Notes
        -----
        This method uses an ellipsoid based model of the Earth to calculate the
        ellipsoid height in an iterative manner described in "Coordinate
        Systems in Geodesy" by E. J. Krakiwsky and D.E. Wells. [KW98b]_

        References
        ----------
        .. [KW98b] E. J. Krakiwsky and D. E. Wells. Coordinate Systems in Geodesy.
           Jan. 1998, pp. 31–33.

        Examples
        --------
        >>> time = Time(julian=np.array([2454545]))
        >>> position = np.array([[6343.82, -2640.87, -11.26]])
        >>> coor = Coordinate(position=position, type="itrs", time=time)
        >>> coor.altitude()
        np.array([504.1269764])
        """

        if self._ITRS is None:
            self.itrs()

        itrs_data = self._ITRS
        x = itrs_data[:, 0]
        y = itrs_data[:, 1]
        z = itrs_data[:, 2]

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

        return h

    def distance(self, location: Any) -> np.ndarray:
        """Return the distance to a ground location.

        Parameters
        ----------
        location : GroundPosition
            Calculate distances to the specified ground location.

        Returns
        -------
        np.ndarray
            Array of shape (n,) containing distances in kilometres between the
            satellite and ground location for each satellite position.
        """

        if self._ITRS is None:
            self.itrs()

        lat, lon = location.lat, location.lon
        gnd_itrs = self._geo_to_itrs(np.array([[lat, lon, 0]]))
        gnd_itrs = np.repeat(gnd_itrs, self.length, 0)

        # Find LOS vector norm.
        LOS = self._ITRS - gnd_itrs
        distances = np.linalg.norm(LOS, axis=1)

        return distances
