

from celest.satellite._angle_representations import _ISO6709_representation
from celest.satellite.nutation_precession import (
    bias_matrix, precession_matrix, nutation_matrix
)
from celest.satellite.time import Time
from polare import Stroke
from typing import Tuple
import numpy as np


class Coordinate(Time):
    """Coordinate(position, frame, julian, offset=0)

    Satellite coordinate transformations.

    The input position and time can be converted into different coordinate
    representations useful for satellite applications. Here, `julian + offset`
    is the Julian time in J2000 epoch and is used in various transformations.

    If `frame=="geo"`, the `position` input can have two columns, (latitude,
    longitude) or three columns (latitude, longitude, altitude). If no
    geodetic altitude is provided, it is assumed to be zero. Other frames
    require three columns.

    Parameters
    ----------
    position : array_like
        2-D array containing position coordinates. Axis 0 is the time series
        axis.

        Supported coordinate frames include the Geocentric Celestial Reference
        System (gcrs), International Terrestrial Reference System (itrs), and
        geographical (geo) system.
    frame : {"gcrs", "geo", "itrs"}
        Frame specifier for the input position.
    julian : array_like
        1-D array containing time data in Julian days.

        If times are not in the J2000 epoch, a non-zero offset must be pased
        in to add to the julian times.
    offset : float, optional
        Offset to convert input time data to the J2000 epoch, default is zero.

    Examples
    --------
    Initialize `Coordinate` using gcrs positions:

    >>> julian = [30462.50, 30462.50]
    >>> position = [[-4681.50824149 -5030.09119386     0.        ]
    ...             [-4714.35351825 -4978.74325953   454.41492765]]
    >>> c = Coordinate(position=position, frame="gcrs", julian=julian, offset=0)

    Get horizontal coordinates for the ground location at (43.65, -79.38):

    >>> location = GroundPosition(latitude=43.65, longitude=-79.38)
    >>> altitude, azimuth = c.horizontal(location=location)
    """

    def __init__(self, position, frame, julian, offset=0) -> None:

        super().__init__(julian, offset)

        time, position = np.array(julian), np.array(position)

        if frame not in ["gcrs", "geo", "itrs"]:
            raise ValueError(f"{frame} is not a valid frame.")

        if position.ndim != 2:
            raise ValueError("position input must be 2-D.")

        if time.ndim != 1:
            raise ValueError("time input must be 1-D.")

        if position.shape[0] != time.size:
            raise ValueError(f"position and time data lengths are mismatched: "
                             f"lengths are {position.shape[0]} and {len(time)}")

        kind_dict = {1: "zero", 2: "linear", 3: "quadratic"}
        if position.shape[0] > 3:
            self._kind = "cubic"
        else:
            self._kind = kind_dict[position.shape[0]]

        self._GCRS = None
        self._ITRS = None
        self._length = None

        self._set_base_position(position, frame)

    def __len__(self) -> int:

        return self._length

    def _set_base_position(self, position, frame) -> None:
        """Initialize base position.

        Parameters
        ----------
        position : array_like
            2-D array containing position coordinates.

            Supported coordinate frames include the Geocentric Celestial
            Reference System (gcrs), International Terrestrial Reference
            System (itrs), and geographical (geo) system.
        frame : {"gcrs", "geo", "itrs"}
            Frame specifier for the input position.

        Notes
        -----
        GCRS and ITRS positions shall have 3 columns representing XYZ
        cartesian points. GEO positions can have either 2 or 3 columns where
        the first two represent geodetic latitude and terrestrial longitude;
        the last column represent geodetic altitude and is assumed zero if not
        passed in.
        """

        self._length = position.shape[0]

        if position.shape[1] == 2:
            height = np.zeros((self._length, 1))
            position = np.concatenate((position, height), axis=1)

        if frame == "gcrs":
            x, y, z = [Stroke(self._julian, position[:, i], self._kind) for i in range(3)]
            self._GCRS = np.array([x, y, z])
        elif frame == "itrs":
            x, y, z = [Stroke(self._julian, position[:, i], self._kind) for i in range(3)]
            self._ITRS = np.array([x, y, z])
        else:
            itrs = self._geo_to_itrs(position)
            x, y, z = [Stroke(self._julian, itrs[:, i], self._kind) for i in range(3)]
            self._ITRS = np.array([x, y, z])

    def _geo_to_itrs(self, position) -> np.ndarray:
        """Geographical to itrs transformation.

        Parameters
        ----------
        position : np.ndarray
            2-D array with columns of geodetic latitude, terrestrial longitude,
            and geodetic altitude given in decimal degrees and kilometres.

        Returns
        -------
        np.ndarray
            2-D array with columns of x, y, z cartesian coordinates in ITRS.

        See Also
        --------
        _itrs_to_geo : Itrs to geographical transformation.

        Notes
        -----
        An Earth ellipsoid model is used for the geographical to itrs
        conversion using the methods described in "Coordinate Systems in
        Geodesy" by E. J. Krakiwsky and D.E. Wells as presented by Christopher
        Lum. [KW98a]_ [Lum20]_

        References
        ----------
        .. [KW98a] E. J. Krakiwsky and D. E. Wells. Coordinate Systems in
           Geodesy. Jan. 1998.
        .. [Lum20] Christopher Lum. Geodetic Coordinates: Computing Latitude and
           Longitude. June 2020.url:https://www.youtube.com/watch?v=4BJ-GpYbZlU.
        """

        a = 6378.137  # WGS84 major axis in km.
        b = 6356.752314245  # WGS84 minor axis in km.

        lat, lon = np.copy(np.radians(position[:, 0:2])).T
        clat, slat = np.cos(lat), np.sin(lat)
        clon, slon = np.cos(lon), np.sin(lon)

        e = np.sqrt(1 - b ** 2 / a ** 2)
        n = a / np.sqrt(1 - e ** 2 * slat ** 2)
        h = position[:, 2]

        x = ((n + h) * clat * clon).reshape((-1, 1))
        y = ((n + h) * clat * slon).reshape((-1, 1))
        z = ((n * (1 - e ** 2) + h) * slat).reshape((-1, 1))

        return np.concatenate((x, y, z), axis=1)

    def _itrs_to_geo(self, position) -> np.ndarray:
        """Itrs to geographical transformation.

        Parameters
        ----------
        position : np.ndarray
            2-D array with columns of x, y, z cartesian coordinates in ITRS.

        Returns
        -------
        np.ndarray
            2-D array with columns of latitude, longitude, and altitude given
            in decimal degrees and kilometres.

        See Also
        --------
        _geo_to_itrs : Geographical to itrs transformation.
        """

        # Calculate latitude and longitude.
        x, y, z = position
        lat = np.degrees(np.arctan(z / np.sqrt(x ** 2 + y ** 2)))
        lon = np.degrees(np.arctan2(y, x))

        return lat, lon

    def geo(self, iso=False, stroke=False) -> np.ndarray:
        """Return geographical coordinates.

        The method can not return the data as both a stroke and iso
        formatted strings. If both `iso=True` and `stroke=True`, the method
        will return a tuple of stroke objects.

        Parameters
        ----------
        iso : bool, optional
            Formats output as ISO6709 sexagesimal position strings if true.
        stroke : bool, optional
            Formats output as a tuple of stroke objects if true. Otherwise,
            the output is a tuple of 1-D arrays.

        Returns
        -------
        tuple
            Tuple of latitude, longitude, and altitude given in decimal
            degrees.

            If `iso=True`, the output is a 1-D array containing formatted
            strings.

            If `stroke=True`, the output is a tuple of stroke objects.

        See Also
        --------
        itrs : Return itrs coordinates.
        gcrs : Return gcrs coordinates.

        Examples
        --------
        Initialize object and get geographical data:

        >>> julian=[2454545]
        >>> position = [[6343.82, -2640.87, -11.26]]
        >>> c = Coordinate(position=position, frame="itrs", julian=julian)
        >>> c.geo()
        np.array([[-9.38870528e-02, -2.26014826e+01, 5.04126976e+02]])

        Return formatted output strings by setting `iso=True`:

        >>> c.geo(iso=True)
        np.array(['00°05′37.99″S 22°36′05.34″W 504.12697'])
        """

        if self._ITRS is None:
            x, y, z = self._GCRS
            x, y, z = [i(self._julian).reshape((-1, 1)) for i in [x, y, z]]
            gcrs = np.concatenate((x, y, z), axis=1)

            itrs = self._gcrs_and_itrs(gcrs, "gcrs")
            x, y, z = [Stroke(self._julian, itrs[:, i], self._kind) for i in range(3)]
            self._ITRS = np.array([x, y, z])

        lat, lon = self._itrs_to_geo(self._ITRS)
        alt = self.altitude(stroke=True)

        if stroke:
            return lat, lon, alt

        lat, lon, alt = lat(self._julian), lon(self._julian), alt(self._julian)

        if iso:
            return _ISO6709_representation(lat, lon, alt)
        else:
            return lat, lon, alt

    def era(self, stroke=False) -> np.ndarray:
        """Return Earth rotation angle in decimal degrees.

        Parameters
        ----------
        stroke : bool, optional
            Formats output as a stroke objects if true. Otherwise, the output
            is a 1-D array.

        Returns
        -------
        Union[np.ndarray, Stroke]
            Earth rotation angle in decimal degrees.

            If `stroke=True`, the output is a stroke object. Otherwise, the
            output is a 1-D array.

        Notes
        -----
        The Earth rotation angle is calculated using the following:

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
        >>> julian = [2454545]
        >>> position = [[6343.82, -2640.87, -11.26]]
        >>> c = Coordinate(position=position, type="itrs", julian=julian)
        >>> c.era()
        np.array([6.2360075])
        """

        # Multiply time elapsed since J2000 by Earth rotation rate and add
        # J2000 orientation.
        dJulian = self._julian - 2451545
        ang = (360.9856123035484 * dJulian + 280.46) % 360

        if stroke:
            return Stroke(self._julian, ang, self._kind)
        else:
            return ang

    def _gcrs_and_itrs(self, position, frame) -> np.ndarray:
        """Transform between gcrs and itrs coordinates.

        Parameters
        ----------
        position : np.ndarray
            2-D array containing x, y, z cartesian coordinates.
        frame : {"itrs", "gcrs"}
            Input data frame.

        Returns
        -------
        np.ndarray
            2-D array containing the transformed x, y, z cartesian data.

        Notes
        -----
        The conversion between the gcrs and itrs frames accounts for Earth
        rotation, precession, and nutation. Polar motion effects are ignored
        due to their poor predictive nature.
        """

        pos_data = np.copy(position)

        # Get dependencies.
        era = self.era()
        bias_mtrx = bias_matrix()
        precession_mtrx = precession_matrix(self._julian)
        nutation_mtrx = nutation_matrix(self._julian)

        # Apply nutation-precession model before ERA rotation.
        if frame == "gcrs":
            pos_data = np.einsum('ij, kj -> ki', bias_mtrx, pos_data)
            pos_data = np.einsum('ijk, ik -> ij', precession_mtrx, pos_data)
            pos_data = np.einsum('ijk, ik -> ij', nutation_mtrx, pos_data)

        # Pre-process ERA angle and sign.
        theta = np.radians(era)
        theta = -theta if frame == "gcrs" else theta

        # Rotate position by the ERA rotational matrix.
        c, s = np.cos(theta), np.sin(theta)
        c1, c2, c3 = np.copy(pos_data).T

        pos_data[:, 0] = c * c1 - s * c2
        pos_data[:, 1] = s * c1 + c * c2
        pos_data[:, 2] = c3

        # Apply nutation-precession model after ERA rotation if frame is ITRS.
        if frame == "itrs":
            pos_data = np.einsum('ijk, ik -> ij',
                                 np.linalg.inv(nutation_mtrx), pos_data)
            pos_data = np.einsum('ijk, ik -> ij',
                                 np.linalg.inv(precession_mtrx), pos_data)
            pos_data = np.einsum('ij, kj -> ki',
                                 np.linalg.inv(bias_mtrx), pos_data)

        return pos_data

    def gcrs(self, stroke=False) -> np.ndarray:
        """Return gcrs coordinates.

        Parameters
        ----------
        stroke : bool, optional
            Formats output as a tuple of stroke objects if true. Otherwise,
            the output is a tuple of 1-D arrays.

        Returns
        -------
        tuple
            Tuple containing x, y, z cartesian coordinates in the gcrs frame.

            The coordinates are stroke objects if `stroke=True`. Otherwise,
            1-D arrays are returned.

        See Also
        --------
        itrs : Return itrs coordinates.
        geo : Return geographical coordinates.

        Examples
        --------
        Calculate gcrs coordinates from itrs coordinates:

        >>> juluian = [2454545]
        >>> position = [[6343.82, -2640.87, -11.26]]
        >>> c = Coordinate(position=position, type="itrs", julian=julian)
        >>> c.gcrs()
        np.array([[6212.21719598, -2937.10811161, -11.26]])
        """

        if self._GCRS is None:
            itrs = np.array([self._ITRS[i](self._julian) for i in range(3)]).T
            gcrs = self._gcrs_and_itrs(itrs, "itrs")
            x, y, z = [Stroke(self._julian, gcrs[:, i], self._kind) for i in range(3)]
            self._GCRS = np.array([x, y, z])

        x, y, z = self._GCRS[0], self._GCRS[1], self._GCRS[2]
        if stroke:
            return x, y, z
        else:
            return x(self._julian), y(self._julian), z(self._julian)

    def itrs(self, stroke=False) -> np.ndarray:
        """Return itrs coordinates.

        Parameters
        ----------
        stroke : bool, optional
            Formats output as a tuple of stroke objects if true. Otherwise,
            the output is a tuple of 1-D arrays.

        Returns
        -------
        tuple
            Tuple containing x, y, z cartesian data in the itrs frame.

            The coordinates are stroke objects if `stroke=True`. Otherwise,
            1-D arrays are returned.

        See Also
        --------
        gcrs : Return gcrs coordinates.
        geo : Return geographical coordinates.

        Examples
        --------
        Calculate itrs coordinates from gcrs coordinates:

        >>> julian = [2454545]
        >>> position = [[6343.82, -2640.87, -11.26]]
        >>> c = Coordinate(position=position, type="gcrs", julian=julian)
        >>> c.itrs()
        np.array([[6461.30569276, -2338.75507354, -11.26]])
        """

        if self._ITRS is None:
            gcrs = np.array([self._GCRS[i](self._julian) for i in range(3)]).T
            itrs = self._gcrs_and_itrs(gcrs, "gcrs")
            x, y, z = [Stroke(self._julian, itrs[:, i], self._kind) for i in range(3)]
            self._ITRS = np.array([x, y, z])

        x, y, z = self._ITRS[0], self._ITRS[1], self._ITRS[2]
        if stroke:
            return x, y, z
        else:
            return x(self._julian), y(self._julian), z(self._julian)

    def _get_ang(self, u, v) -> np.ndarray:
        """Return the degree angle bewteen two vectors.

        Parameters
        ----------
        u, v : np.ndarray
            1-D or 2-D arrays containing row vectors.

            If the dimensions of `u` and `v` do not match, the 1-D array will
            be broadcast to have the same number of rows as the 2-D array.

        Returns
        -------
        np.ndarray
            1-D array containing the degree angle between to vector arrays.
        """

        ua = None if u.ndim == 1 else 1
        va = None if v.ndim == 1 else 1

        n = np.sum(u * v, axis=(ua or va))
        d = np.linalg.norm(u, axis=ua) * np.linalg.norm(v, axis=va)
        ang = np.degrees(np.arccos(n / d))

        return ang

    def _altitude(self, location) -> Stroke:
        """Return the altitude angle of a satellite.

        Parameters
        ----------
        location : GroundPosition
            Eath bound origin of the horizontal system.

        Returns
        -------
        Stroke
            Stroke containing the altitude of the satellite.
        """

        # Get origin of horizontal system.
        lat, lon, radius = location.lat, location.lon, location.radius
        gnd_itrs = self._geo_to_itrs(np.array([[lat, lon, 0]])).reshape((3,))

        return 90 - self._get_ang(self._ITRS - gnd_itrs, gnd_itrs)

    def _azimuth(self, location) -> Stroke:
        """Return the azimuth angle of the satellite.

        Parameters
        ----------
        location : GroundPosition
            Eath bound origin of the horizontal system.

        Returns
        -------
        Stroke
            Stroke containing the azimuth of the satellite.
        """

        # Get origin of horizontal system.
        lat, lon, radius = location.lat, location.lon, location.radius
        gnd_itrs = self._geo_to_itrs(np.array([[lat, lon, 0]])).reshape((3,))

        los = self._ITRS - gnd_itrs

        # Find Earth tangent passing through horizontal origin and z-axis.
        st = [0, 0, radius / np.sin(np.radians(lat))] - gnd_itrs

        # Find LOS projection on tangent plane.
        los_on_n = np.sum(los * gnd_itrs)
        los_proj = los - np.array([los_on_n * gnd_itrs[i] for i in range(3)]) / radius ** 2

        # Determing azimuth.
        direc = np.cross(st, gnd_itrs)
        neg_ind = np.sum(los_proj * direc) < 0
        az = self._get_ang(st, los_proj)
        az = 360 * neg_ind - 1 * (2 * neg_ind - 1) * az

        return az

    def horizontal(self, location, stroke=False) -> Tuple:
        """Return horizontal coordinates in decimal degrees.

        The azimuth angle ranges from 0 to 360 degrees, measured clockwise
        from North. The altitude angle ranges from 0 to 90 degrees measured
        above the local horizon.

        Parameters
        ----------
        location : GroundPosition
            Eath bound origin of the horizontal system.
        stroke : bool, optional
            Formats output as a tuple of stroke objects if true. Otherwise,
            the output is a tuple of 1-D arrays.

        Returns
        -------
        tuple
            Tuple containing the altitude and azimuth angles.

            The angles are stroke objects if `stroke=True`. Otherwise, 1-D
            arrays are returned.

        Examples
        --------
        Calculate the horizontal coordinates at (52.1579, -106.6702):

        >>> julian=[2454545]
        >>> position = [[6343.82, -2640.87, -11.26]]
        >>> location = GroundPosition(52.1579, -106.6702)
        >>> c = Coordinate(position=position, frame="ecef", julian=julian)
        >>> c.horizontal(location=location)
        (array([-40.8786098]), array([94.73615482]))
        """

        if self._ITRS is None:
            self.itrs()

        # Calculate the altitude and azimuth angles.
        alt = self._altitude(location)
        az = self._azimuth(location)

        if stroke:
            return alt, az
        else:
            return alt(self._julian), az(self._julian)

    def off_nadir(self, location, stroke=False) -> np.ndarray:
        """Return off-nadir angle in decimal degrees.

        The off-nadir angle is the angular distance of a ground location from
        the satellites nadir.

        Parameters
        ----------
        location : GroundPosition
            Ground location of interest.
        stroke : bool, optional
            Formats output as a stroke objects if true. Otherwise, the output
            is a 1-D arrays.

        Returns
        -------
        Union[Stroke, np.ndarray]
            Off-nadir angles in decimal degrees.

            The off-nadir angles are returned as a stroke object if
            `stroke=True`. Otherwise, a 1-D array is returned.

        Examples
        --------
        Calculate the off-nadir angle for the location (52.1579, -106.6702):

        >>> toronto = GroundPosition(52.1579, -106.6702)

        >>> julian = [30462.70517808 30462.70583972]
        >>> position = [[-1148.75741119 -5527.65458256  3930.33125383]
        ...             [-1186.47076862 -5258.03277545  4276.52904538]]
        >>> c = Coordinate(position=position, frame="itrs", julian=julian, offset=2430000)
        >>> c.off_nadir(location=toronto)
        array([67.08711357 65.27748945])
        """

        if self._ITRS is None:
            self.itrs()

        lat, lon = location.lat, location.lon
        loc = self._geo_to_itrs(np.array([[lat, lon, 0]]))
        ang = self._get_ang((self._ITRS - loc).reshape((3,)), self._ITRS)

        if stroke:
            return ang
        else:
            return ang(self._julian)

    def _WGS84_radius(self, lattitude) -> np.ndarray:
        """Return Earth's geocentric radius using WGS84.

        Parameters
        ----------
        latitude : np.ndarray
            1-D array containing latitude in decimal degrees.

        Returns
        -------
        np.ndarray
            1-D array containing the geocentric radius in kilometres.

        Notes
        -----
        The radius is calculated using the WGS84 Earth ellipsoid model as
        detailed by "Earth Radius by Latitude" by Timur. [Tim18]_

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

        c, s = np.cos(phi), np.sin(phi)

        num = (a ** 2 * c) ** 2 + (b ** 2 * s) ** 2
        denom = (a * c) ** 2 + (b * s) ** 2
        radius = np.sqrt(num / denom)

        return radius

    def altitude(self, stroke=False) -> np.ndarray:
        """Return the geodetic altitude in kilometers.

        This method uses the WGS84 reference ellipsoid to calculate the
        geodetic altitude above the Earth's surface.

        Parameters
        ----------
        stroke : bool, optional
            Formats output as a stroke object if true. Otherwise, the output
            is a 1-D array.

        Returns
        -------
        Union[Stroke, np.ndarray]
            Geodetic altitude in kilometres.

            The geodetic altitude is returned as a stroke object if
            `stroke=True`. Otherwise, a 1-D array is returned.

        Notes
        -----
        This method uses an ellipsoid based model of the Earth to calculate
        the ellipsoid height in an iterative manner described in "Coordinate
        Systems in Geodesy" by E. J. Krakiwsky and D.E. Wells. [KW98b]_

        References
        ----------
        .. [KW98b] E. J. Krakiwsky and D. E. Wells. Coordinate Systems in
           Geodesy. Jan. 1998, pp. 31–33.

        Examples
        --------
        >>> julian=[2454545]
        >>> position = [[6343.82, -2640.87, -11.26]]
        >>> c = Coordinate(position=position, type="itrs", julian=julian)
        >>> c.altitude()
        np.array([504.1269764])
        """

        if self._ITRS is None:
            self.itrs()

        # Define WGS84 Parameters.
        a, b = 6378.137, 6356.752314245

        tol = 10e-10

        x, y, z = self._ITRS

        e = np.sqrt(1 - b ** 2 / a ** 2)
        p = np.sqrt(x ** 2 + y ** 2)

        h = np.sqrt(x ** 2 + y ** 2 + z ** 2) - np.sqrt(a * b)
        phi = np.arctan((z / p) * (1 - (e ** 2 * a) / (a + h)) ** -1)

        def new_vals(a, b, e, p, N, h, phi, z):

            N = a / np.sqrt(np.cos(phi) ** 2 + b ** 2 / a ** 2 * np.sin(phi) ** 2)
            h = p / np.cos(phi) - N
            phi = np.arctan((z / p) * (1 - (e ** 2 * N) / (N + h)) ** -1)

            return N, h, phi

        Np, hp, phip = new_vals(a, b, e, p, a, h, phi, z)

        while (np.mean((hp - h)(self._julian)) > a * tol) or (np.mean((phip - phi)(self._julian)) > tol):

            n, h, phi = Np, hp, phip
            Np, hp, phip = new_vals(a, b, e, p, n, h, phi, z)

        if stroke:
            return h
        else:
            return h(self._julian)

    def distance(self, location, stroke=False) -> np.ndarray:
        """Return distance to the ground location.

        Parameters
        ----------
        location : GroundPosition
        stroke : bool, optional
            Formats output as a stroke object if true. Otherwise, the output
            is a 1-D array.

        Returns
        -------
        Union[Stroke, np.ndarray]
            Distance in kilometres.

            The distance is returned as a stroke object if `stroke=True`.
            Otherwise, a 1-D array is returned.

        Examples
        --------
        Calculate the satellite distances from position (52.1579, -106.6702):

        >>> times = [30462.5, 30462.50069444]
        >>> position = [[6343.81620221, -2640.87223125, -11.25541802],
        ...             [6295.64583763, -2718.09271472, 443.08232543]]
        >>> location = GroundPosition(52.1579, -106.6702)
        >>> c = Coordinate(position=position, frame="itrs", julian=julian, offset=2430000)
        >>> c.distance(location=location)
        np.array([9070.49268746, 8776.7179543])
        """

        if self._ITRS is None:
            self.itrs()

        lat, lon = location.lat, location.lon
        gnd_itrs = self._geo_to_itrs(np.array([[lat, lon, 0]]))
        distance = np.linalg.norm(self._ITRS - gnd_itrs)

        if stroke:
            return distance
        else:
            return distance(self._julian)
