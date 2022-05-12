

from celest.satellite._angle_representations import _ISO6709_representation
from celest.satellite.nutation_precession import (
    bias_matrix,
    precession_matrix,
    nutation_matrix
)
from celest.satellite.time import Time
from polare import Stroke
from typing import Tuple, Union
import numpy as np


class Coordinate(Time):
    """Coordinate(position, velocity, frame, julian, offset=0)

    Satellite coordinate transformations.

    The input position, velocity, and time can be converted into different
    reference frames useful for satellite applications. Here, `julian + offset`
    is the Julian time in J2000 epoch.

    The `position` and `velocity` inputs require three columns corresponding to
    the x, y, z cartesian coordinates.

    Parameters
    ----------
    position : array_like
        2-D array containing position coordinates in km. Axis 0 is the time
        series axis.

        Supported coordinate frames include the Geocentric Celestial Reference
        System (gcrs) and International Terrestrial Reference System (itrs).
    velocity : array_like
        2-D array containing velocity coordinates in m/s. Axis 0 is the time
        series axis.

        See `position` input for supported coordinate frames.
    frame : {"gcrs", "itrs"}
        Frame specifier for the input position.
    julian : array_like
        1-D array containing time data in Julian days.

        If times are not in the J2000 epoch, a non-zero offset must be pased
        in to add to the julian times.
    offset : float, optional
        Offset to convert input time data to the J2000 epoch, default is zero.

    Examples
    --------
    Initialize `Coordinate` using gcrs data:

    >>> julian = [30462.50000, 30462.50069]
    >>> position = [[-4681.50824, -5030.09119, 000.00000]
    ...             [-4714.35352, -4978.74326, 454.41493]]
    >>> velocity = [[-0.72067, 0.67072, 7.57919]
    ...             [-0.37378, 1.04024, 7.56237]]
    >>> c = Coordinate(position=position, velocity=velocity, frame="gcrs",
    ...                julian=julian, offset=2430000)

    Get horizontal coordinates for the ground location at (43.65, -79.38):

    >>> location = GroundPosition(latitude=43.65, longitude=-79.38, height=0.076)
    >>> elevation, azimuth = c.horizontal(location=location)
    """

    def __init__(self, position, velocity, frame, julian, offset=0) -> None:

        super().__init__(julian, offset)

        position = np.array(position)
        velocity = np.array(velocity)
        time = np.array(julian)

        if frame not in ["gcrs", "itrs"]:
            raise ValueError(f"{frame} is not a valid frame.")

        if position.ndim != 2:
            raise ValueError("position input must be 2-D.")

        if velocity.ndim != 2:
            raise ValueError("velocity input must be 2-D.")

        if time.ndim != 1:
            raise ValueError("time input must be 1-D.")

        if position.shape[0] != time.size:
            raise ValueError(f"position and time data lengths are mismatched:"
                             f" lengths are {len(position)} and {len(time)}")

        if velocity.shape[0] != time.size:
            raise ValueError(f"velocity and time data lengths are mismatched:"
                             f" lengths are {len(velocity)} and {len(time)}")

        kind_dict = {1: "zero", 2: "linear", 3: "quadratic"}
        self._kind = "cubic" if len(position) > 3 else kind_dict[len(position)]

        self._GCRS = None
        self._ITRS = None
        self._length = None

        self._set_base_data(position, velocity, frame)

    def __len__(self) -> int:

        return self._length

    def _stroke_init(self, data) -> Stroke:
        """Initialize data Stroke.

        Parameters
        ----------
        data : np.ndarray
            1-D array containing data.

        Returns
        -------
        Stroke
        """

        return Stroke(self._julian, data, self._kind)

    def _set_base_data(self, position, velocity, frame) -> None:
        """Initialize base data.

        The `position` and `velocity` inputs require three columns
        corresponding to the x, y, z cartesian coordinates.

        Parameters
        ----------
        position : array_like
            2-D array containing position coordinates.

            Supported coordinate frames include the Geocentric Celestial
            Reference System (gcrs) and International Terrestrial Reference
            System (itrs).
        velocity : array_like
            2-D array containing velocity coordinates.

            See `position` input for supported coordinate frames.
        frame : {"gcrs", "itrs"}
            Frame specifier for the input position.
        """

        self._length = len(position)

        if frame == "gcrs":
            x, y, z = [self._stroke_init(position[:, i]) for i in range(3)]
            vx, vy, vz = [self._stroke_init(velocity[:, i]) for i in range(3)]
            self._GCRS = np.array([x, y, z, vx, vy, vz])
        elif frame == "itrs":
            x, y, z = [self._stroke_init(position[:, i]) for i in range(3)]
            vx, vy, vz = [self._stroke_init(velocity[:, i]) for i in range(3)]
            self._ITRS = np.array([x, y, z, vx, vy, vz])
        else:
            raise ValueError(f"{frame} is not a valid frame.")

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

    def geo(self, iso=False, stroke=False) -> Union[Tuple, np.ndarray]:
        """Return geographical coordinates.

        The method cannot return data as both a stroke and iso
        formatted strings. If both `iso=True` and `stroke=True`, the Stroke
        return will take precedence.

        Parameters
        ----------
        iso : bool, optional
            Formats output as ISO6709 sexagesimal position strings if true.
        stroke : bool, optional
            Formats output as a tuple of stroke objects if true. Otherwise,
            the output is a tuple of 1-D arrays.

        Returns
        -------
        Union[Tuple, np.ndarray]
            Tuple of latitude, longitude, and altitude given in decimal
            degrees and km.

            If `iso=True`, the output is a 1-D array containing formatted
            strings.

            If `stroke=True`, the output is a tuple of stroke objects.

        See Also
        --------
        itrs : Return itrs coordinates.
        gcrs : Return gcrs coordinates.

        Examples
        --------
        Initialize `Coordinate` using gcrs data:

        >>> julian = [30462.50000, 30462.50069]
        >>> position = [[-4681.50824, -5030.09119, 000.00000]
        ...             [-4714.35352, -4978.74326, 454.41493]]
        >>> velocity = [[-0.72067, 0.67072, 7.57919]
        ...             [-0.37378, 1.04024, 7.56237]]
        >>> c = Coordinate(position=position, velocity=velocity, frame="gcrs",
        ...                julian=julian, offset=2430000)

        Get geographical data:

        >>> c.geo()
        (array([-0.09540, 3.69543]),
         array([-22.60391, -23.35431]),
         array([493.42306, 493.59707]))

        Return formatted output strings by setting `iso=True`:

        >>> c.geo(iso=True)
        array(['00°05′37.86″S 22°36′05.44″W 493.42km'
               '03°42′10.85″N 23°21′06.73″W 493.60km'])
        """

        if self._ITRS is None:
            x, y, z, vx, vy, vz = self._GCRS

            gcrs_pos = np.array([i(self._julian) for i in [x, y, z]]).T
            gcrs_vel = np.array([i(self._julian) for i in [vx, vy, vz]]).T

            itrs_pos = self._gcrs_and_itrs(gcrs_pos, "gcrs")
            x, y, z = [self._stroke_init(itrs_pos[:, i]) for i in range(3)]

            itrs_vel = self._gcrs_and_itrs(gcrs_vel, "gcrs")
            vx, vy, vz = [self._stroke_init(itrs_vel[:, i]) for i in range(3)]

            self._ITRS = np.array([x, y, z, vx, vy, vz])

        lat, lon = self._itrs_to_geo(self._ITRS[:3])
        alt = self.altitude(stroke=True)

        if stroke:
            return lat, lon, alt

        lat, lon, alt = lat(self._julian), lon(self._julian), alt(self._julian)

        if iso:
            return _ISO6709_representation(lat, lon, alt)
        else:
            return lat, lon, alt

    def era(self, stroke=False) -> Union[np.ndarray, Stroke]:
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
        Initialize `Coordinate` using gcrs data:

        >>> julian = [30462.50000, 30462.50069]
        >>> position = [[-4681.50824, -5030.09119, 000.00000]
        ...             [-4714.35352, -4978.74326, 454.41493]]
        >>> velocity = [[-0.72067, 0.67072, 7.57919]
        ...             [-0.37378, 1.04024, 7.56237]]
        >>> c = Coordinate(position=position, velocity=velocity, frame="gcrs",
        ...                julian=julian, offset=2430000)

        Get era data:

        >>> c.era()
        array([249.65772, 249.90840])
        """

        # Multiply time elapsed since J2000 by Earth rotation rate and add
        # J2000 orientation.
        dJulian = self._julian - 2451545
        ang = (360.9856123035484 * dJulian + 280.46) % 360

        return Stroke(self._julian, ang, self._kind) if stroke else ang

    def _gcrs_and_itrs(self, data, frame) -> np.ndarray:
        """Transform between gcrs and itrs coordinates.

        Parameters
        ----------
        data : np.ndarray
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
        due to their poor predictability.
        """

        pos_data = np.copy(data)

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

    def gcrs(self, stroke=False) -> Tuple:
        """Return gcrs coordinates.

        Parameters
        ----------
        stroke : bool, optional
            Formats output as a tuple of stroke objects if true. Otherwise,
            the output is a tuple of 1-D arrays.

        Returns
        -------
        Tuple
            Tuple containing x, y, z cartesian positions and velocities in the
            gcrs frame.

            The coordinates are stroke objects if `stroke=True`. Otherwise,
            1-D arrays are returned.

        See Also
        --------
        itrs : Return itrs coordinates.
        geo : Return geographical coordinates.

        Examples
        --------
        Initialize `Coordinate` using itrs data:

        >>> julian = [30462.50000, 30462.50069]
        >>> position = [[6343.70553, -2641.13728, -11.44111]
        ...             [6295.54131, -2718.36512, 442.89640]]
        >>> velocity = [[-0.37176, -0.92574, 7.57747]
        ...             [-0.84200, -0.72521, 7.56150]]
        >>> c = Coordinate(position=position, velocity=velocity, frame="gcrs",
        ...                julian=julian, offset=2430000)

        Get era data:

        >>> c.gcrs()
        (array([-4681.50824, -4714.35352]),
         array([-5030.09119, -4978.74326]),
         array([0, 454.41493]),
         array([-0.72067, -0.37378]))
         array([0.67072, 1.04024]),
         array([7.57919, 7.56237]))
        """

        if self._GCRS is None:
            x, y, z, vx, vy, vz = self._ITRS

            itrs_pos = np.array([i(self._julian) for i in [x, y, z]]).T
            itrs_vel = np.array([i(self._julian) for i in [vx, vy, vz]]).T

            gcrs_pos = self._gcrs_and_itrs(itrs_pos, "itrs")
            x, y, z = [self._stroke_init(gcrs_pos[:, i]) for i in range(3)]

            gcrs_vel = self._gcrs_and_itrs(itrs_vel, "itrs")
            vx, vy, vz = [self._stroke_init(gcrs_vel[:, i]) for i in range(3)]

            self._GCRS = np.array([x, y, z, vx, vy, vz])

        x, y, z, vx, vy, vz = self._GCRS

        if stroke:
            return x, y, z, vx, vy, vz
        else:
            return x(self._julian), y(self._julian), z(self._julian), \
                   vx(self._julian), vy(self._julian), vz(self._julian)

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
            Tuple containing x, y, z cartesian positions and velocities in the
            itrs frame.

            The coordinates are stroke objects if `stroke=True`. Otherwise,
            1-D arrays are returned.

        See Also
        --------
        gcrs : Return gcrs coordinates.
        geo : Return geographical coordinates.

        Examples
        --------
        Initialize `Coordinate` using gcrs data:

        >>> julian = [30462.50000, 30462.50069]
        >>> position = [[-4681.50824, -5030.09119, 000.00000]
        ...             [-4714.35352, -4978.74326, 454.41493]]
        >>> velocity = [[-0.72067, 0.67072, 7.57919]
        ...             [-0.37378, 1.04024, 7.56237]]
        >>> c = Coordinate(position=position, velocity=velocity, frame="gcrs",
        ...                julian=julian, offset=2430000)

        Get itrs data:

        >>> c.itrs()
        (array([[6343.70553, 6295.54131]),
         array([-2641.13728, -2718.36512]),
         array([-11.44111, 442.89640]),
         array([-0.37176, -0.84200]),
         array([-0.92574, -0.72521]),
         array([7.57747, 7.56150]))
        """

        if self._ITRS is None:
            x, y, z, vx, vy, vz = self._GCRS

            gcrs_pos = np.array([i(self._julian) for i in [x, y, z]]).T
            gcrs_vel = np.array([i(self._julian) for i in [vx, vy, vz]]).T

            itrs_pos = self._gcrs_and_itrs(gcrs_pos, "gcrs")
            x, y, z = [self._stroke_init(itrs_pos[:, i]) for i in range(3)]

            itrs_vel = self._gcrs_and_itrs(gcrs_vel, "gcrs")
            vx, vy, vz = [self._stroke_init(itrs_vel[:, i]) for i in range(3)]

            self._ITRS = np.array([x, y, z, vx, vy, vz])

        x, y, z, vx, vy, vz = self._ITRS

        if stroke:
            return x, y, z, vx, vy, vz
        else:
            return x(self._julian), y(self._julian), z(self._julian), \
                   vx(self._julian), vy(self._julian), vz(self._julian)

    def _gcrs_to_lvlh_matrix(self, stroke=False) -> np.ndarray:
        """Return the gcrs to lvlh rotation matrix.

        Parameters
        ----------
        stroke : bool, optional
            Formats output as a 2-D array of Stroke objects if true.
            Otherwise, the output is a 3-D array.

        Returns
        -------
        np.ndarray
            2-D rotation matrix containing Stroke objects if `stroke=True`.
            Otherwise, 3-D rotation matrix.
        """

        gcrs_x, gcrs_y, gcrs_z, gcrs_vx, gcrs_vy, gcrs_vz = self._GCRS
        r = np.array([gcrs_x, gcrs_y, gcrs_z])
        v = np.array([gcrs_vx, gcrs_vy, gcrs_vz])

        norm_r = np.linalg.norm(r)
        r_cross_v = np.cross(r, v)
        norm_r_cross_v = np.linalg.norm(r_cross_v)

        lvlh_z = - np.array([gcrs_x / norm_r,
                             gcrs_y / norm_r,
                             gcrs_z / norm_r])
        lvlh_y = - np.array([r_cross_v[0] / norm_r_cross_v,
                             r_cross_v[1] / norm_r_cross_v,
                             r_cross_v[2] / norm_r_cross_v])
        lvlh_x = np.cross(lvlh_y, lvlh_z)

        Aoi = np.array([lvlh_x, lvlh_y, lvlh_z])

        if stroke:
            return Aoi
        else:
            a11 = Aoi[0, 0](self._julian)
            a12 = Aoi[0, 1](self._julian)
            a13 = Aoi[0, 2](self._julian)
            a21 = Aoi[1, 0](self._julian)
            a22 = Aoi[1, 1](self._julian)
            a23 = Aoi[1, 2](self._julian)
            a31 = Aoi[2, 0](self._julian)
            a32 = Aoi[2, 1](self._julian)
            a33 = Aoi[2, 2](self._julian)
            return np.array([[a11, a12, a13],
                             [a21, a22, a23],
                             [a31, a32, a33]])

    def lvlh(self, stroke=False) -> Tuple:
        """Return the LVLH (Hill frame) coordinates.

        The local-vertical local-horizontal frame (also known as the Hill
        frame) is a body frame where the z-axis is algigned with the negative
        of the geocentric position vector, the y-axis is aligned with the
        negative orbit normal, and the x-axis completes the right handed triad.

        Parameters
        ----------
        stroke : bool, optional
            Formats output as a tuple of stroke objects if true. Otherwise,
            the output is a tuple of 1-D arrays.

        Returns
        -------
        Tuple
            Tuple containing x, y, z cartesian positions and velocities in the
            VLVH frame.

            The coordinates are stroke objects if `stroke=True`. Otherwise,
            1-D arrays are returned.

        Notes
        -----
        The LVLH frame definition was taken from NASA's technical memorandum
        on coordinate frames for the space shuttle program. [NASA1974]_

        References
        ----------
        .. [NASA1974] Coordinate Systems for the Space Shuttle Program, Lyndon
           B. Johnson Space Center, Houston, Texas 77058, Oct 1974, no. NASA
           TM X-58153.

        Examples
        --------
        Initialize `Coordinate` using gcrs data:

        >>> julian = [30462.50000, 30462.50069]
        >>> position = [[-4681.50824, -5030.09119, 000.00000]
        ...             [-4714.35352, -4978.74326, 454.41493]]
        >>> velocity = [[-0.72067, 0.67072, 7.57919]
        ...             [-0.37378, 1.04024, 7.56237]]
        >>> c = Coordinate(position=position, velocity=velocity, frame="gcrs",
        ...                julian=julian, offset=2430000)

        Get lvlh data:

        >>> c.lvlh()
        (array([0, -1.13687e-13]),
         array([-4.54747e-13, 1.06581e-13]),
         array([-6871.56000, -6871.64511]),
         array([7.64286, 7.64272]),
         array([-2.22045e-16, -1.11022e-16]),
         array([5.55112e-17, -2.83622e-03]))
        """

        if self._GCRS is None:
            self.gcrs()

        gcrs_x, gcrs_y, gcrs_z, gcrs_vx, gcrs_vy, gcrs_vz = self._GCRS
        r = np.array([gcrs_x, gcrs_y, gcrs_z])
        v = np.array([gcrs_vx, gcrs_vy, gcrs_vz])

        Aoi = self._gcrs_to_lvlh_matrix(stroke=True)

        lvlh_r = np.dot(Aoi, r)
        lvlh_v = np.dot(Aoi, v)

        x, y, z = lvlh_r
        vx, vy, vz = lvlh_v

        if stroke:
            return x, y, z, vx, vy, vz
        else:
            return x(self._julian), y(self._julian), z(self._julian), \
                   vx(self._julian), vy(self._julian), vz(self._julian)

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

    def _elevation(self, location) -> Stroke:
        """Return the elevation angle of a satellite.

        Parameters
        ----------
        location : GroundPosition
            Eath bound origin of the horizontal system.

        Returns
        -------
        Stroke
            Stroke containing the elevation of the satellite.
        """

        if self._ITRS is None:
            self.itrs()

        # Get origin of horizontal system.
        lat, lon, height = location.latitude, location.longitude, location.height
        gnd_itrs = self._geo_to_itrs(np.array([[lat, lon, height]])).reshape((3,))

        return 90 - self._get_ang(self._ITRS[:3] - gnd_itrs, gnd_itrs)

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

        if self._ITRS is None:
            self.itrs()

        # Get origin of horizontal system.
        lat, lon = location.latitude, location.longitude
        radius, height = location.radius, location.height
        gnd_itrs = self._geo_to_itrs(np.array([[lat, lon, height]])).reshape((3,))

        los = self._ITRS[:3] - gnd_itrs

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
        from North. The elevation angle ranges from 0 to 90 degrees measured
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
            Tuple containing the elevation and azimuth angles.

            The angles are stroke objects if `stroke=True`. Otherwise, 1-D
            arrays are returned.

        Examples
        --------
        Initialize `Coordinate` using gcrs data:

        >>> julian = [30462.50000, 30462.50069]
        >>> position = [[-4681.50824, -5030.09119, 000.00000]
        ...             [-4714.35352, -4978.74326, 454.41493]]
        >>> velocity = [[-0.72067, 0.67072, 7.57919]
        ...             [-0.37378, 1.04024, 7.56237]]
        >>> c = Coordinate(position=position, velocity=velocity, frame="gcrs",
        ...                julian=julian, offset=2430000)

        Calculate the horizontal coordinates at (52.16, -106.67):

        >>> location = GroundPosition(latitude=52.16, longitude=-106.67,
        ...                           height=0.482)
        >>> c.horizontal(location=location)
        (array([-40.88076, -39.01006]),
         array([94.73900, 92.99101]))
        """

        # Calculate the elevation and azimuth angles.
        alt = self._elevation(location)
        az = self._azimuth(location)

        if stroke:
            return alt, az
        else:
            return alt(self._julian), az(self._julian)

    def off_nadir(self, location, stroke=False) -> Union[Stroke, np.ndarray]:
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
        Initialize `Coordinate` using gcrs data:

        >>> julian = [30462.50000, 30462.50069]
        >>> position = [[-4681.50824, -5030.09119, 000.00000]
        ...             [-4714.35352, -4978.74326, 454.41493]]
        >>> velocity = [[-0.72067, 0.67072, 7.57919]
        ...             [-0.37378, 1.04024, 7.56237]]
        >>> c = Coordinate(position=position, velocity=velocity, frame="gcrs",
        ...                julian=julian, offset=2430000)

        Calculate the off-nadir data at (52.16, -106.67):

        >>> location = GroundPosition(latitude=52.16, longitude=-106.67,
        ...                           height=0.482)
        >>> c.off_nadir(location=location)
        array([67.10070, 65.29457])
        """

        if self._ITRS is None:
            self.itrs()

        lat, lon, height = location.latitude, location.longitude, location.height
        loc = self._geo_to_itrs(np.array([[lat, lon, height]]))
        ang = self._get_ang((self._ITRS[:3] - loc).reshape((3,)), self._ITRS[:3])

        return ang if stroke else ang(self._julian)

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

    def altitude(self, stroke=False) -> Union[Stroke, np.ndarray]:
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
        Initialize `Coordinate` using gcrs data:

        >>> julian = [30462.50000, 30462.50069]
        >>> position = [[-4681.50824, -5030.09119, 000.00000]
        ...             [-4714.35352, -4978.74326, 454.41493]]
        >>> velocity = [[-0.72067, 0.67072, 7.57919]
        ...             [-0.37378, 1.04024, 7.56237]]
        >>> c = Coordinate(position=position, velocity=velocity, frame="gcrs",
        ...                julian=julian, offset=2430000)

        Calculate altitude data:

        >>> c.altitude()
        array([493.42306, 493.59714])
        """

        if self._ITRS is None:
            self.itrs()

        # Define WGS84 Parameters.
        a, b = 6378.137, 6356.752314245

        tol = 10e-10

        x, y, z = self._ITRS[:3]

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

    def distance(self, location, stroke=False) -> Union[Stroke, np.ndarray]:
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
        Initialize `Coordinate` using gcrs data:

        >>> julian = [30462.50000, 30462.50069]
        >>> position = [[-4681.50824, -5030.09119, 000.00000]
        ...             [-4714.35352, -4978.74326, 454.41493]]
        >>> velocity = [[-0.72067, 0.67072, 7.57919]
        ...             [-0.37378, 1.04024, 7.56237]]
        >>> c = Coordinate(position=position, velocity=velocity, frame="gcrs",
        ...                julian=julian, offset=2430000)

        Calculate the satellite distances from position (52.16, -106.67):

        >>> location = GroundPosition(latitude=52.16, longitude=-106.67,
        ...                           height=0.482)
        >>> c.distance(location=location)
        array([9070.80824, 8777.02141])
        """

        if self._ITRS is None:
            self.itrs()

        lat, lon, height = location.latitude, location.longitude, location.height
        gnd_itrs = self._geo_to_itrs(np.array([[lat, lon, height]]))
        distance = np.linalg.norm(self._ITRS[:3] - gnd_itrs)

        return distance if stroke else distance(self._julian)
