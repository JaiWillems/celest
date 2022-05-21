

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


WGS84_MAJOR_AXIS = 6378.137
WGS84_MINOR_AXIS = 6356.752314245

EARTH_ROTATION_RATE = 360.9856123035484
ERA_AT_JD2000 = 280.46


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
        julian = np.array(julian)

        if frame not in ["gcrs", "itrs"]:
            raise ValueError(f"{frame} is not a valid frame.")
        if position.ndim != 2:
            raise ValueError("position input must be 2-D.")
        if velocity.ndim != 2:
            raise ValueError("velocity input must be 2-D.")
        if julian.ndim != 1:
            raise ValueError("time input must be 1-D.")
        if position.shape[0] != julian.size:
            raise ValueError("position and time lengths are mismatched")
        if velocity.shape[0] != julian.size:
            raise ValueError("velocity and time lengths are mismatched")

        self._set_interpolataion_kind(julian.size)
        self._set_base_data(position, velocity, frame)

    def __len__(self) -> int:

        return self._length

    def _set_interpolataion_kind(self, data_length):

        length_to_kind = {1: "zero", 2: "linear", 3: "quadratic"}
        self._kind = "cubic" if data_length > 3 else length_to_kind[data_length]

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

        base_coor = np.concatenate((position, velocity), axis=1)

        converted_position = self._gcrs_and_itrs(position, frame)
        converted_velocity = self._gcrs_and_itrs(velocity, frame)
        converted_coor = np.concatenate((converted_position, converted_velocity),
                                        axis=1)

        if frame == "gcrs":
            self._GCRS = self._stroke_array_from_2D_array_cols(base_coor)
            self._ITRS = self._stroke_array_from_2D_array_cols(converted_coor)
        elif frame == "itrs":
            self._GCRS = self._stroke_array_from_2D_array_cols(converted_coor)
            self._ITRS = self._stroke_array_from_2D_array_cols(base_coor)

    def _stroke_array_from_2D_array_cols(self, data):

        return np.array([self._stroke_init(i) for i in data.T])

    def _stroke_init(self, data):

        return Stroke(self._julian, data, self._kind)

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

        latitude, longitude = self._itrs_to_geo(self._ITRS[:3])
        altitude = self.altitude(True)

        if stroke:
            return latitude, longitude, altitude

        latitude = latitude(self._julian)
        longitude = longitude(self._julian)
        altitude = altitude(self._julian)

        if iso:
            return _ISO6709_representation(latitude, longitude, altitude)
        else:
            return latitude, longitude, altitude

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

        x, y, z = position
        latitude = np.arctan(z / np.sqrt(x ** 2 + y ** 2))
        longitude = np.arctan2(y, x)

        return np.degrees(latitude), np.degrees(longitude)

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

        latitude, longitude, height = position.T
        latitude, longitude = np.radians(latitude), np.radians(longitude)

        e = np.sqrt(1 - WGS84_MINOR_AXIS ** 2 / WGS84_MAJOR_AXIS ** 2)
        n = WGS84_MAJOR_AXIS / np.sqrt(1 - e ** 2 * np.sin(latitude) ** 2)

        x = (n + height) * np.cos(latitude) * np.cos(longitude)
        y = (n + height) * np.cos(latitude) * np.sin(longitude)
        z = (n * (1 - e ** 2) + height) * np.sin(latitude)

        return np.concatenate([i.reshape((-1, 1)) for i in [x, y, z]], axis=1)

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

        time_since_JD2000 = self._julian - 2451545
        rotation_angle = (EARTH_ROTATION_RATE * time_since_JD2000 +
                          ERA_AT_JD2000) % 360

        return self._stroke_init(rotation_angle) if stroke else rotation_angle

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

        return self._gcrs_to_itrs(data) if frame == "gcrs" else self._itrs_to_gcrs(data)

    def _gcrs_to_itrs(self, gcrs) -> np.ndarray:

        gcrs = np.copy(gcrs)

        earth_rotation_angle = self.era()
        bias_mtrx = bias_matrix()
        precession_mtrx = precession_matrix(self._julian)
        nutation_mtrx = nutation_matrix(self._julian)

        gcrs = np.einsum('ij, kj -> ki', bias_mtrx, gcrs)
        gcrs = np.einsum('ijk, ik -> ij', precession_mtrx, gcrs)
        gcrs = np.einsum('ijk, ik -> ij', nutation_mtrx, gcrs)

        gcrs = self._rotate_by_earth_rotation_angle(gcrs, -earth_rotation_angle)

        return gcrs

    def _rotate_by_earth_rotation_angle(self, data, degree_era):

        radian_era = np.radians(degree_era)
        c, s = np.cos(radian_era), np.sin(radian_era)
        c1, c2, c3 = np.copy(data).T

        data[:, 0] = c * c1 - s * c2
        data[:, 1] = s * c1 + c * c2
        data[:, 2] = c3

        return data

    def _itrs_to_gcrs(self, itrs) -> np.ndarray:

        itrs = np.copy(itrs)

        earth_rotation_angle = self.era()
        bias_mtrx = bias_matrix()
        precession_mtrx = precession_matrix(self._julian)
        nutation_mtrx = nutation_matrix(self._julian)

        itrs = self._rotate_by_earth_rotation_angle(itrs, earth_rotation_angle)

        itrs = np.einsum('ijk, ik -> ij', np.linalg.inv(nutation_mtrx), itrs)
        itrs = np.einsum('ijk, ik -> ij', np.linalg.inv(precession_mtrx), itrs)
        itrs = np.einsum('ij, kj -> ki', np.linalg.inv(bias_mtrx), itrs)

        return itrs

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

        return tuple(self._GCRS if stroke else [i(self._julian) for i in self._GCRS])

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

        return tuple(self._ITRS if stroke else [i(self._julian) for i in self._ITRS])

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

        gcrs_position, gcrs_velocity = self._GCRS[:3], self._GCRS[3:]

        Aoi = self._gcrs_to_lvlh_matrix(stroke=True)
        lvlh_position = np.dot(Aoi, gcrs_position)
        lvlh_velocity = np.dot(Aoi, gcrs_velocity)

        lvlh = np.concatenate((lvlh_position, lvlh_velocity), axis=0)

        return tuple(lvlh if stroke else [i(self._julian) for i in lvlh.T])

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

        position, velocity = self._GCRS[:3], self._GCRS[3:]

        norm_position = np.linalg.norm(position)
        position_cross_velocity = np.cross(position, velocity)
        n = np.linalg.norm(position_cross_velocity)

        lvlh_z = - np.array([i / norm_position for i in position])
        lvlh_y = - np.array([i / n for i in position_cross_velocity])
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

        numerator = np.sum(u * v, axis=(ua or va))
        denominator = np.linalg.norm(u, axis=ua) * np.linalg.norm(v, axis=va)
        ang = np.degrees(np.arccos(numerator / denominator))

        return ang

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

        el_az = [self._elevation(location), self._azimuth(location)]

        return tuple(el_az if stroke else [i(self._julian) for i in el_az])

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

        ground_itrs = self._calculate_ground_location_itrs(location)
        return 90 - self._get_ang(self._ITRS[:3] - ground_itrs, ground_itrs)
    
    def _calculate_ground_location_itrs(self, location):

        ground_geographic = np.array([[location.latitude, location.longitude,
                                       location.height]])

        return self._geo_to_itrs(ground_geographic).reshape((3,))


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

        latitude, radius = location.latitude, location.radius
        ground_itrs = self._calculate_ground_location_itrs(location)

        ground_to_satellite = self._ITRS[:3] - ground_itrs

        z_axis_vector = [0, 0, radius / np.sin(np.radians(latitude))]
        surface_tangent = z_axis_vector - ground_itrs

        sat_to_ground_on_normal = np.sum(ground_to_satellite * ground_itrs)
        misc_vector = np.array([sat_to_ground_on_normal * ground_itrs[i] for i in range(3)]) / radius ** 2
        sat_to_ground_on_plane = ground_to_satellite - misc_vector

        direction = np.cross(surface_tangent, ground_itrs)
        negative_indices = np.sum(sat_to_ground_on_plane * direction) < 0
        azimuth = self._get_ang(surface_tangent, sat_to_ground_on_plane)
        azimuth = 360 * negative_indices - 1 * (2 * negative_indices - 1) * azimuth

        return azimuth

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

        ground_itrs = self._calculate_ground_location_itrs(location)
        ground_to_satellite = (self._ITRS[:3] - ground_itrs).reshape((3,))
        off_nadir = self._get_ang(ground_to_satellite, self._ITRS[:3])

        return off_nadir if stroke else off_nadir(self._julian)

    def _WGS84_radius(self, latitude) -> np.ndarray:
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

        c, s = np.cos(np.radians(latitude)), np.sin(np.radians(latitude))

        num = (WGS84_MAJOR_AXIS ** 2 * c) ** 2 + (WGS84_MINOR_AXIS ** 2 * s) ** 2
        denom = (WGS84_MAJOR_AXIS * c) ** 2 + (WGS84_MINOR_AXIS * s) ** 2
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

        a, b = WGS84_MAJOR_AXIS, WGS84_MINOR_AXIS

        tol = 10e-10

        x, y, z, _, _, _ = self._ITRS

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

        return h if stroke else h(self._julian)

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

        ground_itrs = self._calculate_ground_location_itrs(location)
        distance = np.linalg.norm(self._ITRS[:3] - ground_itrs)

        return distance if stroke else distance(self._julian)
