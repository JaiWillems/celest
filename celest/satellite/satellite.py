

from celest.satellite.coordinate import Coordinate
from typing import Tuple
import numpy as np
import pandas as pd


class Satellite(Coordinate):
    """Satellite(position, velocity, frame, julian, offset=0)

    Satellite abstraction for satellite-ground encounters.

    `julian + offset` is the Julian time in the J2000 epoch associated with
    input positions.

    Parameters
    ----------
    position : array_like
        2-D array containing position coordinates.

        Supported coordinate frames include the Geocentric Celestial Reference
        System (gcrs) and International Terrestrial Reference System (itrs).
    velocity : array_like
        2-D array containing velocity coordinates.

        Refer to the `position` parameter for supported coordinate frames.
    frame : {"gcrs", "geo", "itrs"}
        Frame specifier for the input position.
    julian : array_like
        1-D array containing time data in Julian days.

        If times are not in the J2000 epoch, a non-zero offset must be passed
        in to add to the julian times.
    offset : float, optional
        Offset to convert input time data to the J2000 epoch, default is zero.

    Examples
    --------
    Initialize `Satellite` using gcrs data:

    >>> julian = [30462.50000, 30462.50069]
    >>> position = [[-4681.50824, -5030.09119, 000.00000]
    ...             [-4714.35352, -4978.74326, 454.41493]]
    >>> velocity = [[-0.72067, 0.67072, 7.57919]
    ...             [-0.37378, 1.04024, 7.56237]]
    >>> s = Satellite(position=position, velocity=velocity, frame="gcrs",
    ...               julian=julian, offset=2430000)
    """

    def __init__(self, position, velocity, frame, julian, offset=0) -> None:

        super().__init__(position, velocity, frame, julian, offset)

    def __len__(self) -> int:

        return self._length

    def attitude(self, location, stroke=False) -> Tuple:
        """Return satellite roll, pitch, and yaw angles in decimal degrees.

        This method returns the attitude angles required to align the
        satellite's nadir with a ground location for ground imaging.

        Parameters
        ----------
        location : GroundPosition
            Location for satellite attitude pointing.
        stroke : bool, optional
            Formats output as a tuple of stroke objects if true. Otherwise,
            the output is a tuple of 1-D arrays.

        Returns
        -------
        Tuple
            Tuple containing roll, pitch, and yaw angles in decimal degrees.

        Notes
        -----
        The methods used were taken from [adcs1]_.

        References
        ----------
        .. [adcs1] G. H. J. van Vuuren, “The design and simulation analysis of
           an attitude determination and control system for a small earth
           observation satellite,” Master of Engineering, Stellenbosch
           University, Stellenbosch, South Africa, Mar 2015.

        Examples
        --------
        Initialize `Satellite` using gcrs data:

        >>> julian = [30462.50000, 30462.50069]
        >>> position = [[-4681.50824, -5030.09119, 000.00000]
        ...             [-4714.35352, -4978.74326, 454.41493]]
        >>> velocity = [[-0.72067, 0.67072, 7.57919]
        ...             [-0.37378, 1.04024, 7.56237]]
        >>> s = Satellite(position=position, velocity=velocity, frame="gcrs",
        ...               julian=julian, offset=2430000)

        Get attitude data:

        >>> location = GroundPosition(latitude=43.65, longitude=-79.38,
        ...                           height=0.076)
        >>> roll, pitch, yaw = s.attitude(location)
        """

        x, y, z, _, _, _ = [i.reshape((-1, 1)) for i in self.lvlh()]
        satellite_position = np.concatenate((x, y, z), axis=1)

        ground_itrs = self._calculate_ground_location_itrs(location).reshape((1, -1))
        repeated_ground_itrs = np.repeat(ground_itrs, self._length, axis=0)
        ground_gcrs = self._itrs_to_gcrs(repeated_ground_itrs)

        Aoi = self._gcrs_to_lvlh_matrix()
        ground_lvlh = np.einsum('jki, ik -> ij', Aoi, ground_gcrs)

        satellite_position_norm = np.linalg.norm(satellite_position, axis=1)
        s = - satellite_position / satellite_position_norm[:, None]

        ground_norm = np.linalg.norm(ground_lvlh - satellite_position, axis=1)
        g = (ground_lvlh - satellite_position) / ground_norm[:, None]

        v, c = np.cross(s, g, axis=1), np.sum(s * g, axis=1)

        a32 = v[:, 0] + v[:, 1] * v[:, 2] / (1 + c)
        a31 = - v[:, 1] + v[:, 0] * v[:, 2] / (1 + c)
        a33 = 1 - (v[:, 0] ** 2 + v[:, 1] ** 2) / (1 + c)
        a12 = - v[:, 2] + v[:, 0] * v[:, 1] / (1 + c)
        a22 = 1 - (v[:, 0] ** 2 + v[:, 2] ** 2) / (1 + c)

        roll = - np.degrees(np.arcsin(a32))
        pitch = np.degrees(np.arctan2(a31, a33))
        yaw = np.degrees(np.arctan2(a12, a22))

        if stroke:
            return tuple([self._stroke_init(i) for i in (roll, pitch, yaw)])
        else:
            return roll, pitch, yaw

    def save(self, times=("julian",), positions=("gcrs",), path=None,
             sep=",", float_format=None) -> None:
        """Export satellite data.

        Parameters
        ----------
        times : Tuple, optional
            Tuple containing time representations to save.

            Available times include "julian", "ut1", "gmst", and "gast".
        positions : Tuple, optional
            Tuple containing position representations to save.

            Available positions include "gcrs", "itrs", and "geo".
        path : str, optional
            String or path object to save data to. If not provided, data is
            returned as a string.
        sep : str, optional
            String of length 1 representing the field delimiter for the
            output file.
        float_format : str, optional
            Format string for floating point numbers.

        Examples
        --------
        >>> positions = ["gcrs", "itrs"]
        >>> finch.save(fname="data.csv", positions=positions, delimiter=",")
        """

        key_mapping = {
            "julian": (self.julian, "Julian.J2000"),
            "ut1": (self.ut1, "UT1.hr"),
            "gmst": (self.gmst, "GMST.hr"),
            "gast": (self.gast, "GAST.hr"),
            "itrs": (self.itrs, "ITRS.x.km", "ITRS.y.km", "ITRS.z.km",
                     "ITRS.vx.m/s", "ITRS.vy.m/s", "ITRS.vz.m/s"),
            "gcrs": (self.gcrs, "GCRS.x.km", "GCRS.y.km", "GCRS.z.km",
                     "GCRS.vx.m/s", "GCRS.vy.m/s", "GCRS.vz.m/s"),
            "geo": (self.geo, "GEO.lat.deg", "GEO.lon.deg", "GEO.height.km")
        }

        for time in times:
            if time not in key_mapping:
                raise ValueError(f"{time} is not a valid time type.")

        for position in positions:
            if position not in key_mapping:
                raise ValueError(f"{position} is not a valid position type.")

        data = {}

        if isinstance(times, str):
            times = (times,)

        for time in times:
            time_information = key_mapping[time]
            time_data = time_information[0]()
            column = time_information[1]

            data[column] = time_data

        if isinstance(positions, str):
            positions = (positions,)

        for position in positions:
            position_information = key_mapping[position]
            position_data = position_information[0]()
            columns = position_information[1:]

            for i, column in enumerate(columns):
                data[column] = position_data[i]

        df = pd.DataFrame(data)
        df.to_csv(path_or_buf=path, sep=sep, float_format=float_format)
