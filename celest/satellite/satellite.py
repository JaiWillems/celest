

from celest.satellite.coordinate import Coordinate
from typing import Literal, Tuple
import pandas as pd
import numpy.typing as npt


class Satellite(Coordinate):
    """Satellite(position, frame, julian, offset=0)

    Satellite abstraction for satellite-ground encounters.

    `julian + offset` is the Julian time in the J2000 epoch associated with
    input positions. If `frame=="geo"`, the `position` input can have 2
    columns (latitude, longitude) or three columns (latitude, longitude,
    altitude). If no geodetic altitude is provided, it is assumed zero. Other
    frames require three columns.

    Parameters
    ----------
    position : array_like
        2-D array containing position coordinates.

        Supported coordinate frames include the Geocentric Celestial Reference
        System (gcrs), International Terrestrial Reference System (itrs), and
        geographical (geo) system.
    frame : {"gcrs", "geo", "itrs"}
        Frame specifier for the input position.
    julian : array_like
        1-D array containing time data in Julian days.

        If times are not in the J2000 epoch, a non-zero offset must be passed
        in to add to the julian times.
    offset : float, optional
        Offset to convert input time data to the J2000 epoch, default is zero.

    Methods
    -------
    save_data(fname, delimiter, times, positions)
        Save satellite time and position data.
    """

    def __init__(self, position: npt.ArrayLike, frame: Literal["gcrs", "geo",
                 "itrs"], julian: npt.ArrayLike, offset=0) -> None:

        super().__init__(position, frame, julian, offset)

    def __len__(self):

        return self._length

    def save_data(self, times: Tuple=("julian",), positions: Tuple=("gcrs",),
                  path: str=None, sep=",", float_format: str=None) -> None:
        """Save satellite data.

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
        >>> finch.save_data(fname="data.csv", positions=positions, delimiter=",")
        """

        key_mapping = {
            "julian": (self.julian, "Julian.J2000"),
            "ut1": (self.ut1, "UT1.hr"),
            "gmst": (self.gmst, "GMST.hr"),
            "gast": (self.gast, "GAST.hr"),
            "itrs": (self.itrs, "ITRS.x.km", "ITRS.y.km", "ITRS.z.km"),
            "gcrs": (self.gcrs, "GCRS.x.km", "GCRS.y.km", "GCRS.z.km"),
            "geo": (self.geo, "GEO.lat.deg", "GEO.lon.deg", "GEO.height.km")
        }

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
