"""Satellite orbital representations and coordinate conversions.

This module contains the `Satellite` class to localize satellite-related
information and functionality.
"""


from celest.core.decorators import set_module
from celest.core.interpolation import _interpolate
from typing import Any, Literal, Tuple
import pandas as pd
import numpy as np


@set_module('celest.satellite')
class Satellite(object):
    """Localize satellite information and functionality.

    The `Satellite` class represents a satellite, be it artificial or natural,
    and allows for the position to be represented with time through multiple
    representations.

    Parameters
    ----------
    position : Coordinate
        Satellite time-evolving positions.

    Attributes
    ----------
    time : Time
        Times associated with the satellite positions.
    position : Coordinate
        Position of the satellite.

    Methods
    -------
    interpolate(windows, factor)
        Interpolate satellite data for window times and positions.
    save_data(fname, delimiter, times, positions)
        Save satellite time and position data.
    """

    def __init__(self, position: Any) -> None:
        """Initialize attributes."""

        self.time = position.time
        self.position = position
        self.length = len(position)
    
    def __len__(self):
        """Return length of satellite position and time data."""

        return self.length
    
    def interpolate(self, windows: Any, factor: int=5) -> None:
        """Interpolate satellite data for window times and positions.
        
        This method takes a series of windows and will interpolate the
        satellites position and time data in the regions corresponding with
        window times.

        Parameters
        ----------
        windows : WindowList
            Windows defining the regions to interpolate.
        factor : int, optional
            The interpolation factor will increase the number of points in
            window regions by the factor `factor`. A strictly positive value
            should be used.
        
        Notes
        -----
        The desire for window defined interpolation stems from the need for
        more detailed data around encounter opportunities to get more precise
        start and end times as well as positions for satellite orientations.
        
        Examples
        --------
        If `IMG_windows` be a series of imaging encounters, we can interpolate
        within such regions as follows:

        >>> satellite.interpolate(windows=IMG_windows, factor=5)
        """
        
        indices = np.zeros((0,), dtype=int)
        time = self.time.julian()
        position = self.position.gcrs()

        for window in windows:
            start = window.start
            end = window.end

            temp_ind = np.where((start < time) & (time < end))[0]
            indices = np.concatenate((indices, temp_ind))
        
        indices = np.split(indices, np.where(np.diff(indices) > 1)[0])
        indices = np.array(indices, dtype=object)

        time_interp = _interpolate(data=time, factor=factor, dt=2, indices=indices)
        position_interp = _interpolate(data=position, factor=factor, dt=2, indices=indices)

        self.time._julian = time_interp
        self.time._length = time_interp.size

        self.position.time = self.time
        self.position._GCRS = position_interp
        self.position._ITRS = None
        self.position._GEO = None
        self.position.length = time_interp.size

        self.length = len(self.position)

    def save_data(self, fname: str, delimiter: Literal[",", "\\t"], times:
                  Tuple=("julian",), positions: Tuple=("gcrs",)) -> None:
        """Save satellite time and position data.

        Parameters
        ----------
        fname : str
            File name of the output file as either a .txt or .csv file.
        delimiter : str
            String of length 1 representing the field delimiter for the output
            file.
        times : Tuple, optional
            List containing the types of time data to save. Possible list
            values include "julian", "ut1", "gmst", and "gast".
        positions : Tuple, optional
            List containing the types of position data to save. Possible
            list values include "gcrs", "itrs", and "geo".

        Notes
        -----
        It is recommended to use a tab delimiter for .txt files and comma
        delimiters for .csv files. The method will return an error if the
        fileName already exists in the current working directory.

        Examples
        --------
        >>> positions = ["gcrs", "itrs"]
        >>> finch.save_data(fname="data.csv", delimiter=",", positions=positions)
        """

        key_mapping = {
            "julian": (self.time.julian, "Julian.J2000"),
            "ut1": (self.time.ut1, "UT1.hr"),
            "gmst": (self.time.gmst, "GMST.hr"),
            "gast": (self.time.gast, "GAST.hr"),
            "itrs": (self.position.itrs, "ITRS.x.km", "ITRS.y.km", "ITRS.z.km"),
            "gcrs": (self.position.gcrs, "GCRS.x.km", "GCRS.y.km", "GCRS.z.km"),
            "geo": (self.position.geo, "GEO.lat.deg", "GEO.lon.deg", "GEO.height.km")
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
                data[column] = position_data[:, i]

        df = pd.DataFrame(data)
        df.to_csv(fname, sep=delimiter)
