

from celest.core.interpolation import _interpolate
from celest.satellite.coordinate import Coordinate
from typing import Any, Literal, Tuple
import pandas as pd
import numpy as np
import numpy.typing as npt


class Satellite(Coordinate):
    """Satellite(position, frame, julian, offset=0)
    
    Satellite abstraction for satellite-ground encounters.

    `julian + offset` is the Julian time in the J2000 epoch associated with
    input positions. If `frame=="geo"`, the `position` input can have 2
    columns, (latitude, longitude) or three columns (lattitude, longitude,
    altitude). If no geodetic altitude is provided, it is assumed zero. Other
    frames require 3 columns.

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

        If times are not in the J2000 epoch, a non-zero offset must be pased
        in to add to the julian times.
    offset : float, optional
        Offset to convert input time data to the J2000 epoch, default is zero.

    Methods
    -------
    interpolate(windows, factor)
        Interpolate satellite data for window times and positions.
    save_data(fname, delimiter, times, positions)
        Save satellite time and position data.
    """

    def __init__(self, position: npt.ArrayLike, frame: Literal["gcrs", "geo",
                 "itrs"], julian: npt.ArrayLike, offset=0) -> None:

        super().__init__(position, frame, julian, offset)
    
    def __len__(self):

        return self._length
    
    def interpolate(self, windows: Any, factor: int=5) -> None:
        """Interpolate coordinate information around window times.

        This method interpolates position and time data around input encounter
        times to gain greater precision where needed.

        Parameters
        ----------
        windows : WindowList
            Windows defining interpolation regions.
        factor : int, optional
            The interpolation factor will increase the number of points in
            window regions by the factor `factor`.
        
        Notes
        -----
        Window defined interpolation is necessary for more precise coordinate
        information around encounter regions as might be used for satellite
        orientation calculations.
        
        Examples
        --------
        If `IMG_windows` be a series of imaging encounters, we can interpolate
        within such regions:

        >>> satellite.interpolate(windows=IMG_windows, factor=5)
        """

        if not isinstance(factor, int):
            raise ValueError("factor must be integer type.")
        
        if factor <= 0:
            raise ValueError("factor must be strictly positive.")
        
        ind = np.zeros((0,), dtype=int)
        time = self._julian
        position = self.gcrs()

        for window in windows:
            start = window.start
            end = window.end

            temp_ind = np.where((start < time) & (time < end))[0]
            ind = np.concatenate((ind, temp_ind))
        
        ind = np.split(ind, np.where(np.diff(ind) > 1)[0])
        ind = np.array(ind, dtype=object)

        time_interp = _interpolate(time, factor, 2, ind)
        position_interp = _interpolate(position, factor, 2, ind)

        self._julian = time_interp
        self._GCRS = position_interp
        self._ITRS = None
        self._GEO = None
        self._length = time_interp.size

    def save_data(self, fname: str, times: Tuple=("julian",), positions:
                  Tuple=("gcrs",), fmt="%.18e", delimiter=" ", newline="\n",
                  header="", footer="", comments="# ", encoding=None) -> None:
        """Save satellite data.

        Parameters
        ----------
        fname : str
        times : Tuple, optional
            Tuple containing time representations to save.
            
            Available times include "julian", "ut1", "gmst", and "gast".
        positions : Tuple, optional
            Tuple containing position representations to save.
            
            Available positions include "gcrs", "itrs", and "geo".
        fmt : str or sequence of strings, optional
            A single format or sequence of formats, or a multiformat string
            (in which case the delimiter is ignored).
        delimiter : str, optional
            String or character separating columns.
        newline : str, optional
            String or character separating lines.
        header : str, optional
            String to be written at the start of the file.
        footer : str, optional
            String to be written at the end of the file.
        comments : str, optional
            String to preappend to the header and footer to mark them as
            comments.
        encoding : {None, str}, optional
            Encoding used to encode the output file.

            Refer to NumPy's `savetxt` documentation for more information on
            encoding options.

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
                data[column] = position_data[:, i]

        df = pd.DataFrame(data)
        df.to_csv(fname, sep=delimiter)
