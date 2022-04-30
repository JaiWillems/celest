

from celest.encounter.groundposition import GroundPosition
from typing import Any
import numpy as np
import pandas as pd


class Window:
    """Window(satellite, location, start, end, enc, ang, lighting)

    Encounter information.

    Parameters
    ----------
    satellite : Satellite
        Satellite for ground location encounter.
    location : GroundPosition
        Ground location for satellite encounter.
    start, end: float
        Window start and end times in J2000 Julian.
    enc : {"data link", "image"}
        Encounter type.
    ang : float
        Encounter constraint angle.
    lighting : {-1, 0, 1}
        Lighting constraint for night only, all time, or day only encounters.

    Attributes
    ----------
    satellite : Satellite
        Satellite for ground location encounter.
    coor : tuple
        Ground location (lat, lon) coordinates.
    start, end : float
        Window start and end times in J2000 Julian.
    dt : float
        Window duration in seconds.
    type : {"data link", "image"}
        Encounter type.
    ang : float
        Encounter constraint angle in degrees.
    lighting : {-1, 0, 1}
        Lighting constraint for night only, all time, or day only encounters.
    """

    def __init__(self, satellite, location, start, end, enc, ang, lighting) -> None:

        self.satellite = satellite
        self.coor = (location.lat, location.lon)

        self.start = start
        self.end = end
        self.duration = 86400 * (end - start)

        self.type = enc
        self.ang = ang
        self.lighting = lighting

    def __str__(self) -> str:

        str_1 = str(self.coor) + " "
        str_2 = self.type + " "
        str_3 = "ang:" + str(self.ang) + " "
        str_4 = "lighting:" + str(self.lighting) + " "
        str_5 = "start:" + str(self.start) + " "
        str_6 = "end:" + str(self.end)

        str_final = str_1 + str_2 + str_3 + str_4 + str_5 + str_6

        return str_final

    def copy(self) -> Any:
        """Return `Window` copy.

        Returns
        -------
        Window

        Examples
        --------
        Let `window_old` be a `Window` object.

        >>> window_new = window_old.copy()
        """

        sat = self.satellite
        gnd = GroundPosition(self.coor[0], self.coor[1])
        start = self.start
        end = self.end
        enc = self.type
        ang = self.ang
        lighting = self.lighting

        return_window = Window(sat, gnd, start, end, enc, ang, lighting)

        return return_window


class Windows(object):
    """Data structure to hold `Window` objects.

    This class is designed to hold encounter information and facilitate
    user-data interactions.
    """

    def __init__(self) -> None:

        self.windows = pd.Series(dtype="object")
        self.index_head = 0

    def __iter__(self) -> Any:

        return self

    def __next__(self) -> Any:

        if self.index_head == len(self.windows):
            self.index_head = 0
            raise StopIteration
        else:
            index = self.windows.index[self.index_head]
            self.index_head += 1

            return self.windows[index]

    def __getitem__(self, key) -> Window:
        """Return window closest to key time.

        This method allows for window indexing. If a scalar index is provided,
        the window with the closest start time will be returned. If the index
        is a range, all windows that have start times within the range will be
        returned. If the index is a tuple, then all the windows with start
        times closes to each tuple entry will be returned.
        """

        inds = self.windows.index.to_numpy()

        if isinstance(key, (int, float)):

            ind = np.argmin(np.abs(inds - key))
            rtn = self.windows[inds[ind]]

        elif isinstance(key, slice):

            rtn = self.windows[key].to_numpy()

        elif isinstance(key, tuple):

            inds_new = []
            for k in key:
                ind = np.argmin(np.abs(inds - k))
                inds_new.append(inds[ind])

            inds_new = list(set(inds_new))
            rtn = self.windows[inds_new].to_numpy()

        return rtn

    def __len__(self) -> int:

        return len(self.windows)

    def _add_window(self, window) -> None:
        """Add new window to data base.

        Parameters
        ----------
        window : Window

        Notes
        -----
        This method adds a new window to the data base by appending it to the
        internal Pandas Series object indexed by the window's start time.
        """

        temp_series = pd.Series(window, index=[window.start])
        self.windows = pd.concat([self.windows, temp_series])
        self.windows.sort_index(inplace=True)

    def stats(self) -> pd.DataFrame:
        """Return enocunter statistics.

        This method returns statistics on the stored encounters. The statistics
        include the 5th, 50th, and 95th percentile of encounter durations and
        enocounter frequency.

        Return
        ------
        DataFrame
            Data base of encounter statistics.

        Examples
        --------
        Assuming `toronto_img` is a `Windows` object returned from the windows
        generation function.

        >>> stats = toronto_img.stats()
        >>> print(stats)
                                       (45, -73), image
        Q5 Duration (s)                      108.468446
        Q50 Duration (s)                     167.687443
        Q75 Duration (s)                     177.506299
        Encounter Frequency (per day)          1.340795
        """

        if len(self.windows) == 0:
            return None

        enc_dict = {}
        start_times = np.empty((0,))
        end_times = np.empty((0,))

        for window in self.windows:

            if (window.coor, window.type) not in enc_dict.keys():
                enc_dict[window.coor, window.type] = np.array(
                    [window.duration])
            else:
                arr = enc_dict[window.coor, window.type]
                new_arr = np.append(arr, window.duration)
                enc_dict[window.coor, window.type] = new_arr

            start_times = np.append(start_times, window.start)
            end_times = np.append(end_times, window.end)

        out_data = {}
        num_days = np.max(end_times) - np.min(start_times)

        for key in enc_dict.keys():

            dt_data = enc_dict[key]
            enc_str = f"{key[0]}, {key[1]}"
            qs = np.percentile(dt_data, [5, 50, 75])
            freq = len(dt_data) / num_days

            out_data[enc_str] = pd.Series([qs[0], qs[1], qs[2], freq])

        df = pd.DataFrame.from_dict(out_data)
        df.index = ["Q5 Duration (s)",
                    "Q50 Duration (s)",
                    "Q75 Duration (s)",
                    "Encounter Frequency (per day)"]

        return df

    def to_numpy(self) -> np.ndarray:
        """Return windows within a NumPy array.

        Returns
        -------
        np.ndarray
            Array containing window data.

        Examples
        --------
        Assuming `toronto_img` is a `Windows` object returned from the windows
        generation function.

        >>> toronto_img.to_numpy()
        array([<celest.encounter._window_handling.Window object at 0x0000023FA0D8DFD0>,
               <celest.encounter._window_handling.Window object at 0x0000023FF0757940>,
               <celest.encounter._window_handling.Window object at 0x0000023FA10D5E50>,
               <celest.encounter._window_handling.Window object at 0x0000023FA10D5CA0>],
               dtype=object)
        """

        return self.windows.to_numpy()

    def save(self, fname, delimiter) -> None:
        """Save encounter information.

        Parameters
        ----------
        fname : str
            File name to save encounter information to.
        delimiter : {",", "\\t"}
            String or character separating columns.
        """

        np_data = self.to_numpy()

        out_data = np.empty(shape=(np_data.shape[0] + 1, 6), dtype="<U25")
        out_data[0, :] = [
            "latitude (deg)",
            "longitude (deg)",
            "start (jul)",
            "end (jul)",
            "duration (sec)",
            "type"
        ]

        for i, window in enumerate(np_data):

            out_data[i + 1, 0] = window.coor[0]
            out_data[i + 1, 1] = window.coor[1]
            out_data[i + 1, 2] = window.start
            out_data[i + 1, 3] = window.end
            out_data[i + 1, 4] = window.duration
            out_data[i + 1, 5] = window.type

        np.savetxt(fname, out_data, delimiter=delimiter,
                   fmt="%s,%s,%s,%s,%s,%s")
