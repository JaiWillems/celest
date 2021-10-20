"""Window information localization and handling.

This module contains the `Window` class to isolate window information in
addition to the `Windows` class to handle `Window` objects in a linked
list structure.
"""


from typing import Any, Literal
from celest.encounter.groundposition import GroundPosition
import pandas as pd
import numpy as np


class Window(object):
    """Isolate window information.

    Parameters
    ----------
    satellite : Satellite
        Satellite associated with the encounter window.
    location : GroundPosition
        Ground location for imaging or ground station communcations.
    start : float
        Window start time in Julian.
    end : float
        Window end time in Julian.
    enc : {"data link", "image"}
        Type of encounter.
    ang : float
        Encounter constraint angle.
    lighting : {-1, 0, 1}
        Lighting constraint to specify night only, any time, or day only
        encounters.
    sca : float
        Solar constraint angle.

    Attributes
    ----------
    satellite : Satellite
        Satellite associated with the encounter window.
    coor : tuple
        Ground location coordinate tuple.
    start : float
        Window start time in Julian.
    end : float
        Window end time in Julian.
    dt : float
        Window duration in seconds.
    type : {"data link", "image"}
        Type of encounter.
    angle_type : {"A", "N"}
        Employed constraint angle as either being the off-nadir or altitude
        angle.
    ang : float
        Encounter constraint angle.
    lighting : {-1, 0, 1}
        Lighting constraint to specify night only, any time, or day only
        encounters.
    sca : float
        Solar constraint angle.

    Methods
    -------
    copy()
        Copy local window data to a new `Window` object.
    """

    def __init__(self, satellite: Any, location: Any, start: float, end: float,
                 enc: Literal["data link", "image"], ang: float, lighting:
                 Literal[-1, 0, 1], sca: float) -> None:
        """Initialize attributes."""

        self.satellite = satellite
        self.coor = (location.lat, location.lon)

        self.start = start
        self.end = end
        self.duration = 86400 * (end - start)

        self.type = enc
        self.angle_type = "N" if enc == "image" else "A"
        self.ang = ang
        self.lighting = lighting
        self.sca = sca

    def __str__(self) -> str:
        """Define `Window` information string."""

        str_1 = str(self.coor) + "_"
        str_2 = self.type + "_"
        str_3 = "ang:" + self.angle_type + str(self.ang) + "_"
        str_4 = "lighting:" + str(self.lighting) + "_"
        str_5 = "sca:" + str(self.sca) + "_"
        str_6 = "start:" + str(self.start) + "_"
        str_7 = "end:" + str(self.end)

        str_final = str_1 + str_2 + str_3 + str_4 + str_5 + str_6 + str_7

        return str_final

    def copy(self) -> Any:
        """Return `Window` copy.

        Returns
        -------
        Window
            Return `Window` with the same attributes as `self`.

        Examples
        --------
        Let `window_old` be a `Window` object.

        >>> windoe_new = window_old.copy()
        """

        sat = self.satellite
        gnd = GroundPosition(self.coor[0], self.coor[1])
        start = self.start
        end = self.end
        enc = self.type
        ang = self.ang
        lighting = self.lighting
        sca = self.sca

        return_window = Window(sat, gnd, start, end, enc, ang, lighting, sca)

        return return_window


class Windows(object):

    def __init__(self) -> None:
        """Initialize attributes."""

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

    def __len__(self) -> int:
        """Return number of windows."""

        return len(self.windows)

    def _add_window(self, window: Window) -> None:
        """Add new window to data base.
        
        Parameters
        ----------
        window : Window
            Window to add to data base.
        """
        
        temp_series = pd.Series(window, index=[window.start])
        self.windows = self.windows.append(temp_series)
        self.windows.sort_index(inplace=True)

    def encounter_stats(self) -> pd.DataFrame:
        """Return enocunter statistics.
        
        Return
        ------
        DataFrame
            Data base of encounter statistics.
        """

        enc_dict = {}
        start_times = np.empty((0,))
        end_times = np.empty((0,))
        
        for window in self.windows:

            if (window.coor, window.type) not in enc_dict.keys():
                enc_dict[window.coor, window.type] = np.array([window.duration])
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

    def to_numpy(self) -> np.array:
        """Return windows within a NumPy array.
        
        Returns
        -------
        np.ndarray
            Array containing window data.
        """

        return self.windows.to_numpy()

    def save_encounters(self, fname: str, delimiter: Literal[",", "\t"]) -> None:
        """Save encounter information.
        
        Parameters
        ----------
        fname : str
            File name to save encounter information to.
        delimiter : {",", "\t"}
            String or character separating columns.
        """
        
        np_data = self.to_numpy()

        out_data = np.empty(shape=(np_data.shape[0], 5), dtype="<U25")

        for i, window in enumerate(np_data):
            
            out_data[i, 0] = window.coor[0]
            out_data[i, 1] = window.coor[1]
            out_data[i, 2] = window.start
            out_data[i, 3] = window.end
            out_data[i, 4] = window.duration
            out_data[i, 5] = window.type

        np.savetxt(out_data, fname=fname, delimiter=delimiter)