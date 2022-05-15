

from typing import Any, Literal
import numpy as np
import pandas as pd


class VTW:
    """VTW(rise_time, set_time, roll=None, pitch=None, yaw=None)

    Visible window time.

    A visibile window time describes an event where the satellite is in line of
    sight of the ground station within an elevation constraint.

    Parameters
    ----------
    rise_time : float
        Rise time of the visible time window in julian days.
    set_time : float
        Set time of the visible time window in julian days.
    roll : Stroke, optional
        Roll angle of the satellite in radians.
    pitch : Stroke, optional
        Pitch angle of the satellite in radians.
    yaw : Stroke, optional
        Yaw angle of the satellite in radians.

    Attributes
    ----------
    rise_time : float
        Rise time of the visible time window in julian days.
    set_time : float
        Set time of the visible time window in julian days.
    duration : float
        Duration of the visible time window in seconds.
    roll : Stroke
        Roll angle of the satellite in radians.
    pitch : Stroke
        Pitch angle of the satellite in radians.
    yaw : Stroke
        Yaw angle of the satellite in radians.
    """

    def __init__(self, rise_time, set_time, roll=None, pitch=None, yaw=None) -> None:

        self.rise_time = rise_time
        self.set_time = set_time
        self.duration = 86400 * (set_time - rise_time)
        self.roll = roll
        self.pitch = pitch
        self.yaw = yaw

    def __str__(self) -> str:

        return f"Rise: {self.rise_time}, Set:{self.set_time}, Duration: {self.duration}s"


class OW:
    """OW(start_time, duration, location, deadline, roll, pitch, yaw)

    Observation window.

    The observation window describes a selected time window in which the
    satellite-to-ground encounter is to be executed.

    Parameters
    ----------
    start_time : float
        Start time of the observing window in julian days.
    duration : float
        Duration of the observing window in seconds.
    location : GroundPosition
        Location of the satellite-ground encounter.
    deadline : float
        Deadline of the observing window in julian days.
    roll : float
        Roll angle of the satellite-ground encounter.
    pitch : float
        Pitch angle of the satellite-ground encounter.
    yaw : float
        Yaw angle of the satellite-ground encounter.

    Attributes
    ----------
    start : float
        Start time of the observing window in julian days.
    duration : float
        Duration of the observing window in seconds.
    location : GroundPosition
        Location of the satellite-ground encounter.
    deadline : float
        Deadline of the observing window in julian days.
    roll : float
        Roll angle of the satellite-ground encounter.
    pitch : float
        Pitch angle of the satellite-ground encounter.
    yaw : float
        Yaw angle of the satellite-ground encounter.
    """

    def __init__(self, start_time, duration, location, deadline, roll, pitch, yaw) -> None:

        self.start_time = start_time
        self.duration = duration
        self.location = location
        self.deadline = deadline
        self.roll = roll
        self.pitch = pitch
        self.yaw = yaw

    def __str__(self) -> str:

        return f"Location: {(self.location.latitude, self.location.longitude)}, Start: {self.start_time}, Duration:{self.duration}s"


class WindowHandler:
    """WindowHandler()

    Data structure to hold `VTW` or `OW` objects.
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

    def __getitem__(self, key) -> Literal[VTW, OW]:
        """Return window closest to key time.

        For a scalar index, the window with the closest start time will be
        returned. For a range index, all windows with start times within the
        range will be returned. For a tuple index, all the windows with start
        times closest to each tuple entry will be returned.
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

    def _window_start(self, window) -> float:

        return window.rise_time if isinstance(window, VTW) else window.start_time

    def _add_window_base(self, window) -> None:
        """Add new window to data base.

        Parameters
        ----------
        window : VTW or OW
        """

        temp_series = pd.Series(window, index=[self._window_start(window)])
        self.windows = pd.concat([self.windows, temp_series])
        self.windows.sort_index(inplace=True)

    def get_window(self, julian) -> Literal[VTW, OW]:
        """Return window closest to julian time.

        Parameters
        ----------
        julian : float
            Julian time in julian days.

        Returns
        -------
        VTW or OW
        """

        return self[julian]

    def get_windows_in_range(self, start, end) -> Literal[VTW, OW]:
        """Return all windows with start times within the range.

        Parameters
        ----------
        start : float
            Start time of the range in julian days.
        end : float
            End time of the range in julian days.

        Returns
        -------
        VTW or OW
        """

        return self[start:end]

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


class VTWHandler(WindowHandler):
    """VTWHandler()

    Data structure to hold `VTW` objects.
    """

    def __init__(self) -> None:

        super().__init__()
    
    def _add_window(self, window) -> None:
        """Add new window to data base.

        Parameters
        ----------
        window : VTW
        """

        if isinstance(window, VTW):
            self._add_window_base(window)
        else:
            raise TypeError("Expected VTW object.")

    def stats(self) -> pd.DataFrame:
        """Return visible time window statistics.

        This method returns statistics on the visible time windows, including
        the 5th, 50th, and 95th percentile of encounter durations and
        enocounter frequency.

        Return
        ------
        DataFrame
            Data base of encounter statistics.

        Examples
        --------
        Assuming `toronto_img` is a `VTWHandler` object returned from the
        windows generation function.

        >>> stats = toronto_img.stats()
        >>> print(stats)
                                              Statistic
        Minimum Duration (s)                        0.0
        Q5 Duration (s)                      108.468446
        Q50 Duration (s)                     167.687443
        Q75 Duration (s)                     177.506299
        Maximum Duration (s)                      524.0
        Encounter Frequency (per day)          1.340795
        """

        if len(self.windows) == 0:
            return None

        duration = []
        rise_times = []
        set_times = []

        for window in self.windows:
            duration.append(window.duration)
            rise_times.append(window.rise_time)
            set_times.append(window.set_time)

        data = {}
        data['Min Duration (s)'] = np.min(duration)
        data['Q5 Duration (s)'] = np.percentile(duration, 5)
        data['Q50 Duration (s)'] = np.percentile(duration, 50)
        data['Q95 Duration (s)'] = np.percentile(duration, 95)
        data['Max Duration (s)'] = np.max(duration)
        data['Encounter Frequency (per day)'] = len(self) / (np.max(set_times) - np.min(rise_times))

        df = pd.DataFrame(pd.Series(data), columns=["Statistic"])

        return df

    def save(self, fname, delimiter) -> None:
        """Save visible time window information.

        Parameters
        ----------
        fname : str
            File name to save encounter information under.
        delimiter : {",", "\\t"}
            String or character separating columns.
        """

        np_data = self.to_numpy()

        out_data = np.empty(shape=(np_data.shape[0] + 1, 3), dtype="<U25")
        out_data[0, :] = ["Rise Time (JD2000)", "Set Time (JD2000)", "Duration (s)"]

        for i, window in enumerate(np_data):

            out_data[i + 1, 0] = window.rise_time
            out_data[i + 1, 1] = window.set_time
            out_data[i + 1, 2] = window.duration

        np.savetxt(fname, out_data, delimiter=delimiter,
                   fmt="%s,%s,%s")


class OWHandler(WindowHandler):
    """OWHandler()

    Data structure to hold `OW` objects.
    """

    def __init__(self) -> None:

        super().__init__()
    
    def _add_window(self, window) -> None:
        """Add new window to data base.

        Parameters
        ----------
        window : OW
        """

        if isinstance(window, OW):
            self._add_window_base(window)
        else:
            raise TypeError("Expected OW object.")

    def save(self, fname, delimiter) -> None:
        """Save observation window information.

        Parameters
        ----------
        fname : str
            File name to save encounter information under.
        delimiter : {",", "\\t"}
            String or character separating columns.
        """

        np_data = self.to_numpy()

        out_data = np.empty(shape=(np_data.shape[0] + 1, 6), dtype="<U25")
        out_data[0, :] = [
            "Latitude (deg)",
            "Longitude (deg)",
            "Start Time (JD2000)",
            "Duration (s)",
            "Quality (1-10)",
            "Deadline (JD2000)"
        ]

        for i, window in enumerate(np_data):

            out_data[i + 1, 0] = window.location.latitude
            out_data[i + 1, 1] = window.location.longitude
            out_data[i + 1, 2] = window.start_time
            out_data[i + 1, 3] = window.duration
            out_data[i + 1, 4] = window.quality
            out_data[i + 1, 5] = window.deadline

        np.savetxt(fname, out_data, delimiter=delimiter,
                   fmt="%s,%s,%s,%s,%s,%s")
