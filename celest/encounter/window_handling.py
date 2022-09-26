

from celest.coordinates.frames.attitude import Attitude
from celest.coordinates.ground_location import GroundLocation
from celest.file_save import TextFileWriter
from celest.units.quantity import Quantity
from celest import units as u
from typing import Union
import numpy as np


class VisibleTimeWindow:
    """VisibleTimeWindow(rise_time, set_time, attitude)

    A visible time window is the time at which a satellite is visible from a
    ground location. It is characterized by the time at which the satellite
    first becomes visible to when it is last visible.

    Parameters
    ----------
    rise_time : float
        The time at which the satellite is first visible in the jd2000 epoch.
    set_time : float
        The time at which the satellite is no longer visible in the jd2000
        epoch.
    attitude : Attitude
        The attitude of the satellite.

    Attributes
    ----------
    rise_time : Quantity
        The time at which the satellite is visible in the JD2000 frame.
    set_time : Quantity
        The time at which the satellite is no longer visible in the JD2000
        frame.
    attitude : Attitude
        The attitude of the satellite.
    """

    def __init__(self, rise_time: float, set_time: float, attitude:
                 Attitude) -> None:
        self._rise_time = Quantity(rise_time, u.jd2000)
        self._set_time = Quantity(set_time, u.jd2000)
        self._attitude = attitude

    def __str__(self) -> str:
        return f"Rise time: {self._rise_time}, Set time: {self._set_time}, "\
               f"Attitude: {self._attitude}"

    def __repr__(self) -> str:
        return f"VisibleTimeWindow({self._rise_time.data}, "\
               f"{self._set_time.data}, {self._attitude})"

    @property
    def rise_time(self) -> Quantity:
        return self._rise_time

    @property
    def set_time(self) -> Quantity:
        return self._set_time

    @property
    def attitude(self) -> Attitude:
        return self._attitude


class ObservationWindow:
    """ObservationWindow(start_time, duration, deadline, location, attitude)

    An observation window is a subset of a visible time window and represents
    the actual encounter to fulfill some mission request.

    Parameters
    ----------
    start_time : Quantity
        The start time of the encounter.
    duration : Quantity
        The duration of the encounter.
    location : GroundLocation
        The ground location involved in the encounter.
    attitude : Attitude
        The attitude of the satellite during the encounter.

    Attributes
    ----------
    start_time : Quantity
        The start time of the encounter in the JD2000 frame.
    duration : Quantity
        The duration of the encounter in seconds.
    location : GroundLocation
        The ground location involved in the encounter.
    attitude : Attitude
        The attitude of the satellite during the encounter.
    """

    def __init__(self, start_time: Quantity, duration: Quantity,
                 location: GroundLocation, attitude: Attitude) -> None:
        self._start_time = start_time
        self._duration = duration
        self._location = location
        self._attitude = attitude

    def __str__(self) -> str:
        return f"Start time: {self._start_time}, Duration: {self._duration}, "\
            f"Location: {self._location}, Attitude: {self._attitude}"

    def __repr__(self) -> str:
        return f"ObservationWindow({self._start_time}, {self._duration}, "\
            f"{self.location}, {self._attitude})"

    @property
    def start_time(self) -> Quantity:
        return self._start_time

    @property
    def duration(self) -> Quantity:
        return self._duration

    @property
    def location(self) -> GroundLocation:
        return self._location

    @property
    def attitude(self) -> Attitude:
        return self._attitude


class WindowCollection:
    """WindowCollection()

    Container to hold window data.

    The `WindowCollection` class can only hold one type of window data (i.e.
    either visible time windows or observation windows) at a time due to their
    differing attributes and physical significance.

    Methods
    -------
    add_window(window)
        Add a window to the container.
    save_text_file(filename)
        Save the window data to a text file.

    Raises
    ------
    TypeError
        If attempting to add both `VisibleTimeWindow` and `ObservationWindow`
        objects.
    """

    def __init__(self) -> None:
        self._window_data = []
        self._current_window_index = 0

    def __str__(self) -> str:
        return str(self._window_data)

    def __repr__(self) -> str:
        return "WindowCollection(" +\
               ", ".join([repr(i) for i in self._window_data]) + ")"

    def __len__(self) -> int:
        return len(self._window_data)

    def __iter__(self) -> Union[VisibleTimeWindow, ObservationWindow]:
        return self

    def __next__(self) -> Union[VisibleTimeWindow, ObservationWindow]:
        if self._current_window_index < len(self._window_data):
            self._current_window_index += 1
            return self._window_data[self._current_window_index - 1]
        else:
            self._current_window_index = 0
            raise StopIteration

    def __getitem__(self, key) -> Union[list, VisibleTimeWindow,
                                        ObservationWindow]:
        return self._window_data[key]

    def add_window(self, window: Union[VisibleTimeWindow,
                   ObservationWindow]) -> None:
        """Add a window to the list of windows.

        Parameters
        ----------
        window : {VisibleTimeWindow, ObservationWindow}
            The window to add to the list.

            All windows in the window handler must be of the same type.
        """

        if not isinstance(window, (VisibleTimeWindow, ObservationWindow)):
            raise TypeError("Window must be of type VisibleTimeWindow or "
                            "ObservationWindow.")
        if len(self._window_data) == 0:
            self._window_data.append(window)
        elif isinstance(window, self._window_data[0].__class__):
            self._window_data.append(window)
        else:
            raise TypeError("All windows added must be consistent types.")

    def save_text_file(self, file_name: str) -> None:
        """Save data as a pretty text file.

        Parameters
        ----------
        file_name : str
            Name of the text file to save data to.
        """

        if len(self._window_data) == 0:
            raise Exception("No windows to save.")

        header = self._window_data[0].__class__.__name__ + " Data"

        if isinstance(self._window_data[0], VisibleTimeWindow):
            data = self._get_visible_time_window_data()
        else:
            data = self._get_observation_window_data()

        writer = TextFileWriter(file_name, header)
        writer.add_layer(data=data)
        writer.save()

    def _get_visible_time_window_data(self) -> list:
        rise_time = []
        set_time = []

        for window in self._window_data:
            rise_time.append(window.rise_time.data)
            set_time.append(window.set_time.data)

        return [
            ["Rise Time", Quantity(np.array(rise_time), u.jd2000)],
            ["Set Time", Quantity(np.array(set_time), u.jd2000)]
        ]

    def _get_observation_window_data(self) -> list:
        start_time = []
        duration = []
        location = []
        roll = []
        pitch = []
        yaw = []

        for window in self._window_data:
            start_time.append(window.start_time.to(u.jd2000))
            duration.append(window.duration.to(u.s))
            location.append(window.location)
            roll.append(window.attitude.roll.to(u.deg)[0])
            pitch.append(window.attitude.pitch.to(u.deg)[0])
            yaw.append(window.attitude.yaw.to(u.deg)[0])

        return [
            ["Start Time", Quantity(np.array(start_time), u.jd2000)],
            ["Duration", Quantity(np.array(duration), u.s)],
            ["Location", np.array(location)],
            ["Roll", Quantity(np.array(roll), u.deg)],
            ["Pitch", Quantity(np.array(pitch), u.deg)],
            ["Yaw", Quantity(np.array(yaw), u.deg)]
        ]
