"""Satellite encounter planning.

The encounter module contains the `Encounter` class that is used to compute,
store, and schedule Earth-satellite encounters.
"""


from celest.satellite.coordinate import Coordinate
from celest.astronomy import Sun
from celest.core.decorators import set_module
from celest.encounter import GroundPosition
from celest.satellite import Satellite
from datetime import datetime
from typing import Literal
import julian
import numpy as np
import pandas as pd


@set_module('celest.encounter')
class Encounter(object):
    """Computes, schedules, and store satellite encounter data.

    The `Encounter` class extends off the information stored in `Satellite` and
    `GroundPosition` objects to compute, store, and schedule Earth-Satellite
    encounters.

    Parameters
    ----------
    satellite : Satellite
        Satellite of interest for ground based encounters.

    Attributes
    ----------
    gs : Dict
        Dictionary of `GroundPosition` objects with encounters initialized.

    Methods
    -------
    add_position(groundPos)
        Add a grouns position for encounter calculations.
    windows(interp=True, factor=5, buffer=10, dt=1)
        Initialize the `GroundPosition.encounters.windows` attribute.
    window_encounter_indices(buffer=0)
        Initialize `GroundPosition.encounters.encounter_indices` attributes.
    encounter_stats()
        Return encounter statistics.
    save_windows(fileName, delimiter)
        Save window data to local directory.
    """

    def __init__(self, satellite: Satellite) -> None:
        """Initialize attribites."""

        self._satellite = satellite
        self._sun_position = None
        self.gs = {}

    def add_position(self, groundPos: GroundPosition) -> None:
        """Add a ground position for encounter calculations.

        This function adds `groundPos` to the `gs` attribute dictionary where
        the key is `groundPos.name` and the key is `groundPos` itself. The
        `gs` attribute will be sifted through when determining encounter
        windows.

        Parameters
        ----------
        groundPos : GroundPosition
            `GroundPosition` object instantiated as per its documentation.
        """

        self.gs[groundPos.name, groundPos]

    def _solar_angles(self, ECEFdata: np.array, timeData: np.array, groundPos: GroundPosition) -> np.array:
        """Return solar angles.

        Parameters
        ----------
        ECEFdata : np.array
            Satellite ECEF position data.
        timeData : np.array
            Satellite julian time data.
        groundPos : GroundPosition
            Location of the communcations ground station.

        Returns
        -------
        np.array
            Array of shape (n,) containing degree solar angles.
        """

        gnd_GEO = np.array([[groundPos.coor[0], groundPos.coor[1], groundPos.radius]])
        gnd_GEO = np.repeat(gnd_GEO, timeData.size, axis=0)
        gnd_ECEF = Coordinate(gnd_GEO, "GEO", timeData).ECEF()

        LOS = ECEFdata - gnd_ECEF
        sun_ECEF = Sun().position(timeData).ECEF()

        dividend = np.sum(LOS, sun_ECEF, axis=1)
        divisor = np.linalg.norm(LOS, axis=1) * np.linalg.norm(sun_ECEF, axis=1)
        arg = np.divide(dividend, divisor)
        ang = np.degrees(np.arccos(arg))

        return ang

    def windows(self, interp: bool=True, factor: int=5, buffer: float=10, dt: int=1) -> None:
        """Initialize the `GroundPosition.encounters.windows` attribute.

        This function iterates through all `gs` values (`GroundPosition`
        objects) and the values in their `encounters` attributes
        (`EncounterSpec` objects) and initializes the `windows` attribute.

        Parameters
        ----------
        interp : {True, False}, optional
            If `interp=True` then the windows will be generated from data
            interpolated around valid encounter regions.
        factor : int, optional
            The factor increase in the number of steps in interpolated regions.
        buffer : float, optional
            Positive angular buffer to interpolate data that comes within the
            buffer from the encounter constraint angle.
        dt : int, optional
            The number of data points adjacent to a valid encounter region to
            interpolate within.

        Exampls:
        --------
        >>> encounters.windows(finch)

        See documentation for more detail on instantiated dependencies.
        """

        for pos in self.gs:
            for enc in self.gs[pos].encounters:

                ground_pos = self.gs[pos]
                enc = ground_pos.encounters[enc]
                ang = enc.ang
                ang_type = enc.ang_type
                solar = enc.solar
                SCA = enc.SCA

                if interp:
                    encounter_ind = self.window_encounter_indices(buffer=buffer)
                    times = self._satellite.times.julian(factor=factor, dt=dt, indices=encounter_ind)
                    ECEFdata = self._satellite.position.ECEF(factor=factor, dt=dt, indices=encounter_ind)
                    alt = self._satellite.position.horizontal(ground_pos, factor=factor, dt=dt, indices=encounter_ind)[:, 0]
                    nadir = self._satellite.position.off_nadir(ground_pos, factor=factor, dt=dt, indices=encounter_ind)
                else:
                    times = self._satellite.times.julian()
                    ECEFdata = self._satellite.position.ECEF()
                    alt = self._satellite.position.horizontal(ground_pos)[:, 0]
                    nadir = self._satellite.position.off_nadir(ground_pos)

                if not ang_type:
                    sun_angs = self._solar_angles(ECEFdata, times, ground_pos)
                    win_ind = np.where((alt > ang) & (sun_angs > SCA))[0]
                elif ang_type:
                    win_ind = np.where((nadir < ang) & (alt >= 0))[0]

                # Get sun position vector.
                if solar != 0:
                    sun_ECEF = Sun().position(timeData=times).ECEF()
                    gnd_ECEF = np.concatenate(np.array(list(ground_pos.coor), np.array(ground_pos.radius)))
                    gnd_ECEF = np.full((sun_ECEF.shape[0], 3), gnd_ECEF)
                    dividend = np.einsum("ij, ij->i", sun_ECEF, gnd_ECEF)
                    divisor = np.linalg.norm(sun_ECEF, axis=1) * np.linalg.norm(gnd_ECEF, axis=1)
                    arg = np.divide(dividend, divisor)
                    ang = np.degrees(np.arccos(arg))

                # Find intersection of night indices and window indices.
                if solar == -1:
                    night_ind = np.where(ang >= 90)
                    win_ind = np.intersect1d(win_ind, night_ind)

                # Find intersection of day indices and window indices.
                if solar == 1:
                    day_ind = np.where(ang < 90)
                    win_ind = np.intersect1d(win_ind, day_ind)

                win_ind_arr = np.split(win_ind, np.where(np.diff(win_ind) != 1)[0]+1)
                windows = np.empty((len(win_ind_arr), 5), dtype="U25")

                for j in range(len(win_ind_arr)):

                    if j == 0:
                        break

                    start = julian.from_jd(times[win_ind_arr[j][0]])
                    end = julian.from_jd(times[win_ind_arr[j][-1]])
                    max_alt = np.max(alt[win_ind_arr[j]])
                    min_nadir = np.min(nadir[win_ind_arr[j]])

                    windows[j, 0] = start
                    windows[j, 1] = end
                    windows[j, 2] = (end-start).total_seconds()
                    windows[j, 3] = max_alt
                    windows[j, 4] = min_nadir

                enc.windows = windows
                enc.length = windows.shape[0]

    def window_encounter_indices(self, buffer: float=0) -> None:
        """Initialize `GroundPosition.encounters.encounter_indices` attributes.

        The encounter indices are the indices of position and time data that
        are involved within an encounter. This function iterates through all
        `gs` values (`GroundPosition` objects`) and the values in their
        `encounters` attributes (`EncounterSpec` objects) and initializes the
        `encounter_indices` attribute.

        Parameters
        ----------
        buffer : float
            Positive angular buffer to interpolate data that comes within the
            buffer from the encounter constraint angle.
        """

        for pos in self.gs:
            for enc in self.gs[pos].encounters:

                # Get encounter information.
                ang_type = self.gs[pos].encounters[enc].ang_type
                ang = self.gs[pos].encounters[enc].ang
                is_max = self.gs[pos].encounters[enc].max_ang
                ground_pos = self.gs[pos]

                # Get ground position information.
                alt = ground_pos.alt
                nadir = ground_pos.nadirAng

                # Derive interpolation regions.
                if not ang_type and is_max:
                    regions = np.where((alt < ang + buffer) & (alt > 0))[0]
                elif not ang_type and not is_max:
                    regions = np.where(alt > ang - buffer)[0]
                elif ang_type and is_max:
                    regions = np.where((nadir < ang + buffer) & (alt > 0))[0]
                elif ang_type and not is_max:
                    regions = np.where((nadir > ang - buffer) & (alt > 0))[0]

                regions = np.split(regions, np.where(np.diff(regions) != 1)[0] + 1)

                self.gs[pos].encounters[enc].encounter_indices = regions

    def encounter_stats(self) -> pd.DataFrame:
        """Return encounter statistics.

        This method produces various statistics for each `EncounterSpec` object
        located within each `GroundPosition.encounters` dictionary of the
        `Encounter.gs` attribute. Statistics include the raw number of viable
        passes, cumulative encounter time, daily average counts, and the daily
        average time for each encounter and encounter type.

        Returns
        -------
        pd.DataFrame
            Pandas DataFrame containing the statistics for each encounter and
            encounter type.

        Notes
        -----
        The returned pandas DataFrame can be printed for easy viewing.

        Examples:
        ---------
        >>> stats = encounters.encounter_stats()

        See documentation for more detail on instantiated dependencies.
        """

        data = {}

        num_DL = 0
        time_DL = 0
        avg_num_DL = 0
        avg_time_DL = 0

        num_IMG = 0
        time_IMG = 0
        avg_num_IMG = 0
        avg_time_IMG = 0

        for pos in self.gs:
            for enc in self.gs[pos].encounters:

                enc_type = self.gs[pos].encounters[enc].type
                windows = self.gs[pos].encounters[enc].windows
                length = round(self.gs[pos].encounters[enc].length, 2)
                time = round(np.sum(windows[:, 2].astype(float)), 2)

                start = datetime.strptime(windows[0, 0], "%Y-%m-%d %H:%M:%S.%f")
                end = datetime.strptime(windows[-1, 1], "%Y-%m-%d %H:%M:%S.%f")
                deltaT = (end-start).days

                avg_num = round(length / deltaT, 2)
                avg_time = round(time / length, 2)

                if enc_type == "DL":
                    num_DL += length
                    time_DL += time
                    avg_num_DL += avg_num
                    avg_time_DL += avg_time

                if enc_type == "IMG":
                    num_IMG += length
                    time_IMG += time
                    avg_num_IMG += avg_num
                    avg_time_IMG += avg_time

                data[enc] = pd.Series([length, avg_num, time, avg_time])

        data["Total DL"] = pd.Series([num_DL, avg_num_DL, time_DL, avg_time_DL])
        data["Total IMG"] = pd.Series([num_IMG, avg_num_IMG, time_IMG, avg_time_IMG])
        df = pd.DataFrame.from_dict(data)
        df.index = ["Number of Viable Encounters",
                    "Average Encounters per Day",
                    "Viable Encounters Duration (s)",
                    "Average Encounter Duration (s)"]

        return df

    def save_windows(self, fileName: str, delimiter: Literal[",", "\\t"]) -> None:
        """Save window data to local directory.

        Parameters
        ----------
        fileName : str
            File name of output file. Can be a .txt or .csv file.
        delimiter : {",", "\\t"}
            String of length 1 representing the feild delimiter for the output
            file.

        Notes
        -----
        It is recommended to use a tab delimiter for .txt files and comma
        delimiters for .csv files. The method will return an error if fileName
        already exists in the current working directory. The windows
        function must be called prior to the save_windows method call.

        Examples
        --------
        >>> encounters.windows(finch)
        >>> encounters.save_windows("EncounterWindows.txt", "\\t")

        See documentation for more detail on instantiated dependencies.
        """
        encounter_types = np.array([])
        encounter_start = np.array([])
        encounter_end = np.array([])
        elapsed_sec = np.array([])
        max_alt = np.array([])
        min_nadir = np.array([])

        for pos in self.gs:
            for enc in self.gs[pos].encounters:

                enc = self.gs[pos].encounters[enc]
                encounter_types = np.append(encounter_types, np.full((enc.length,), enc))
                encounter_start = np.append(encounter_start, enc.windows[:, 0])
                encounter_end = np.append(encounter_end, enc.windows[:, 1])
                elapsed_sec = np.append(elapsed_sec, enc.windows[:, 2])
                max_alt = np.append(max_alt, enc.windows[:, 3])
                min_nadir = np.append(min_nadir, enc.windows[:, 4])

        data = {}
        data["Encounter Start"] = pd.Series(encounter_start)
        data["Encounter End"] = pd.Series(encounter_end)
        data["Elapsed Seconds"] = pd.Series(elapsed_sec)
        data["Maximum Altitude"] = pd.Series(max_alt)
        data["Minimum Nadir"] = pd.Series(min_nadir)
        data["Name"] = pd.Series(encounter_types)

        df = pd.DataFrame(data)
        df.to_csv(fileName, sep=delimiter)

        return df
