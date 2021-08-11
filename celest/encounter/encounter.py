"""Satellite pass analysis and encounter planning.

The encounter module contains the `Encounter` class that is used to compute,
store, and schedule Earth-satellite encounters.
"""


from celest.core.decorators import set_module
from celest.encounter import analytical_encounter_ind
from celest.satellite import Coordinate, Time
from datetime import datetime
from typing import Any, Dict, Literal
from jplephem.spk import SPK
import julian
import numpy as np
import pkg_resources
import pandas as pd


@set_module('celest.encounter')
class Encounter(object):
    """Compute, schedule, and store satellite encounter data.

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
        Add a ground position for encounter calculations.
    windows(interp=True, factor=5, buffer=10, dt=1)
        Initialize the `GroundPosition.encounters.windows` attribute.
    window_encounter_ind(groundPos, ang, angType, sca=0, solar=0, **kwargs)
        Return window encounter indices.
    encounter_stats()
        Return encounter statistics.
    save_windows(fileName, delimiter)
        Save window data to local directory.
    """

    def __init__(self, satellite: Any) -> None:
        """Initialize attribites."""

        self._satellite = satellite
        self._sun_position = None
        self.gs = {}

    def add_position(self, groundPos: Any) -> None:
        """Add a ground position for encounter calculations.

        This function adds `groundPos` to the `gs` attribute dictionary where
        the key is `groundPos.name` and the value is `groundPos` itself. The
        `gs` attribute will be iterated through when determining encounter
        windows.

        Parameters
        ----------
        groundPos : GroundPosition
            `GroundPosition` object instantiated as per its documentation.
        """

        self.gs[groundPos.name] = groundPos

    def _sun_ECEF(self, timeData: np.array) -> np.array:
        """Return Sun ECEF position.

        Parameters
        ----------
        timeData : np.array
            Array of shape (n,) containing satellite Julian time data.

        Returns
        -------
        np.array
            Array of shape (n, 3) containing Sun ECEF positions.
        """

        ephem = pkg_resources.resource_filename(__name__, '../data/de421.bsp')
        self.kernal = SPK.open(ephem)

        ssb2sunb = self.kernal[0, 10].compute(timeData)
        ssb2eb = self.kernal[0, 3].compute(timeData)
        eb2e = self.kernal[3, 399].compute(timeData)
        e2sun = (ssb2sunb - ssb2eb - eb2e).T

        sun_ECEF = Coordinate(e2sun, "ECI", Time(timeData)).ECEF()

        return sun_ECEF

    def _get_ang(self, vecOne: np.array, vecTwo: np.array) -> float:
            """Calculate degree angle bewteen two vectors.

            Parameters
            ----------
            vecOne, vecTwo : np.array
                Arrays of shape (n,3) with rows of ECEF data.

            Returns
            -------
            float
                Degree angle between the two arrays.
            """

            # Use simple linalg formula.
            dividend = np.sum(vecOne * vecTwo, axis=1)
            norm_one = np.linalg.norm(vecOne, axis=1)
            norm_two = np.linalg.norm(vecTwo, axis=1)
            divisor = norm_one * norm_two
            arg = np.divide(dividend, divisor)
            ang = np.degrees(np.arccos(arg))

            return ang

    def window_encounter_ind(self, groundPos: Any, ang: float, angType:
                             Literal["A", "T"], sca: float=0, solar:
                             Literal[-1, 0, 1]=0, **kwargs:
                             Dict[str, Any]) -> np.array:
        """Return window encounter indices.

        This method calculates the window encounter indices for the encounter
        defined by the `ang`, `angType`, `sca`, and `solar` parameters with
        the `groundPos` location.

        Parameters
        ----------
        groundPos : GroundPosition
            Specifying the ground location encounter.
        ang : float
            Angluar constraint for the encounter in degrees.
        angType : {"A", "N"}
            String specifying the constraint angle as either the altitude, or
            off-nadir angle type. Refer to notes for greater detail in the
            angle definitions.
        SCA : float, optional
            Float specifying the minimum angle between a satellite's position
            vector and the Sun's position vector in a ground-location-centric
            reference system. Refer to notes for more information.
        solar : {-1, 0, 1}, optional
            Defines sunlight constraint where -1 gets windows at night, 0 gets
            windows at day or night, and 1 gets windows at day.
        kwargs : Dict, optional
            Optional keyword arguments passed into the `Interpolation`
            classes `_interp()` method.

        Returns
        -------
        np.array
            Array of arrays of indices defining encounter windows.

        Notes
        -----
        The altitude angle is defined as the angle of a celestial object
        measured in increasing degrees from the horizon. A valid encounter
        region with the altitude angle are the regions where `alt_angle > ang`.
        The off-nadir angle is defined as the angle of an Earth based location
        measured in increasing degrees from a satellite's nadir. A valid
        encounter region with the off-nadir angle are the regions where
        `nadir_ang < ang`.

        In some instances of ground to satellite communication, hardware damage
        can be incurred when the ground station is within a certain angle of
        the sun. The solar constraint angle allows encounters to be calculated
        out of direct alignment of the sun by invalidating encounter regions
        where the sun is behind or close to the satellite as seen from the
        ground station assuming the ground station is actively tracking the
        satellite.
        """

        sat_ECEF = self._satellite.position.ECEF(**kwargs)
        time_data = self._satellite.times.julian(**kwargs)

        ground_GEO = [groundPos.coor[0], groundPos.coor[1], groundPos.radius]
        ground_GEO = np.repeat(np.array([ground_GEO]), time_data.size, axis=0)
        gnd_ECEF = Coordinate(ground_GEO, "GEO", Time(time_data)).ECEF()

        enc_params = [sat_ECEF, gnd_ECEF, ang, angType]
        analytical_ind = analytical_encounter_ind(*enc_params)

        sun_ECEF = self._sun_ECEF(time_data)
        solar_angs = self._get_ang(sat_ECEF - gnd_ECEF, sun_ECEF)
        sca_ind = np.where(solar_angs > sca)[0]

        enc_ind = np.intersect1d(analytical_ind, sca_ind)

        if solar != 0:
            sun_angs = self._get_ang(gnd_ECEF, sun_ECEF)

            if solar == 1:
                # Isolate indices occuring at day.
                day_ind = np.where(sun_angs < 90)
                enc_ind = np.intersect1d(enc_ind, day_ind)
            else:
                # Isolate indices occuring at night.
                night_ind = np.where(sun_angs >= 90)
                enc_ind = np.intersect1d(enc_ind, night_ind)

        window_ind = np.split(enc_ind, np.where(np.diff(enc_ind) != 1)[0] + 1)

        return np.array(window_ind, dtype=object)

    def windows(self, interp: bool=True, factor: int=5, buffer: float=10, dt:
                int=1) -> None:
        """Initialize the `GroundPosition.encounters.windows` attribute.

        This function iterates through all encounters stored in groundposition
        objects stored within the satellite to calculate the windows associated
        with the encounters.

        Parameters
        ----------
        interp : bool, optional
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
                encounter = ground_pos.encounters[enc]
                ang_type = encounter.ang_type
                ang = encounter.ang
                sca = encounter.solar_constraint_ang
                solar = encounter.solar

                if interp:

                    if ang_type == "N":
                        off_nadir = self._satellite.position.off_nadir(ground_pos)
                        ind = np.where(off_nadir < ang + buffer)[0]

                    elif ang_type == "A":
                        alt, _ = self._satellite.position.horizontal(ground_pos)
                        ind = np.where(alt > ang - buffer)[0]

                    ind = np.split(ind, np.where(np.diff(ind) != 1)[0] + 1)

                    times = self._satellite.times.julian(factor=factor, dt=dt, indices=ind)
                    kwargs = {"factor": factor, "dt": dt, "indices": ind}

                else:
                    times = self._satellite.times.julian()
                    kwargs = {}

                enc_params = [ground_pos, ang, ang_type, sca, solar]
                window_ind = self.window_encounter_ind(*enc_params, **kwargs)

                if window_ind.size != 0:

                    n = len(window_ind)
                    windows = np.empty((n, 3), dtype="U25")

                    for j in range(n):

                        start = julian.from_jd(times[window_ind[j][0]])
                        end = julian.from_jd(times[window_ind[j][-1]])

                        windows[j, 0] = start
                        windows[j, 1] = end
                        windows[j, 2] = (end - start).total_seconds()

                    encounter.windows = windows
                    encounter.length = n

    def encounter_stats(self) -> pd.DataFrame:
        """Return encounter statistics.

        This method produces various statistics for each defined encounter.
        Statistics include the raw number of viable passes, cumulative
        encounter time, daily average counts, and the 5th, 50th, and 95th
        percentile encounter durations for each encounter and encounter type.

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
        window_len_DL = []

        num_IMG = 0
        time_IMG = 0
        avg_num_IMG = 0
        avg_time_IMG = 0
        window_len_IMG = []

        for pos in self.gs:
            for enc in self.gs[pos].encounters:

                enc_type = self.gs[pos].encounters[enc].type
                windows = self.gs[pos].encounters[enc].windows

                if windows is not None:
                    length = round(self.gs[pos].encounters[enc].length, 2)
                    time = round(np.sum(windows[:, 2].astype(float)), 2)
                    Q5 = np.percentile(windows[:, 2].astype(float), 5)
                    Q95 = np.percentile(windows[:, 2].astype(float), 95)

                    start = datetime.strptime(windows[0, 0], "%Y-%m-%d %H:%M:%S.%f")
                    end = datetime.strptime(windows[-1, 1], "%Y-%m-%d %H:%M:%S.%f")
                    deltaT = (end-start).total_seconds() / 3600 / 24

                    avg_num = round(length / deltaT, 2)
                    avg_time = round(time / length, 2)

                    if enc_type == "T":
                        num_DL += length
                        time_DL += time
                        avg_num_DL += avg_num
                        avg_time_DL += avg_time

                        window_len_DL.extend(windows[:, 2].astype(float).tolist())

                    if enc_type == "I":
                        num_IMG += length
                        time_IMG += time
                        avg_num_IMG += avg_num
                        avg_time_IMG += avg_time

                        window_len_IMG.extend(windows[:, 2].astype(float).tolist())

                    data[enc] = pd.Series([length, avg_num, time, Q5, avg_time, Q95])

        if len(window_len_DL) == 0:
            window_len_DL.append(0)
        if len(window_len_IMG) == 0:
            window_len_IMG.append(0)
        DL_Q5 = np.percentile(window_len_DL, 5)
        DL_Q95 = np.percentile(window_len_DL, 95)
        IMG_Q5 = np.percentile(window_len_IMG, 5)
        IMG_Q95 = np.percentile(window_len_IMG, 95)

        DL_params = [num_DL, avg_num_DL, time_DL, DL_Q5, avg_time_DL, DL_Q95]
        data["Total DL"] = pd.Series(DL_params)

        IMG_params = [num_IMG, avg_num_IMG, time_IMG, IMG_Q5, avg_time_IMG, IMG_Q95]
        data["Total IMG"] = pd.Series(IMG_params)

        df = pd.DataFrame.from_dict(data)
        df.index = ["Number of Viable Encounters",
                    "Average Encounters per Day",
                    "Viable Encounters Duration (s)",
                    "5th Percentile Encounter Duration (s)",
                    "50th Percentile Encounter Duration (s)",
                    "95th Percentile Encounter Duration (s)"]

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
        delimiters for .csv files. The windows function must be called prior
        to the save_windows method call.

        Examples
        --------
        >>> encounters.windows(finch)
        >>> encounters.save_windows("EncounterWindows.txt", "\\t")

        See documentation for more detail on instantiated dependencies.
        """

        encounter_names = np.array([])
        encounter_types = np.array([])
        encounter_start = np.array([])
        encounter_end = np.array([])
        elapsed_sec = np.array([])

        for pos in self.gs:
            for enc in self.gs[pos].encounters:

                enc = self.gs[pos].encounters[enc]

                if enc.windows is not None:
                    encounter_names = np.append(encounter_names,
                                                np.full((enc.length,),
                                                enc.name))
                    encounter_types = np.append(encounter_types,
                                                np.full((enc.length,),
                                                enc.type))
                    encounter_start = np.append(encounter_start,
                                                enc.windows[:, 0])
                    encounter_end = np.append(encounter_end,
                                              enc.windows[:, 1])
                    elapsed_sec = np.append(elapsed_sec, enc.windows[:, 2])

        data = {}
        data["Encounter Start"] = pd.Series(encounter_start)
        data["Encounter End"] = pd.Series(encounter_end)
        data["Elapsed Seconds"] = pd.Series(elapsed_sec)
        data["Type"] = pd.Series(encounter_types)
        data["Name"] = pd.Series(encounter_names)

        df = pd.DataFrame(data)
        df.to_csv(fileName, sep=delimiter)

        return df
