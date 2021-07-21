"""Satellite encounter planning.

The encounter module contains the `Encounter` class that is used to compute,
store, and schedule Earth-satellite encounters.
"""


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
    gs : dict
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

    def windows(self, interp: bool=True, factor: int=5, buffer: float=10, dt: int=1) -> None:
        """Initialize the `GroundPosition.encounters.windows` attribute.

        This function iterates through all `gs` values (`GroundPosition`
        objects`) and the values in their `encounters` attributes
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

                groundPos = self.gs[pos]
                enc = self.gs[pos].encounters[enc]
                ang = enc.ang
                angType = enc.angType
                maxAng = enc.maxAng
                solar = enc.solar

                if interp:
                    encounter_ind = self.window_encounter_indices(buffer=buffer)
                    times = self._satellite.times
                    alt = self._satellite.position.horizontal(groundPos, factor=factor, dt=dt)[:, 0]
                    nadir = self._satellite.position.off_nadir(groundPos, factor=factor, dt=dt)
                else:
                    times = self._satellite.times
                    alt = self._satellite.position.horizontal(groundPos)[:, 0]
                    nadir = self._satellite.position.off_nadir(groundPos)

                if not angType and maxAng:
                    winInd = np.where((0 < alt) & (alt < ang))[0]
                elif not angType and not maxAng:
                    winInd = np.where(alt > ang)[0]
                elif angType and maxAng:
                    winInd = np.where((nadir < ang) & (alt >= 0))[0]
                elif angType and not maxAng:
                    winInd = np.where((nadir > ang) & (alt >= 0))[0]

                # Get sun position vector.
                if solar != 0:
                    sunECEFpos = Sun().position(timeData=times).ECEF()
                    gndECEFpos = np.concatenate(np.array(list(groundPos.coor), np.array(groundPos.radius)))
                    gndECEFpos = np.full((sunECEFpos.shape[0], 3), gndECEFpos)
                    dividend = np.einsum("ij, ij->i", sunECEFpos, gndECEFpos)
                    divisor = np.multiply(np.linalg.norm(sunECEFpos, axis=1), np.linalg.norm(gndECEFpos, axis=1))
                    arg = np.divide(dividend, divisor)
                    ang = np.degrees(np.arccos(arg))

                # Find intersection of night indices and window indices.
                if solar == -1:
                    nightInd = np.where(ang >= 90)
                    winInd = np.intersect1d(winInd, nightInd)

                # Find intersection of day indices and window indices.
                if solar == 1:
                    dayInd = np.where(ang < 90)
                    winInd = np.intersect1d(winInd, dayInd)

                winIndArr = np.split(winInd, np.where(np.diff(winInd) != 1)[0]+1)
                windowTimes = np.empty((len(winIndArr), 5), dtype="U25")

                for j in range(len(winIndArr)):

                    if j == 0:
                        break

                    start = julian.from_jd(times[winIndArr[j][0]])
                    end = julian.from_jd(times[winIndArr[j][-1]])
                    maxAlt = np.max(alt[winIndArr[j]])
                    minNadir = np.min(nadir[winIndArr[j]])

                    windowTimes[j, 0] = start
                    windowTimes[j, 1] = end
                    windowTimes[j, 2] = (end-start).total_seconds()
                    windowTimes[j, 3] = maxAlt
                    windowTimes[j, 4] = minNadir

                enc.windows = windowTimes
                enc.length = windowTimes.shape[0]
    
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
                angType = self.gs[pos].encounters[enc].angType
                ang = self.gs[pos].encounters[enc].ang
                isMax = self.gs[pos].encounters[enc].maxAng
                groundPos = self.gs[pos].encounters[enc].groundPos

                # Get ground position information.
                alt = groundPos.alt
                nadir = groundPos.nadirAng

                # Derive interpolation regions.
                if not angType and isMax:
                    regions = np.where((alt < ang + buffer) & (alt > 0))[0]
                elif not angType and not isMax:
                    regions = np.where(alt > ang - buffer)[0]
                elif angType and isMax:
                    regions = np.where((nadir < ang + buffer) & (alt > 0))[0]
                elif angType and not isMax:
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

        numDL = 0
        timeDL = 0
        avgNumDL = 0
        avgTimeDL = 0

        numIMG = 0
        timeIMG = 0
        avgNumIMG = 0
        avgTimeIMG = 0

        for pos in self.gs:
            for enc in self.gs[pos].encounters:

                encType = self.gs[pos].encounters[enc].type
                windows = self.gs[pos].encounters[enc].windows
                length = round(self.gs[pos].encounters[enc].length, 2)
                time = round(np.sum(windows[:, 2].astype(float)), 2)

                start = datetime.strptime(windows[0, 0], "%Y-%m-%d %H:%M:%S.%f")
                end = datetime.strptime(windows[-1, 1], "%Y-%m-%d %H:%M:%S.%f")
                deltaT = (end-start).days

                avgNum = round(length / deltaT, 2)
                avgTime = round(time / length, 2)

                if encType == "DL":
                    numDL += length
                    timeDL += time
                    avgNumDL += avgNum
                    avgTimeDL += avgTime

                if encType == "IMG":
                    numIMG += length
                    timeIMG += time
                    avgNumIMG += avgNum
                    avgTimeIMG += avgTime

                data[enc] = pd.Series([length, avgNum, time, avgTime])

        data["Total DL"] = pd.Series([numDL, avgNumDL, timeDL, avgTimeDL])
        data["Total IMG"] = pd.Series([numIMG, avgNumIMG, timeIMG, avgTimeIMG])
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
        encounterTypes = np.array([])
        encounterStart = np.array([])
        encounterEnd = np.array([])
        elapsedSec = np.array([])
        maxAlt = np.array([])
        minNadir = np.array([])

        for pos in self.gs:
            for enc in self.gs[pos].encounters:

                enc = self.gs[pos].encounters[enc]
                encounterTypes = np.append(encounterTypes, np.full((enc.length,), enc))
                encounterStart = np.append(encounterStart, enc.windows[:, 0])
                encounterEnd = np.append(encounterEnd, enc.windows[:, 1])
                elapsedSec = np.append(elapsedSec, enc.windows[:, 2])
                maxAlt = np.append(maxAlt, enc.windows[:, 3])
                minNadir = np.append(minNadir, enc.windows[:, 4])

        data = {}
        data["Encounter Start"] = pd.Series(encounterStart)
        data["Encounter End"] = pd.Series(encounterEnd)
        data["Elapsed Seconds"] = pd.Series(elapsedSec)
        data["Maximum Altitude"] = pd.Series(maxAlt)
        data["Minimum Nadir"] = pd.Series(minNadir)
        data["Name"] = pd.Series(encounterTypes)

        df = pd.DataFrame(data)
        df.to_csv(fileName, sep=delimiter)

        return df
