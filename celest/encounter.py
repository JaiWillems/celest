"""Satellite encounter planning.

The encounter module contains the Encounter class that is used to compute,
store, and schedule Earth-satellite encounters.
"""


import numpy as np
from scipy.interpolate import interp1d
import pandas as pd
import pkg_resources
import julian
from datetime import datetime
from jplephem.spk import SPK
from typing import Literal, Tuple
from celest.satellite import Satellite
from celest.groundposition import GroundPosition


class EncounterSpec(object):
    """Localize encounter information.

    The EncounterSpec object localizes encounter information specific to a
    GroundPosition Object.

    Parameters
    ----------
    name : str
        String acting as the encounter identifier. Used for later indexing.
    encType :  {"IMG", "DL"}
        Specifies encounter category as either imaging or data linking.
    goundPos : GroundPosition
        Ground location associated with the encounter.
    ang : float
        Angluar constraint for the encounter in degrees.
    angType : {"alt", "nadirLOS"}
        String specifying the constraint angle as either the altitude, or
        nadir-LOS angle type.
    maxAng : bool
        Defines the contraint angle as a maximum constraint if True or as
        minimum constraint if False. Note that the nadirLOS angle is
        measured to increase away from nadir.
    solar : {-1, 0, 1}, optional
        Defines sunlight constraint where -1 gets windows at night, 0 gets
        windows at day or night, and 1 gets windows at day.

    Attributes
    ----------
    name : str
        Encounter identifier.
    type : {"IMG", "DL"}
        Specifies encounter category as either imaging or data linking.
    groundPos : GroundPosition
        Ground location associated with the encounter.
    ang : float
        Angluar constraint for the encounter in degrees.
    angType : {"alt", "nadirLOS"}
        String specifying the constraint angle as either the altitude, or
        nadir-LOS angle type.
    maxAng : bool
        Defines the contraint angle as a maximum constraint if True or as
        minimum constraint if False. Note that the nadirLOS angle is
        measured to increase away from nadir.
    solar : {-1, 0, 1}, optional
        Defines sunlight constraint where -1 gets windows at night, 0 gets
        windows at day or night, and 1 gets windows at day.
    windows : np.array
        Array of shape (n,3) of window start, end, and elapsed seconds data.
    length : int
        Length, n, of data attributes.
    """

    def __init__(self, name: str, encType: Literal["IMG", "DL"], groundPos:
                 GroundPosition, ang: float, angType: Literal["alt", "nadirLOS"],
                 maxAng: bool, solar: Literal[-1, 0, 1]=0) -> None:
        """Define instance variables."""
        self.name = name
        self.type = encType
        self.groundPos = groundPos
        self.ang = ang
        self.angType = angType
        self.maxAng = maxAng
        self.solar = solar
        self.windows = None
        self.length = 0

    def __str__(self) -> str:
        """Defines EncounterSpec informaiton string."""
        data = np.array([self.name, self.type, self.groundPos.name, self.ang,
                         self.angType, self.maxAng, self.solar])
        index = np.array(["Name:", "Encounter Type:", "Ground Position Name:",
                          "Constraint Angle:", "Angle Type:", "Is Max:",
                          "Solar Constraint:"])

        df = pd.DataFrame(data=data, index=index)
        return df.to_string()


class Encounter(object):
    """Computes, schedules, and store satellite encounter data.

    The Encounter class extends off the information stored in Satellite and
    GroundPosition objects to compute, store, and schedule Earth-Satellite
    encounters.

    Attributes
    ----------
    encounters : Dict
        Dictionary of EncounterSpec objects where the key is the
        EncounterSpec's name attribute.
    sunPos : np.array
        Array of shape (n,3) containing the sun position in the ECEF
        cartesian frame.

    Methods
    -------
    add_encounter(name, encType, groundPos, ang, angType, maxAng, solar=0)
        Defines and stores EncounterSpec object in the encounters attribute
        dictionary.
    _sun_position(timeData)
        Instantiates sunPos attribute.
    windows(satellite)
        Instantiates windows attribute of the EncounterSpec objects.
    save_windows(fileName, delimiter)
        Saves window data in local directory.
    encounter_stats()
        Get encounter statistics.
    """

    def __init__(self) -> None:
        """Define instance variables."""
        self.encounters = {}
        self.sunPos = None

    def add_encounter(self, name: str, encType: Literal["IMG", "DL"],
                      groundPos: GroundPosition, ang: float, angType:
                      Literal["alt", "nadirLOS"], maxAng: bool, solar:
                      Literal[-1, 0, 1]=0) -> None:
        """Define an encounter type.

        This method uses the input data to create a key/value pair in the
        encounters attribute dictionary.

        Parameters
        ----------
        name : str
            String acting as the encounter identifier. Used for later indexing.
        encType :  {"IMG", "DL"}
            Specifies encounter category as either imaging or data linking.
        goundPos : GroundPosition
            Ground location associated with the encounter.
        ang : float
            Angluar constraint for the encounter in degrees.
        angType : {"alt", "nadirLOS"}
            String specifying the constraint angle as either the altitude, or
            nadir-LOS angle type.
        maxAng : bool
            Defines the contraint angle as a maximum constraint if True or as
            minimum constraint if False. Note that the nadirLOS angle is
            measured to increase away from nadir.
        solar : {-1, 0, 1}, optional
            Defines sunlight constraint where -1 gets windows at night, 0 gets
            windows at day or night, and 1 gets windows at day.

        Examples
        --------
        >>> encounters = Encounter()
        >>> encounters.add_encounter("CYYZ IMG", "IMG", toronto, 30, "nadirLOS",
        ...                          True, solar=1)
        """
        if angType == "nadirLOS":
            angType = 1
        elif angType == "alt":
            angType = 0

        encounter = EncounterSpec(
            name, encType, groundPos, ang, angType, maxAng, solar)
        self.encounters[name] = encounter

    def _sun_position(self, timeData: np.array) -> list:
        """Instantiate sunPos attribute.

        This method uses the de421 ephemeris to calculate the sun"s position in
        the ECEF cartesian frame at each time in the timeData parameter.

        Parameters
        ----------
        timeData : np.array
            Array of shape (n,) containing datetime objects in UTC.

        Returns
        -------
        list
            Array of shape (n,3) with columns of X, Y, Z ECEF sun position
            data.

        Examples
        --------
        >>> encounter = Encounter()
        >>> UTCTimeData = np.array(['2020-06-01 12:00:00.0340', ...,
        ...                        '2020-06-01 12:01:00.0340'])
        >>> sunPos = encounter._sun_position(timeData=UTCTimeData)
        """
        ephem = pkg_resources.resource_filename(__name__, 'data/de421.bsp')
        kernal = SPK.open(ephem)

        # Get sun position in ECI.
        ssb2sun = kernal[0, 10].compute(timeData)
        ssb2eb = kernal[0, 3].compute(timeData)
        eb2e = kernal[3, 399].compute(timeData)
        e2sun = (ssb2sun - ssb2eb - eb2e).T

        # Convert sun position data to ECEF
        sun = Satellite()
        e2sunECEF = sun.ECEF(posData=e2sun, timeData=timeData, jul=True)
        self.sunPos = e2sunECEF
        return self.sunPos

    def _special_interp(self, satellite: Satellite, groundPos: GroundPosition,
                        encounter: EncounterSpec, factor: int, buffer: float,
                        dt: int) -> Tuple[np.array, np.array, np.array, np.array]:
        """Interpolate encounter regions.

        This method identifies the positional data relevent to the inputed
        encounter and interpolates these regions.

        Parameters
        ----------
        satellite : Satellite
            Satellite object related to the encounters.
        grouundPos : GroundPosition
            Specifies ground location associated with the encounters.
        encounter : EncounterSpec
            Defines the valid encounter region.
        factor : int
            The factor increase in the number of steps in interpolated regions.
        buffer : float
            Positive angular buffer to interpolate data that comes within the
            buffer from the encounter constraint angle.
        dt : int
            The number of data points adjacent to a valid encounter region to
            interpolate within.

        returns
        -------
        np.array
            The array of new times corresponding to interpolated position data.
        np.array
            The array of interpolated altitude data.
        np.array
            The array of interpolated azimuth data.
        np.array
            The array of interpolated nadir angle data.

        Notes
        -----
        The inspiration for the buffer is to interpolate regions that have the
        potential to surpass the encounter constraint angle but may perhaps lie
        between data points. The buffer is recommended all times but
        escpecially with coarse data.
        """
        # Get satellite information.
        times = satellite.times
        ECEFdata = satellite.ECEFdata

        # Get encounter information.
        angType = encounter.angType
        ang = encounter.ang
        isMax = encounter.maxAng

        # Get ground position information.
        alt = groundPos.alt
        az = groundPos.az
        nadir = groundPos.nadirAng

        # Set up ECEF interpolation functions.
        x = interp1d(times, ECEFdata[:, 0], kind="cubic")
        y = interp1d(times, ECEFdata[:, 1], kind="cubic")
        z = interp1d(times, ECEFdata[:, 2], kind="cubic")

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

        for i in np.flip(regions):

            # Reform the viable region.
            i = np.append(i, np.arange(i[-1] + 1, i[-1] + 1 + dt, 1))
            i = np.insert(i, 0, np.arange(i[0] - dt, i[0], 1))

            minI = i[0]
            maxI = i[-1]

            # Determine new times and insert into current time data.
            time = np.linspace(times[minI], times[maxI],
                               int(factor * (maxI - minI + 1)))
            times = np.delete(times, i)
            times = np.insert(times, minI, time)

            # Interpolate, insert, and reform ECEF data.
            xInterp = np.delete(ECEFdata[:, 0], i)
            yInterp = np.delete(ECEFdata[:, 1], i)
            zInterp = np.delete(ECEFdata[:, 2], i)
            xN = np.insert(xInterp, minI, x(time)).reshape((-1, 1))
            yN = np.insert(yInterp, minI, y(time)).reshape((-1, 1))
            zN = np.insert(zInterp, minI, z(time)).reshape((-1, 1))
            ECEFdata = np.concatenate((xN, yN, zN), axis=1)

            # Get interpolated angle data.
            groundLoc = GroundPosition(groundPos.name, groundPos.coor)
            finch = Satellite()
            finch.time_data(timeData=times, jul=True)
            finch.position_data(posData=ECEFdata, type="ECEF")
            alt, az = finch.horizontal(groundPos=groundLoc)
            nadir = finch.nadir_ang(groundPos=groundLoc)

        return times, alt, az, nadir

    def windows(self, satellite: Satellite, interp: bool=True, factor: int=5,
                buffer: float=10, dt: int=1) -> None:
        """Instantiates windows attribute of EncounterSpec objects.

        This method determines determines the windows for each EncounterSpec
        object specified in the encounters attribute.

        Parameters
        ----------
        satellite : Satellite
            Satellite object with appropriate dependencies instantiated.
        interp : {True, False}, optional
            If True then the windows will be generated from data interpolated
            around valid encounter regions.
        factor : int, optional
            The factor increase in the number of steps in interpolated regions.
        buffer : float, optional
            Positive angular buffer to interpolate data that comes within the
            buffer from the encounter constraint angle.
        dt : int, optional
            The number of data points adjacent to a valid encounter region to
            interpolate within.

        Notes
        -----
        If an EncounterSpec objects use an "alt" angle type then the Satellite
        must have the corresponding GroundPositions alt and az attributes
        instantiated. If an EncounterSpec objects use a "nadirLOS" angle type
        then the Satellite must have the corresponding GroundPositions nadirAng
        attribute instantiated.

        Exampls:
        --------
        >>> encounters.windows(finch)

        See documentation for more detail on instantiated dependencies.
        """

        for i in self.encounters:

            encounter = self.encounters[i]
            groundPos = encounter.groundPos.name
            ang = encounter.ang
            angType = encounter.angType
            maxAng = encounter.maxAng
            solar = encounter.solar

            if interp:
                params = [satellite, satellite.gs[groundPos],
                          encounter, factor, buffer, dt]
                times, alt, az, nadir = self._special_interp(*params)
            else:
                times = satellite.times
                alt = satellite.gs[groundPos].alt
                nadir = satellite.gs[groundPos].nadirAng

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
                sunECEFpos = self._sun_position(timeData=times)
                gndECEFpos = np.full(
                    (sunECEFpos.shape[0], 3), satellite.gs[groundPos].ECEFpos)
                dividend = np.einsum("ij, ij->i", sunECEFpos, gndECEFpos)
                divisor = np.multiply(np.linalg.norm(
                    sunECEFpos, axis=1), np.linalg.norm(gndECEFpos, axis=1))
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
            windowTimes = np.empty((len(winIndArr), 3), dtype="U25")

            for j in range(len(winIndArr)):

                if j == 0:
                    break

                start = julian.from_jd(times[winIndArr[j][0]])
                end = julian.from_jd(times[winIndArr[j][-1]])

                windowTimes[j, 0] = start
                windowTimes[j, 1] = end
                windowTimes[j, 2] = (end-start).total_seconds()

            self.encounters[i].windows = windowTimes
            self.encounters[i].length = windowTimes.shape[0]

    def save_windows(self, fileName: str, delimiter: str) -> None:
        """Save window data to local directory.

        Parameters
        ----------
        fileName : str
            File name of output file. Can be a .txt or .csv file.
        delimiter : str
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
        >>> encounters.save_windows("EncounterWindows.txt", "\t")

        See documentation for more detail on instantiated dependencies.
        """
        encounterTypes = np.array([])
        encounterStart = np.array([])
        encounterEnd = np.array([])
        elapsedSec = np.array([])

        for key in self.encounters:
            encounterTypes = np.append(encounterTypes, np.full(
                (self.encounters[key].length,), key))
            encounterStart = np.append(
                encounterStart, self.encounters[key].windows[:, 0])
            encounterEnd = np.append(
                encounterEnd, self.encounters[key].windows[:, 1])
            elapsedSec = np.append(
                elapsedSec, self.encounters[key].windows[:, 2])

        data = {}
        data["Encounter Start"] = pd.Series(encounterStart)
        data["Encounter End"] = pd.Series(encounterEnd)
        data["Elapsed Seconds"] = pd.Series(elapsedSec)
        data["Name"] = pd.Series(encounterTypes)

        df = pd.DataFrame(data)
        df.to_csv(fileName, sep=delimiter)

    def encounter_stats(self) -> pd.DataFrame:
        """Generate encounter statistics.

        This method produces various statistics for each encounterSpec object
        Statistics include the raw number of viable passes, cumulative
        encounter time, daily average counts, and the daily average time for
        each encounter and encounter type.

        Parameters
        ----------
        None

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

        for encounter in self.encounters:

            encType = self.encounters[encounter].type
            windows = self.encounters[encounter].windows
            length = round(self.encounters[encounter].length, 2)
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

            data[encounter] = pd.Series([length, avgNum, time, avgTime])

        data["Total DL"] = pd.Series([numDL, avgNumDL, timeDL, avgTimeDL])
        data["Total IMG"] = pd.Series([numIMG, avgNumIMG, timeIMG, avgTimeIMG])
        df = pd.DataFrame.from_dict(data)
        df.index = ["Number of Viable Encounters",
                    "Average Encounters per Day",
                    "Viable Encounters Duration (s)",
                    "Average Encounter Duration (s)"]

        return df
