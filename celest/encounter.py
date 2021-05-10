"""Satellite encounter windows and scheduling.

Notes
-----
This module contains the Encounter class to determine and schedule
satellite-earth encounters.

Class's
-------
Encounter: Object used to calculate and store encounter information.
    addEncounter : Defines and stores EncounterSpec object in encounters
        attribute.
    getSunPos : Instantiates sunPos attribute.
    getWindows : Instantiates windows attribute of EncounterSpec object.
    saveWindows : Saves window data in local directory.

EncounterSpec : Object used to localize GoundPosition encounter information.
"""


import numpy as np
import pandas as pd
import pkg_resources
import julian
from datetime import datetime
from jplephem.spk import SPK
from celest import Satellite



class Encounter:
    """
    Encounter object computes, schedules, and stores satellite encounters.

    Methods
    -------
    addEncounter : Defines and stores EncounterSpec object in encounters
        attribute.
    getSunPos : Instantiates sunPos attribute.
    getWindows : Instantiates windows attribute of EncounterSpec object.
    saveWindows : Saves window data in local directory.

    Instance Variables
    ------------------
    encounters : Dictionary of EncounterSpec objects.
    sunPos : Sun positions in ECEF frame.
    """

    def __init__(self):
        """Define instance variables."""

        self.encounters = {}
        self.sunPos = None

    def addEncounter(self, name, encType, groundPos, ang, angType, maxAng, solar=0):
        """
        Define an encounter type.

        Uses the input data to create a key/value pair in the dictionary of the
        encounters attribute.

        Parameters
        ----------
        name : str
            String acting as the encounter identifier.
        encType :  str
            Specifies encounter as either imaging, "IMG", or data linking,
            "DL".
        goundPos : GroundPosition object
            Ground location associated with the encounter.
        ang : int or float
            Angluar constraint for the encounter.
        angType : str
            String specifying the constraint angle as either the altitude,
            "alt", or nadir-LOS, "nadirLOS", angle type.
        maxAng : bool
            Defines the contraint angle as a maximum constraint if True or as
            minimum constraint if False. Note that the nadirLOS angle is
            measured to increase away from nadir.
        solar : int
            Defines sunlight constraint. -1 -> windows at night, 0 -> windows
            at day or night, 1 -> windows at day.

        Returns
        -------
        None
        """
        if angType == "nadirLOS":
            angType = 1
        elif angType == "alt":
            angType = 0

        encounter = EncounterSpec(name, encType, groundPos, ang, angType, maxAng, solar)
        self.encounters[name] = encounter

    def getSunPos(self, timeData):
        """
        Instantiate sunPos attribute.

        Uses de421 ephemeris to calculate the sun"s position in the ECEF frame
        at each time in timeData.

        Parameters
        ----------
        timeData : ndarray
            Array of shape (n,). Contains datetime objects in UTC.

        Returns
        -------
        self.sunPos : list
            Array of shape (n,3) with columns of X, Y, Z ECEF sun position
            data.
        """
        # Get sun position in ECI.
        julianTimes = np.zeros((timeData.size,))

        for i, time in enumerate(timeData):
            julianTimes[i] = julian.to_jd(datetime.strptime(time[:18],
                                          "%Y-%m-%d %H:%M:%S"))

        ephem = pkg_resources.resource_filename(__name__, 'data/de421.bsp')
        kernal = SPK.open(ephem)

        ssb2sun = kernal[0, 10].compute(julianTimes)
        ssb2eb = kernal[0, 3].compute(julianTimes)
        eb2e = kernal[3, 399].compute(julianTimes)
        e2sun = (ssb2sun - ssb2eb - eb2e).T

        # Convert sun position data to ECEF
        sun = Satellite()
        e2sunECEF = sun.getECEF(posData=e2sun, timeData=timeData)
        self.sunPos = e2sunECEF
        return self.sunPos

    def getWindows(self, satellite):
        """
        Instantiates windows attribute of EncounterSpec objects.

        Determines the windows for all encounter types specified in the
        encounters attribute.

        Parameters
        ----------
        satellite : Satellite object
            Satellite object with appropriate dependencies instantiated.

        Returns
        -------
        self.windows : dict
            Dictionary of key/value pairs where the keys are the encounter
            names and the values are the windows associated with a specific
            encounter type.

        Usage
        -----
        Satellite.times attribute must be instantiated. nadirAng must
        be instantiated for the ground position located in satellite.gs if the
        corresponding encounter is a "nadirLOS" type. alt and az must be
        instantiated for the ground position located in satellite.gs if the
        encounter is an "alt" type.
        """

        for i in self.encounters:

            encounter = self.encounters[i]
            groundPos = encounter.groundPos.name
            ang = encounter.ang
            angType = encounter.angType
            maxAng = encounter.maxAng
            solar = encounter.solar

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
                sunECEFpos = self.getSunPos(timeData=times)
                gndECEFpos = np.full((sunECEFpos.shape[0], 3),
                                     satellite.gs[groundPos].ECEFpos)
                dividend = np.einsum("ij, ij->i", sunECEFpos, gndECEFpos)
                divisor = np.multiply(np.linalg.norm(sunECEFpos, axis=1),
                                      np.linalg.norm(gndECEFpos, axis=1))
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

                start = datetime.strptime(times[winIndArr[j][0]][:18],
                                          "%Y-%m-%d %H:%M:%S")
                end = datetime.strptime(times[winIndArr[j][-1]][:18],
                                        "%Y-%m-%d %H:%M:%S")

                windowTimes[j, 0] = start
                windowTimes[j, 1] = end
                windowTimes[j, 2] = (end-start).total_seconds()

            self.encounters[i].windows = windowTimes
            self.encounters[i].length = windowTimes.shape[0]

    def saveWindows(self, fileName, delimiter):
        """
        Save window data to local directory.

        Parameters
        ----------
        fileName : str
            File name of output file. Can be a .txt or .csv file.
        delimiter : str
            String of length 1 representing the feild delimiter for the output
            file.

        Returns
        -------
        None

        Usage
        -----
        It is recommended to use a tab delimiter for .txt files and comma
        delimiters for .csv files. Will return an error if fileName already
        exists in the current working directory.
        """
        encounterTypes = np.array([])
        encounterStart = np.array([])
        encounterEnd = np.array([])
        elapsedSec = np.array([])

        for key in self.encounters:
            encounterTypes = np.append(encounterTypes,
                                       np.full((self.encounters[key].length,),
                                               key))
            encounterStart = np.append(encounterStart,
                                       self.encounters[key].windows[:, 0])
            encounterEnd = np.append(encounterEnd,
                                     self.encounters[key].windows[:, 1])
            elapsedSec = np.append(elapsedSec,
                                   self.encounters[key].windows[:, 2])

        data = {}
        data["Encounter Start"] = pd.Series(encounterStart)
        data["Encounter End"] = pd.Series(encounterEnd)
        data["Elapsed Seconds"] = pd.Series(elapsedSec)
        data["Name"] = pd.Series(encounterTypes)

        df = pd.DataFrame(data)
        df.to_csv(fileName, sep=delimiter)

    def getStats(self):
        """
        Get encounter statistics.

        This method produces various statistics for each encounter and
        encounter type. The generated statistics include the raw number of
        viable passes, cumulative time, daily average counts, and the daily
        average time for each encounter and encounter type.

        Parameters
        ----------
        None

        Returns
        -------
        data : Pandas DataFrame
            Pandas DataFrame containing the statistics for each encounter and
            encounter type.

        Usage
        -----
        The returned pandas DataFrame can be printed for easy viewing of the
        statistics.
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

            start = datetime.strptime(windows[0, 0], "%Y-%m-%d %H:%M:%S")
            end = datetime.strptime(windows[-1, 1], "%Y-%m-%d %H:%M:%S")
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


class EncounterSpec:
    """
    EncounterSpec object localizes encounter information specific to a
    GroundPosition Object. Units are in degrees and meters.

    Methods
    -------
    None

    Instance Variables
    ------------------
    name : Encounter identifier.
    type : Encounter type.
    groundPos : GroundPosition object related to the encounter.
    ang : Constraint angle for a viable encounter.
    angType : Constraint angle type.
    maxAng : Direction of the angular constraint.
    solar : Encounter daylight requirement.
    windows : Encounter Windows.
    length : Length of data attributes.
    """

    def __init__(self, name, encType, groundPos, ang, angType, maxAng, solar=0):
        """
        Define instance variables.

        Parameters
        ----------
        name : str
            String acting as the encounter identifier.
        encType :  str
            Specifies encounter as either imaging, "IMG", or data linking,
            "DL".
        goundPos : GroundPosition object
            Ground location associated with the encounter.
        ang : int or float
            Angluar constraint for the encounter.
        angType : str
            String specifying the constraint angle as either the altitude,
            "alt", or nadir-LOS, "nadirLOS", angle type.
        maxAng : bool
            Defines the contraint angle as a maximum constraint if True or as
            minimum constraint if False. Note that the nadirLOS angle is
            measured to increase away from nadir.
        solar : int
            Defines sunlight constraint. -1 -> windows at night, 0 -> windows
            at day or night, 1 -> windows at day.
        """

        self.name = name
        self.type = encType
        self.groundPos = groundPos
        self.ang = ang
        self.angType = angType
        self.maxAng = maxAng
        self.solar = solar
        self.windows = None
        self.length = 0

    def __str__(self):
        """Defines EncounterSpec informaiton string."""
        data = np.array([self.name, self.type, self.groundPos.name, self.ang,
                         self.angType, self.maxAng, self.solar])
        index = np.array(["Name:", "Encounter Type:", "Ground Position Name:",
                          "Constraint Angle:", "Angle Type:", "Is Max:",
                          "Solar Constraint:"])

        df = pd.DataFrame(data=data, index=index)
        return df.to_string()
