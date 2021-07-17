"""Satellite orbital representations and coordinate conversions.

The `Satellite` object localizes satellite related information and
functionality, and is passed into the `Encounter` class for encounter
planning.
"""


from celest.core.decorators import set_module
from celest.satellite import Coordinate
from typing import List, Literal
import pandas as pd


@set_module('celest.satellite')
class Satellite(object):
    """Localize satellite information and functionality.

    The `Satellite` class represents a satellite, be it artificial or natural,
    and allows for the position to be represented with time through multiple
    representations.

    Parameters
    ----------
    coordinates : Coordinate
        `Coordinate` object containing the position and time evolution of the
        satellite.

    Attributes
    ----------
    time : Time
        Times associated with the satellite positions.
    position : Coordinate
        Position of the satellite.

    Methods
    -------
    solar_power(pointProfiles)
        Calculate the solar power generated from the satellite solar cells.
    solar_radiation_pressure(pointProfiles)
        Calculate the solar radiation pressure experienced by the satellite.
    save_data(fileName, delimiter, posTypes)
        Save the time and position data of the satellite.
    """
    
    def __init__(self, coordinates: Coordinate) -> None:
        """Initialize attributes."""

        self.time = coordinates.timeData
        self.position = coordinates
    
    def save_data(self, fileName: str, delimiter: Literal[",", "\\t"], posTypes: List) -> None:
        """Save satellite data to local directory.

        Parameters
        ----------
        fileName : str
            File name of the output file as wither a .txt or .csv file.
        delimiter : str
            String of length 1 representing the feild delimiter for the output
            file.
        posTypes : List
            List containing the types of position data to store. Possible
            list values include "GEO", "ECI", "ECEF", and "Altitude".

        Notes
        -----
        It is recommended to use a tab delimiter for .txt files and comma
        delimiters for .csv files. The method will return an error if the
        fileName already exists in the current working directory.

        Examples
        --------
        >>> posTypes = ["ECI", "ECEF"]
        >>> finch.save_data(fileName="data.csv", delimiter=",", posTypes=posTypes)
        """

        data = {}
        data["Time (julian)"] = pd.Series(self.time)
        if "GEO" in posTypes:
            GEOpos = self.position.GEO()
            data["GEO.lat"] = pd.Series(GEOpos[:, 0])
            data["GEO.lon"] = pd.Series(GEOpos[:, 1])
            data["GEO.radius"] = pd.Series(GEOpos[:, 2])
        if "ECI" in posTypes:
            ECIpos = self.position.ECI()
            data["ECI.X"] = pd.Series(ECIpos[:, 0])
            data["ECI.y"] = pd.Series(ECIpos[:, 1])
            data["ECI.z"] = pd.Series(ECIpos[:, 2])
        if "ECEF" in posTypes:
            ECEFpos = self.position.ECEF()
            data["ECEF.X"] = pd.Series(ECEFpos[:, 0])
            data["ECEF.y"] = pd.Series(ECEFpos[:, 1])
            data["ECEF.z"] = pd.Series(ECEFpos[:, 2])
        if "Altitude" in posTypes:
            data["Altitude"] = pd.Series(self.position.altitude())

        df = pd.DataFrame(data)
        df.to_csv(fileName, sep=delimiter)
