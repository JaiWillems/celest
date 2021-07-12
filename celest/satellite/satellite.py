"""Satellite orbital representations and coordinate conversions.

The `Satellite` object localizes satellite related information and
functionality, and is passed into the `Encounter` class for encounter
planning.
"""


from celest.core.decorators import set_module
from celest.satellite import Coordinate
from typing import Literal
import numpy as np


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
    save_data(fileName, delimiter)
        Save the time and position data of the satellite.

    """
    
    def __init__(self, coordinates: Coordinate) -> None:
        """Initialize attributes."""

        self.time = Coordinate.timeData
        self.position = Coordinate
    
    def solar_power(self, pointProfiles: np.array=None) -> float:
        pass
    
    def solar_radiation_pressure(self, pointProfiles: np.array=None) -> float:
        pass
    
    def save_data(self, fileName: str, delimiter: Literal[",", "\t"]) -> None:
        pass

