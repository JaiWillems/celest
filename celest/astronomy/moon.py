"""Determine Moon positional and special event information.

The astronomy module contains the `Moon` class to access positional and special
event information. All time and positional representations are returned as
`Time` and `Coordinate` objects to provide the flexibility of different data
representaitons.
"""

from celest.core.decorators import set_module
from celest.satellite import Coordinate, Time
from celest.astronomy import CelestialObject
from celest.encounter import GroundPosition
import numpy as np


@set_module('celest.astronomy')
class Moon(CelestialObject):
    """Determine Moon positional and special event information.
    
    Methods
    -------
    position(timeData)
        Return the Moon's position.
    rise(timeData, groundPos)
        Return the lunar rise times.
    set(timeData, groundPos)
        Return the lunar set times.
    peak(timeData, groundPos)
        Return the lunar peak times.
    """
    
    def __init__(self):
        """Initialize attributes."""

        super().__init__()
    
    def position(self, timeData: Time) -> Coordinate:
        """Return the Moon's position.

        Parameters
        ----------
        timeData : Time
            Find the positions corresponding to the imputted times.
        
        Returns
        -------
        Coordinate
            Position of Mercury at times, `timeData`, as a `Coordinate` object.
        """

        eb2m = self.kernal[3, 301].compute(timeData.julian)
        eb2e = self.kernal[3, 399].compute(timeData.julian)
        e2m = (eb2m - eb2e).T

        return Coordinate(basePos=e2m, type="ECI", timeData=timeData)
    
    def rise(self, timeData: Time, groundPos: GroundPosition) -> np.array:
        """Return the lunar rise times.

        Parameters
        ----------
        timeData : Time
            Find the rise times encapsulated within the imputted times.
        groundPos : GroundPosition
            Position of the observer.
        
        Returns
        -------
        Time
            `Time` object containing rise times.
        """
        
        posData = self.position(timeData)
        return self._find_rise(posData, timeData, groundPos)
    
    def set(self, timeData: Time, groundPos: GroundPosition) -> np.array:
        """Return the lunar set times.

        Parameters
        ----------
        timeData : Time
            Find the set times encapsulated within the imputted times.
        groundPos : GroundPosition
            Position of the observer.
        
        Returns
        -------
        Time
            `Time` object containing set times.
        """
        
        posData = self.position(timeData)
        return self._find_set(posData, timeData, groundPos)
    
    def peak(self, timeData: Time, groundPos: GroundPosition) -> np.array:
        """Return the lunar peak times.

        Parameters
        ----------
        timeData : Time
            Find the peak times encapsulated within the imputted times.
        groundPos : GroundPosition
            Position of the observer.
        
        Returns
        -------
        Time
            `Time` object containing peak times.
        """
        
        posData = self.position(timeData)
        return self._find_peak(posData, timeData, groundPos)
