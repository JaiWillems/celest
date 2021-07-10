"""
"""

from celest.core.decorators import set_module
from celest.satellite import Coordinate, Time
from celest.astronomy import CelestialObject
from celest.encounter import GroundPosition
import numpy as np


@set_module('celest.astronomy')
class Moon(CelestialObject):
    """
    """
    
    def __init__(self):
        super().__init__()

        self._pos_data = None

        pass
    
    def position(self, timeData: Time) -> Coordinate:
        pass
    
    def rise(self, timeData: Time, groundPos: GroundPosition) -> np.array:
        return self._find_rise(posData, timeData, groundPos)
    
    def set(self, timeData: Time, groundPos: GroundPosition) -> np.array:
        return self._find_set(posData, timeData, groundPos)
    
    def peak(self, timeData: Time, groundPos: GroundPosition) -> np.array:
        return self._find_peak(posData, timeData, groundPos)
    
    def phase(self, timeData: Time, groundPos: GroundPosition) -> np.array:
        pass
    
    def full_moon(self, timeData: Time, groundPos: GroundPosition) -> np.array:
        pass
    
    def lunar_eclipse(self, timeData: Time, groundPos: GroundPosition) -> np.array:
        pass
