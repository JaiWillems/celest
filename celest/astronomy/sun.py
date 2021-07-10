"""
"""


from celest.core.decorators import set_module
from celest.satellite import Coordinate, Time
from celest.astronomy import CelestialObject
from celest.encounter import GroundPosition
import numpy as np


@set_module('celest.astronomy')
class Sun(CelestialObject):

    def __init__(self) -> None:
        super().__init__()

        self._pos_data = None

        pass
    
    def position(self, timeData: Time) -> Coordinate:
        return self._find_bc_position(bcCode=10)
    
    def mean_position(self, timeData: Time) -> Coordinate:
        pass
    
    def twilight_times(self, timeData: Time, type: str, groundPos: GroundPosition) -> np.array:
        """
        Refer to the following:
        https://www.timeanddate.com/astronomy/different-types-twilight.html
        """
        pass
    
    def rise(self, timeData: Time, groundPos: GroundPosition) -> np.array:
        return self._find_rise(posData, timeData, groundPos)
    
    def set(self, timeData: Time, groundPos: GroundPosition) -> np.array:
        return self._find_set(posData, timeData, groundPos)
    
    def peak(self, timeData: Time, groundPos: GroundPosition) -> np.array:
        return self._find_peak(posData, timeData, groundPos)

    def solar_eclipse(self, timeData: Time, groundPos: GroundPosition) -> np.array:
        pass
    
    def solstice(self, timeData: Time) -> np.array:
        pass
    
    def equinox(self, timeData: Time) -> np.array:
        pass
