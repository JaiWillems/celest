"""
"""


from celest.core.decorators import set_module
from celest.satellite import Coordinate, Time
from celest.encounter import GroundPosition
from jplephem import SPK
import pkg_resources
import numpy as np


@set_module('celest.astronomy')
class CelestialObject(object):
    """
    """

    def __init__(self) -> None:
        """Initialize attributes."""
        
        ephem = pkg_resources.resource_filename(__name__, 'data/de421.bsp')
        self.kernal = SPK.open(ephem)
    
    def _find_position(self, timeData: Time, bcCode: int, pCode: int) -> Coordinate:
        """Calculate celestial object position.

        Parameters
        ----------
        timeData : Time
            Times where position calculations are desired.
        bcCode : int
            The ephem code for the objects barycenter.
        pCode : int
            The ephem code for the object.
        
        Returns
        -------
        Coordinate
            Object containing the position representations.
        """

        ssb2posb = self.kernal[0, bcCode].compute(timeData.julian)
        posb2pos = self.kernal[bcCode, pCode].compute(timeData.julian)
        ssb2eb = self.kernal[0, 3].compute(timeData.julian)
        eb2e = self.kernal[3, 399].compute(timeData.julian)
        e2pos = (ssb2posb + posb2pos - ssb2eb - eb2e).T

        return Coordinate(basePos=e2pos, type="GCRS", timeData=timeData)
    
    def _find_bc_position(self, timeData: Time, bcCode: int) -> Coordinate:
        """Calculate celestial object barycenter position.

        Parameters
        ----------
        timeData : Time
            Times where position calculations are desired.
        bcCode : int
            The ephem code for the objects barycenter.
        
        Returns
        -------
        Coordinate
            Object containing the barycenter position representations.
        """

        ssb2posb = self.kernal[0, bcCode].compute(timeData.julian)
        ssb2eb = self.kernal[0, 3].compute(timeData.julian)
        eb2e = self.kernal[3, 399].compute(timeData.julian)
        e2pos = (ssb2posb - ssb2eb - eb2e).T

        return Coordinate(basePos=e2pos, type="GCRS", timeData=timeData)
    
    def _find_rise(self, posData: Coordinate, timeData: Time, groundPos: GroundPosition) -> np.array:
        pass
    
    def _find_set(self, posData: Coordinate, timeData: Time, groundPos: GroundPosition) -> np.array:
        pass
    
    def _find_peak(self, posData: Coordinate, timeData: Time, groundPos: GroundPosition) -> np.array:
        pass