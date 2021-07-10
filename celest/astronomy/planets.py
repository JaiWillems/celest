"""Determine the position of celestial objects.

The astronomy module was created to allow easy access to positions of various
celestial objects for all times from the present to 2050. It can be used
with the `Satellite` class to provide an observer centric coordinate system
such as the horizontal system for easy navigation of the night sky.
"""


from celest.core.decorators import set_module
from celest.satellite import Coordinate, Time
from celest.astronomy import CelestialObject
from celest.encounter import GroundPosition
import numpy as np


@set_module('celest.astronomy')
class Planets(CelestialObject):
    """Determine positions of various celestial objects.

    The `Planets` class allows one to gain the position data of
    various celestial objects using JPL's de421 ephemeris.

    Attributes
    ----------
    kernal : jplephem.daf.DAF
        Object containing the de421 emphemeris information.
    
    Methods
    -------
    mercury_position(timeData)
        Get Mercury's position in the ECI frame.
    venus_position(timeData)
        Get Venus' position in the ECI frame.
    mars_position(timeData)
        Get Mars' position in the ECI frame.
    jupiter_position(timeData)
        Get Jupiter's position in the ECI frame.
    saturn_position(timeData)
        Get Saturn's position in the ECI frame.
    uranus_position(timeData)
        Get Uranus' position in the ECI frame.
    neptunr_position(timeData)
        Get Neptune's position in the ECI frame.
    pluto_position(timeData)
        Get Pluto's position in the ECI frame.
    """

    def __init__(self):
        super().__init__()

        self._mercury_pos_data = None
        self._venus_pos_data = None
        self._mars_pos_data = None
        self._jupiter_pos_data = None
        self._saturn_pos_data = None
        self._uranus_pos_data = None
        self._neptun_pos_data = None
        self._pluto_pos_data = None

        pass
    
    def mercury_position(self, timeData: Time) -> Coordinate:
        return self._find_position(bcCode=1, pCode=199)
    
    def venus_position(self, timeData: Time) -> Coordinate:
        return self._find_position(bcCode=2, pCode=299)
    
    def mars_position(self, timeData: Time) -> Coordinate:
        return self._find_position(bcCode=4, pCode=499)
    
    def jupiter_position(self, timeData: Time) -> Coordinate:
        return self._find_bc_position(bcCode=5)
    
    def saturn_position(self, timeData: Time) -> Coordinate:
        return self._find_bc_position(bcCode=6)
    
    def uranus_position(self, timeData: Time) -> Coordinate:
        return self._find_bc_position(bcCode=7)
    
    def neptun_position(self, timeData: Time) -> Coordinate:
        return self._find_bc_position(bcCode=8)
    
    def pluto_position(self, timeData: Time) -> Coordinate:
        return self._find_bc_position(bcCode=9)

    def rise(self, object: str, groundPos: GroundPosition) -> np.array:
        return self._find_rise(posData, timeData, groundPos)
    
    def set(self, object: str, groundPos: GroundPosition) -> np.array:
        return self._find_set(posData, timeData, groundPos)
    
    def peak(self, object: str, groundPos: GroundPosition) -> np.array:
        return self._find_peak(posData, timeData, groundPos)
