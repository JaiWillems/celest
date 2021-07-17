"""Determine the position of celestial objects.

The astronomy module contains the `Planets` class to access position
information for local planets. All time and positional representations are
returned as `Time` and `Coordinate` objects to provide the flexibility of
different data representaitons.
"""


from celest.astronomy import CelestialObject
from celest.core.decorators import set_module
from celest.encounter import GroundPosition
from celest.satellite import Coordinate, Time
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
        Return Mercury's position.
    venus_position(timeData)
        Return Venus' position.
    mars_position(timeData)
        Return Mars' position.
    jupiter_position(timeData)
        Return Jupiter's barycenter position.
    saturn_position(timeData)
        Return Saturn's barycenter position.
    uranus_position(timeData)
        Return Uranus's barycenter position.
    neptune_position(timeData)
        Return Neptune's barycenter position.
    pluto_position(timeData)
        Return Pluto's barycenter position.
    rise(object, timeData, groundPos)
        Return the rise times for a celestial object.
    set(object, timeData, groundPos)
        Return the set times for a celestial object.
    peak(object, timeData, groundPos)
        Return the peak times for a celestial object.
    """

    def __init__(self):
        """Initialize attributes."""

        super().__init__()
    
    def mercury_position(self, timeData: Time) -> Coordinate:
        """Return Mercury's position.

        Parameters
        ----------
        timeData : Time
            Find the positions corresponding to the imputted times.
        
        Returns
        -------
        Coordinate
            Position of Mercury at times, `timeData`, as a `Coordinate` object.
        """

        return self._find_position(timeData, bcCode=1, pCode=199)
    
    def venus_position(self, timeData: Time) -> Coordinate:
        """Return Venus' position.

        Parameters
        ----------
        timeData : Time
            Find the positions corresponding to the imputted times.
        
        Returns
        -------
        Coordinate
            Position of Mercury at times, `timeData`, as a `Coordinate` object.
        """

        return self._find_position(timeData, bcCode=2, pCode=299)
    
    def mars_position(self, timeData: Time) -> Coordinate:
        """Return Mars' position.

        Parameters
        ----------
        timeData : Time
            Find the positions corresponding to the imputted times.
        
        Returns
        -------
        Coordinate
            Position of Mercury at times, `timeData`, as a `Coordinate` object.
        """

        return self._find_position(timeData, bcCode=4, pCode=499)
    
    def jupiter_position(self, timeData: Time) -> Coordinate:
        """Return Jupiter's barycenter position.

        Parameters
        ----------
        timeData : Time
            Find the positions corresponding to the imputted times.
        
        Returns
        -------
        Coordinate
            Position of Mercury at times, `timeData`, as a `Coordinate` object.
        """

        return self._find_bc_position(timeData, bcCode=5)
    
    def saturn_position(self, timeData: Time) -> Coordinate:
        """Return Saturn's barycenter position.

        Parameters
        ----------
        timeData : Time
            Find the positions corresponding to the imputted times.
        
        Returns
        -------
        Coordinate
            Position of Mercury at times, `timeData`, as a `Coordinate` object.
        """

        return self._find_bc_position(timeData, bcCode=6)
    
    def uranus_position(self, timeData: Time) -> Coordinate:
        """Return Uranus's barycenter position.

        Parameters
        ----------
        timeData : Time
            Find the positions corresponding to the imputted times.
        
        Returns
        -------
        Coordinate
            Position of Mercury at times, `timeData`, as a `Coordinate` object.
        """

        return self._find_bc_position(timeData, bcCode=7)
    
    def neptune_position(self, timeData: Time) -> Coordinate:
        """Return Neptune's barycenter position.

        Parameters
        ----------
        timeData : Time
            Find the positions corresponding to the imputted times.
        
        Returns
        -------
        Coordinate
            Position of Mercury at times, `timeData`, as a `Coordinate` object.
        """

        return self._find_bc_position(timeData, bcCode=8)
    
    def pluto_position(self, timeData: Time) -> Coordinate:
        """Return Pluto's barycenter position.

        Parameters
        ----------
        timeData : Time
            Find the positions corresponding to the imputted times.
        
        Returns
        -------
        Coordinate
            Position of Mercury at times, `timeData`, as a `Coordinate` object.
        """

        return self._find_bc_position(timeData, bcCode=9)
    
    def _get_pos_data(self, object: int, timeData: Time) -> np.array:
        """Return the object's position data.

        Parameters
        ----------
        object : int
            The planets are defined by an integer identifier given by the
            following mapping: (0) Mercury, (1) Venus, (2) Mars, (3) Jupiter,
            (4) Saturn, (5) Uranus, (6) Neptune, and (7) Pluto.
        timeData : Time
            Find the positions corresponding to the imputted times.
        
        Returns
        -------
        Coordinate
            Position of Mercury at times, `timeData`, as a `Coordinate` object.
        """
        
        planet = {
            0: self.mercury_position,
            1: self.venus_position,
            2: self.mars_position,
            3: self.jupiter_position,
            4: self.saturn_position,
            5: self.uranus_position,
            6: self.neptune_position,
            7: self.pluto_position
        }

        return planet[object](timeData)

    def rise(self, object: int, timeData: Time, groundPos: GroundPosition) -> Time:
        """Return the rise times for a celestial object.

        Parameters
        ----------
        object : int
            The planets are defined by an integer identifier given by the
            following mapping: (0) Mercury, (1) Venus, (2) Mars, (3) Jupiter,
            (4) Saturn, (5) Uranus, (6) Neptune, and (7) Pluto.
        timeData : Time
            Find the rise times encapsulated within the imputted times.
        groundPos : GroundPosition
            Position of the observer.
        
        Returns
        -------
        Time
            `Time` object containing rise times.
        """

        posData = self._get_pos_data(self, object, timeData)
        return self._find_rise(posData, groundPos)
    
    def set(self, object: int, timeData: Time, groundPos: GroundPosition) -> Time:
        """Return the set times for a celestial object.

        Parameters
        ----------
        object : int
            The planets are defined by an integer identifier given by the
            following mapping: (0) Mercury, (1) Venus, (2) Mars, (3) Jupiter,
            (4) Saturn, (5) Uranus, (6) Neptune, and (7) Pluto.
        timeData : Time
            Find the set times encapsulated within the imputted times.
        groundPos : GroundPosition
            Position of the observer.
        
        Returns
        -------
        Time
            `Time` object containing set times.
        """

        posData = self._get_pos_data(self, object, timeData)
        return self._find_set(posData, groundPos)
    
    def peak(self, object: int, timeData: Time, groundPos: GroundPosition) -> Time:
        """Return the peak times for a celestial object.

        Parameters
        ----------
        object : int
            The planets are defined by an integer identifier given by the
            following mapping: (0) Mercury, (1) Venus, (2) Mars, (3) Jupiter,
            (4) Saturn, (5) Uranus, (6) Neptune, and (7) Pluto.
        timeData : Time
            Find the peak times encapsulated within the imputted times.
        groundPos : GroundPosition
            Position of the observer.
        
        Returns
        -------
        Time
            `Time` object containing peak times.
        """

        posData = self._get_pos_data(self, object, timeData)
        return self._find_peak(posData, groundPos)
