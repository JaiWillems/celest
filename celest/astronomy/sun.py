"""Determine Sun positional and special event information.

The astronomy module contains the `Sun` class to access positional and special
event information. All time and positional representations are returned as
`Time` and `Coordinate` objects to provide the flexibility of different data
representaitons.
"""


from celest.astronomy import CelestialObject
from celest.core.decorators import set_module
from celest.encounter import GroundPosition
from celest.satellite import Coordinate, Time
from typing import Literal


@set_module('celest.astronomy')
class Sun(CelestialObject):
    """Determine Sun positional and special event information.
    
    Methods
    -------
    position(timeData)
        Return the Sun's position.
    dawn_times(timeData, type, groundPos)
        Return dawn times.
    dusk_times(timeData, type, groundPos)
        Return dusk times.
    rise(timeData, groundPos)
        Return the solar rise times.
    set(timeData, groundPos)
        Return the solar set times.
    peak(timeData, groundPos)
        Return the solar peak times.
    """

    def __init__(self) -> None:
        """Initialize attributes."""

        super().__init__()
    
    def position(self, timeData: Time) -> Coordinate:
        """Return Sun's position.

        Parameters
        ----------
        timeData : Time
            Find the positions corresponding to the imputted times.
        
        Returns
        -------
        Coordinate
            Position of Mercury at times, `timeData`, as a `Coordinate` object.
        """
        
        return self._find_bc_position(timeData, bcCode=10)

    def dawn_times(self, timeData: Time, type: Literal["A", "C", "N"], groundPos: GroundPosition) -> Time:
        """Return dawn times.

        Return the civil, nautical, or astonomical dawn times encapsulated
        within the `timeData` parameter.

        Parameters
        ----------
        timeData : Time
            Time data encapsulating the desired dawn times.
        type : {"A", "C", "N"}
            String specifying the dawn times corresponding to civil, nautical,
            or astronomical dawn using "C", "N", and "A", respectively.
        groundPos : GroundPosition
            The ground location of the observer.
        
        Returns
        -------
        Time
            The dawn times as seen from `groundPos`.
        
        See Also
        --------
        dusk_times : Return dusk times.
        """

        shift = {"C": -6, "N": -12, "A": -18}
        sunPos = self.position(timeData)

        return self._find_altitude_zeros(sunPos, groundPos, slope=1, shift=shift[type])
    
    def dusk_times(self, timeData: Time, type: Literal["A", "C", "N"], groundPos: GroundPosition) -> Time:
        """Return dusk times.

        Return the civil, nautical, or astonomical dusk times encapsulated
        within the `timeData` parameter.

        Parameters
        ----------
        timeData : Time
            Time data encapsulating the desired dusk times.
        type : {"A", "C", "N"}
            String specifying the dusk times corresponding to civil, nautical,
            or astronomical dusk using "C", "N", and "A", respectively.
        groundPos : GroundPosition
            The ground location of the observer.
        
        Returns
        -------
        Time
            The dusk times as seen from `groundPos`.
        
        See Also
        --------
        dawn_times : Return dawn times.
        """
        
        shift = {"C": -6, "N": -12, "A": -18}
        sunPos = self.position(timeData)

        return self._find_altitude_zeros(sunPos, groundPos, slope=-1, shift=shift[type])

    
    def rise(self, timeData: Time, groundPos: GroundPosition) -> Time:
        """Return the solar rise times.

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
        
        See Also
        --------
        set : Return the solar set times.
        """
        
        posData = self.position(timeData)
        return self._find_rise(posData, groundPos)
    
    def set(self, timeData: Time, groundPos: GroundPosition) -> Time:
        """Return the solar set times.

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
        
        See Also
        --------
        rise : Return the solar rise times.
        """
        
        posData = self.position(timeData)
        return self._find_set(posData, groundPos)
    
    def peak(self, timeData: Time, groundPos: GroundPosition) -> Time:
        """Return the solar peak times.

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
        return self._find_peak(posData, groundPos)
