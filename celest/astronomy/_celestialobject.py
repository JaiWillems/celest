"""Determine celestial object positional information.

The astronomy module contains the `CelestialObject` class to access positional
information. All time and position representations are returned as `Time` and
`Coordinate` objects to provide the flexibility of different data
representations.
"""


from celest.core.decorators import set_module
from celest.satellite import Coordinate, Time
from celest.encounter import GroundPosition
from jplephem import SPK
from typing import Literal
import pkg_resources
import numpy as np


@set_module('celest.astronomy')
class CelestialObject(object):
    """Celestial object positional information class.

    The `CelestialObject` class is a parent class for the `Planets`, `Sun`, and
    `Moon` classes to provide certain general positional information.
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

        return Coordinate(basePos=e2pos, type="ECI", timeData=timeData)
    
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

        return Coordinate(basePos=e2pos, type="ECI", timeData=timeData)
    
    def _find_altitude_zeros(self, posData: Coordinate, groundPos: GroundPosition, slope: Literal[-1, 1], shift: int=0) -> Time:
        """Find times corresponding to zero altitude angles.

        Parameters
        ----------
        posData : Coordinate
            Position data of a space bound body.
        groundPos : GroundPosition
            Observer position.
        slope : {-1, 1}
            If `slope=-1`, find the zero positions where the derivative is
            negative. If `slope=1`, find the zero positions where the
            derivative is positive.
        shift : int, optional
            The inputted shift, is the desired "zero" position to find the
            altitude zeros, essentially acting as the vertical shift of the
            altitude vs time plot.
        
        Returns
        -------
        Time
            The times coresponding to zero altitude angles.
        """

        alt = posData.horizontal(groundPos)[:, 0] - shift
        time = posData.time.julian

        ind_1 = np.where(np.diff(np.sign(alt)) == slope * 2)[0]
        ind_2 = ind_1 + 1

        m = (time[ind_2] - time[ind_1]) / (alt[ind_2] - alt[ind_1])
        rise_times = - m * alt[ind_1] + time[ind_1]

        return Time(rise_times)
    
    def _find_rise(self, posData: Coordinate, groundPos: GroundPosition) -> Time:
        """Return rise times.

        Return the rise times of an object with the time evolving position,
        `posData`, from the ground location, `groundPos`. All rise times
        encapsulated in the given position data will be returned.

        Parameters
        ----------
        posData : Coordinate
            Position data of the rising object.
        groundPos : GroundPosition
            Ground location of observer.
        
        Returns
        -------
        Time
            The objects rise times as seen from `groundPos`.
        """

        return self._find_altitude_zeros(posData, groundPos, slope=1)
    
    def _find_set(self, posData: Coordinate, groundPos: GroundPosition) -> Time:
        """Return set times.

        Return the set times of an object with the time evolving position,
        `posData`, from the ground location, `groundPos`. All set times
        encapsulated in the given position data will be returned.

        Parameters
        ----------
        posData : Coordinate
            Position data of the rising object.
        groundPos : GroundPosition
            Ground location of observer.
        
        Returns
        -------
        Time
            The objects set times as seen from `groundPos`.
        """

        return self._find_altitude_zeros(posData, groundPos, slope=-1)
    
    def _find_peak(self, posData: Coordinate, groundPos: GroundPosition) -> Time:
        """Return peak times.

        Return the peak times of an object with the time evolving position,
        `posData`, from the ground location, `groundPos`. All peak times
        encapsulated in the given position data will be returned.

        Parameters
        ----------
        posData : Coordinate
            Position data of the orbiting object.
        groundPos : GroundPosition
            Ground location of observer.
        
        Returns
        -------
        Time
            The objects peak times as seen from `groundPos`.
        """
        alt = posData.horizontal(groundPos, factor=5)[:, 0]
        time = posData.time.julian

        ind = np.where(np.diff(np.sign(np.diff(alt))) == -2)[0] + 1

        peak_times = time[ind]

        return Time(peak_times)
