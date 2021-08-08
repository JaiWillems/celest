"""Determine celestial object positional information.

The astronomy module contains the `CelestialObject` class to access positional
information.
"""


from celest.core.decorators import set_module
from jplephem import SPK
from typing import Literal
import numpy as np
import pkg_resources


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
    
    def _find_position(self, timeData: np.array, bcCode: int, pCode: int) -> np.array:
        """Calculate celestial object position.

        Parameters
        ----------
        timeData : np.array
            Array of shape (n,) containing Julian times where position
            calculations are desired.
        bcCode : int
            The ephem code for the objects barycenter.
        pCode : int
            The ephem code for the object.
        
        Returns
        -------
        np.array
            Array of shape (n, 3) containing the ECI position representation.
        """

        ssb2posb = self.kernal[0, bcCode].compute(timeData)
        posb2pos = self.kernal[bcCode, pCode].compute(timeData)
        ssb2eb = self.kernal[0, 3].compute(timeData)
        eb2e = self.kernal[3, 399].compute(timeData)
        e2pos = (ssb2posb + posb2pos - ssb2eb - eb2e).T

        return e2pos
    
    def _find_bc_position(self, timeData: np.array, bcCode: int) -> np.array:
        """Calculate celestial object barycenter position.

        Parameters
        ----------
        timeData : np.array
            Array of shape (n,) containing Julian times where position
            calculations are desired.
        bcCode : int
            The ephem code for the objects barycenter.
        
        Returns
        -------
        np.array
            Array of shape (n, 3) containing the ECI barycenter position
            representation.
        """

        ssb2posb = self.kernal[0, bcCode].compute(timeData)
        ssb2eb = self.kernal[0, 3].compute(timeData)
        eb2e = self.kernal[3, 399].compute(timeData)
        e2pos = (ssb2posb - ssb2eb - eb2e).T

        return e2pos
    
    def _find_altitude_zeros(self, posData: np.array, timeData: np.array, slope: Literal[-1, 1], shift: int=0) -> np.array:
        """Find times corresponding to zero altitude angles.

        Parameters
        ----------
        posData : np.array
            Array of shape (n,) containing altitude data of a space bound body
            from the ground location of interest.
        timeData : np.array
            Array of shape (n,) containing julian timeData corresponding to the
            inputted `posData`.
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
        np.array
            Array of shape (n,) containing julian times coresponding to zero
            altitude angles.
        """

        alt = posData - shift
        time = timeData

        ind_1 = np.where(np.diff(np.sign(alt)) == slope * 2)[0]
        ind_2 = ind_1 + 1

        m = (time[ind_2] - time[ind_1]) / (alt[ind_2] - alt[ind_1])
        rise_times = - m * alt[ind_1] + time[ind_1]

        return rise_times
    
    def _find_rise(self, posData: np.array, timeData: np.array) -> np.array:
        """Return rise times.

        Return the rise times of an object with the time evolving position,
        `posData` (with reference to a ground location). All rise times
        encapsulated in the given position data will be returned.

        Parameters
        ----------
        posData : np.array
            Array of shape (n,) containing altitude data of a rising object
            from the ground location of interest.
        timeData : np.array
            Array of shape (n,) containing julian timeData corresponding to the
            inputted `posData`.
        
        Returns
        -------
        np.array
            Array of shape (n,) containing julian times coresponding to the
            objects rise times as seen from the positions corresponding ground
            position.
        """

        return self._find_altitude_zeros(posData, timeData, slope=1)
    
    def _find_set(self, posData: np.array, timeData: np.array) -> np.array:
        """Return set times.

        Return the set times of an object with the time evolving position,
        `posData`(with reference to a ground location). All set times
        encapsulated in the given position data will be returned.

        Parameters
        ----------
        posData : np.array
            Array of shape (n,) containing altitude data of a setting object
            from the ground location of interest.
        timeData : np.array
            Array of shape (n,) containing julian timeData corresponding to the
            inputted `posData`.
        
        Returns
        -------
        np.array
            Array of shape (n,) containing julian times coresponding to the
            objects set times as seen from the positions corresponding ground
            position.
        """

        return self._find_altitude_zeros(posData, timeData, slope=-1)
    
    def _find_peak(self, posData: np.array, timeData: np.array) -> np.array:
        """Return peak times.

        Return the peak times of an object with the time evolving position,
        `posData` (with reference to a ground location). All peak times
        encapsulated in the given position data will be returned.

        Parameters
        ----------
        posData : np.array
            Array of shape (n,) containing altitude data of an orbiting object
            from the ground location of interest.
        timeData : np.array
            Array of shape (n,) containing julian timeData corresponding to the
            inputted `posData`.
        
        Returns
        -------
        np.array
            Array of shape (n,) containing julian times coresponding to the
            objects peak times as seen from the positions corresponding ground
            position.
        """

        ind = np.where(np.diff(np.sign(np.diff(posData))) == -2)[0] + 1
        peak_times = timeData[ind]

        return peak_times
