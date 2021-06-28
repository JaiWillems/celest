"""Determine the position of celestial objects.

The astronomy module was created to allow easy access to positions of various
celestial objects for all times from the present to 2050. It can be used
with the `Satellite` class to provide an observer centric coordinate system
such as the horizontal system for easy navigation of the night sky.
"""


import numpy as np
import pkg_resources
from jplephem.spk import SPK


class CelestialObject(object):
    """Determine ECI positions of various celestial objects.

    The `CelestialObject` class allows one to gain the ECI position data of
    various celestial objects using JPL's de421 ephemeris.

    Attributes
    ----------
    kernal : jplephem.daf.DAF
        Object containing the de421 emphemeris information.
    
    Methods
    -------
    sun_position(tiimeData)
        Get the Sun's position in the ECI frame.
    moon_position(tiimeData)
        Get the Moon's position in the ECI frame.
    mercury_position(tiimeData)
        Get Mercury's position in the ECI frame.
    venus_position(tiimeData)
        Get Venus's position in the ECI frame.
    mars_position(tiimeData)
        Get Mars's position in the ECI frame.
    """

    def __init__(self):
        """Initialize the ephemeris kernal."""
        ephem = pkg_resources.resource_filename(__name__, 'data/de421.bsp')
        self.kernal = SPK.open(ephem)

    def sun_position(self, timeData: np.array) -> list:
        """Calculate the Sun's position in the ECI frame.

        Parameters
        ----------
        timeData : np.array
            Array of shape (n,) containing Julian times as floats.
        
        Returns
        -------
        list
            Array of shape (n,3) with columns of X, Y, Z ECI position data.
        """
        # Get the sun's position in ECI.
        ssb2sun = self.kernal[0, 10].compute(timeData)
        ssb2eb = self.kernal[0, 3].compute(timeData)
        eb2e = self.kernal[3, 399].compute(timeData)
        e2sunECI = (ssb2sun - ssb2eb - eb2e).T

        return e2sunECI
    
    def moon_position(self, timeData: np.array) -> np.array:
        """Calculate the Moon's position in the ECI frame.

        Parameters
        ----------
        timeData : np.array
            Array of shape (n,) containing Julian times as floats.
        
        Returns
        -------
        list
            Array of shape (n,3) with columns of X, Y, Z ECI position data.
        """
        # Get the moon's position in ECI.
        eb2m = self.kernal[3, 301].compute(timeData)
        eb2e = self.kernal[3, 399].compute(timeData)
        e2mECI = (eb2m - eb2e).T

        return e2mECI

    def mercury_position(self, timeData: np.array) -> np.array:
        """Calculate Mercury's position in the ECI frame.

        Parameters
        ----------
        timeData : np.array
            Array of shape (n,) containing Julian times as floats.
        
        Returns
        -------
        list
            Array of shape (n,3) with columns of X, Y, Z ECI position data.
        """
        # Get mercury's position in ECI.
        ssb2mercb = self.kernal[0, 1].compute(timeData)
        mercb2merc = self.kernal[1, 199].compute(timeData)
        ssb2eb = self.kernal[0, 3].compute(timeData)
        eb2e = self.kernal[3, 399].compute(timeData)
        e2mercECI = (ssb2mercb + mercb2merc - ssb2eb - eb2e).T

        return e2mercECI

    def venus_position(self, timeData: np.array) -> np.array:
        """Calculate Venus's position in the ECI frame.

        Parameters
        ----------
        timeData : np.array
            Array of shape (n,) containing Julian times as floats.
        
        Returns
        -------
        list
            Array of shape (n,3) with columns of X, Y, Z ECI position data.
        """
        # Get venus's position in ECI.
        ssb2vb = self.kernal[0, 2].compute(timeData)
        vb2v = self.kernal[2, 299].compute(timeData)
        ssb2eb = self.kernal[0, 3].compute(timeData)
        eb2e = self.kernal[3, 399].compute(timeData)
        e2vECI = (ssb2vb + vb2v - ssb2eb - eb2e).T

        return e2vECI

    def mars_position(self, timeData: np.array) -> np.array:
        """Calculate Mar's position in the ECI frame.

        Parameters
        ----------
        timeData : np.array
            Array of shape (n,) containing Julian times as floats.
        
        Returns
        -------
        list
            Array of shape (n,3) with columns of X, Y, Z ECI position data.
        """
        # Get mars's position in ECI.
        ssb2marsb = self.kernal[0, 4].compute(timeData)
        marsb2mars = self.kernal[4, 499].compute(timeData)
        ssb2eb = self.kernal[0, 3].compute(timeData)
        eb2e = self.kernal[3, 399].compute(timeData)
        e2marsECI = (ssb2marsb + marsb2mars - ssb2eb - eb2e).T

        return e2marsECI
    
    def jupiter_position(self, timeData: np.array) -> np.array:
        """Calculate Jupiter's position in the ECI frame.

        Parameters
        ----------
        timeData : np.array
            Array of shape (n,) containing Julian times as floats.
        
        Returns
        -------
        list
            Array of shape (n,3) with columns of X, Y, Z ECI position data.
        """
        # Get jupiter's position in ECI.
        ssb2jb = self.kernal[0, 5].compute(timeData)
        ssb2eb = self.kernal[0, 3].compute(timeData)
        eb2e = self.kernal[3, 399].compute(timeData)
        e2jbECI = (ssb2jb - ssb2eb - eb2e).T

        return e2jbECI
    
    def saturn_position(self, timeData: np.array) -> np.array:
        """Calculate Saturn's position in the ECI frame.

        Parameters
        ----------
        timeData : np.array
            Array of shape (n,) containing Julian times as floats.
        
        Returns
        -------
        list
            Array of shape (n,3) with columns of X, Y, Z ECI position data.
        """
        # Get saturn's position in ECI.
        ssb2sb = self.kernal[0, 6].compute(timeData)
        ssb2eb = self.kernal[0, 3].compute(timeData)
        eb2e = self.kernal[3, 399].compute(timeData)
        e2sbECI = (ssb2sb - ssb2eb - eb2e).T

        return e2sbECI
    
    def uranus_position(self, timeData: np.array) -> np.array:
        """Calculate Uranus's position in the ECI frame.

        Parameters
        ----------
        timeData : np.array
            Array of shape (n,) containing Julian times as floats.
        
        Returns
        -------
        list
            Array of shape (n,3) with columns of X, Y, Z ECI position data.
        """
        # Get uranus's position in ECI.
        ssb2ub = self.kernal[0, 7].compute(timeData)
        ssb2eb = self.kernal[0, 3].compute(timeData)
        eb2e = self.kernal[3, 399].compute(timeData)
        e2ubECI = (ssb2ub - ssb2eb - eb2e).T

        return e2ubECI
    
    def neptune_position(self, timeData: np.array) -> np.array:
        """Calculate Neptune's position in the ECI frame.

        Parameters
        ----------
        timeData : np.array
            Array of shape (n,) containing Julian times as floats.
        
        Returns
        -------
        list
            Array of shape (n,3) with columns of X, Y, Z ECI position data.
        """
        # Get neptune's position in ECI.
        ssb2nb = self.kernal[0, 8].compute(timeData)
        ssb2eb = self.kernal[0, 3].compute(timeData)
        eb2e = self.kernal[3, 399].compute(timeData)
        e2nbECI = (ssb2nb - ssb2eb - eb2e).T

        return e2nbECI
    
    def pluto_position(self, timeData: np.array) -> np.array:
        """Calculate Pluto's position in the ECI frame.

        Parameters
        ----------
        timeData : np.array
            Array of shape (n,) containing Julian times as floats.
        
        Returns
        -------
        list
            Array of shape (n,3) with columns of X, Y, Z ECI position data.
        """
        # Get pluto's position in ECI.
        ssb2pb = self.kernal[0, 9].compute(timeData)
        ssb2eb = self.kernal[0, 3].compute(timeData)
        eb2e = self.kernal[3, 399].compute(timeData)
        e2pbECI = (ssb2pb - ssb2eb - eb2e).T

        return e2pbECI
