"""Localize ground position information."""


from celest.core.decorators import set_module
from celest.encounter import EncounterSpec
from typing import Tuple, Literal
import numpy as np


@set_module('celest.encounter')
class GroundPosition(object):
    """Localize ground position based information.

    The `GoundPosition` class stores ground location data and is used in the
    `Satellite` and `Encounter` class' to manage the groundposition-encounter
    relationship.

    Parameters
    ----------
    name : str
        Name of ground location used for identification and indexing.
    coor : tuple
        Geodetic coordinates for the ground position location given in decimal
        degrees in the (latitude, longitude) format.
    
    Attributes
    ----------
    name : str
        Name of ground location used for identification and indexing.
    coor : tuple
        Geodetic coordinates for the ground position location given in decimal
        degrees in the (latitude, longitude) format.
    radius : float
        Radius of earths surface at the given coordinates using WGS84.
    position : Coordinate
        Coordinate object containing the ground position coordinates.
    encounters : dict
        Dictionary where the keys are the encounter names and the values are
        the corresponding `EncounterSpec` objects.

    Methods
    -------
    _radius(obsCoor)
        Used to instantiate the radius attribute.
    add_encounter(name, encType, ang, angType, maxAng, solar)
        Define a ground position specific encounter.

    Examples
    --------
    >>> toronto = GroundPosition(name="Toronto", coor=(43.662300, -79.394530))
    """
    
    def __init__(self, name: str, coor: Tuple[float, float]) -> None:
        """Initialize instance variables."""

        self.name = name
        self.coor = coor
        self.radius = self._radius(coor)
        self.encounters = {}
    
    def __str__(self) -> str:
        """Defines GroundPosition information string."""

        title = 'Celest.GoundPosition Object\n'
        name = f'Name: {self.name}\t'
        coor = f'Coordinates: {self.coor}\t'
        radius = f'Radius: {self.radius}'

        return title + name + coor + radius
    
    def _radius(self, obsCoor: Tuple[float, float]) -> float:
        """Instantiates radius attribute.

        This method uses the World Geodetic System, WGS84, to calculate the
        Earth's radius at the given coordinates.

        Parameters
        ----------
        obsCoor : tuple
            Geodetic coordinates for the ground position location given in
            decimal degrees in the (latitude, longitude) format.

        Returns
        -------
        float
            The Earths radius at obsCoor given in km.

        Notes
        -----
        The Earth can be modeled as an ellipsoid given by the following
        equation:

        .. math:: r = \sqrt{\\frac{(6378.14)^2(6356.75)^2}{(6378.14)^2\sin{\phi}^2+(6356.75)^2\cos{\phi}^2}}

        where :math:`\phi` is the observers lattitude.
        """

        phi = np.radians(obsCoor[0])

        # Define WGS84 Parameters.
        semiMajor = 6378.137**2
        semiMinor = 6356.752314245**2

        numerator = semiMajor * semiMinor
        denominator = semiMajor * np.sin(phi)**2 + semiMinor * np.cos(phi)**2

        return np.sqrt(numerator / denominator)
    
    def add_encounter(self, name: str, encType: Literal["I", "T"], ang: float,
                      angType: Literal["A", "N"], maxAng: bool, solar:
                      Literal[-1, 0, 1]=0) -> None:
        """Define an encounter for the current ground position.

        Parameters
        ----------
        name : str
            String acting as the encounter identifier. Used for later indexing.
        encType :  {"I", "T"}
            Specifies encounter category as either imaging or transmission.
        ang : float
            Angluar constraint for the encounter in degrees.
        angType : {"A", "N"}
            String specifying the constraint angle as either the altitude, or
            off-nadir angle type.
        maxAng : bool
            Defines the contraint angle as a maximum constraint if True or as
            minimum constraint if False. Note that the off-nadir angle is
            measured to increase away from nadir.
        solar : {-1, 0, 1}, optional
            Defines sunlight constraint where -1 gets windows at night, 0 gets
            windows at day or night, and 1 gets windows at day.
        """

        encounter_object = EncounterSpec(name, encType, ang, angType, maxAng, solar)
        self.encounters[name] = encounter_object