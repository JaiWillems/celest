"""Localize ground position information.

This module contains the `GoundPosition` class which stores ground location
data and is an instrumental part of encounter planning by managing the ground
location and satellite encounter relationship.
"""


from celest.core.decorators import set_module
from celest.encounter import EncounterSpec
from typing import Literal, Tuple
import numpy as np


@set_module('celest.encounter')
class GroundPosition(object):
    """Localize Earth bound location information.

    The `GoundPosition` class stores ground location data and is an
    instrumental part of encounter planning by managing the ground location and
    satellite encounter relationship.

    Parameters
    ----------
    name : str
        Name of the ground location used for identification and indexing.
    coor : tuple
        Geodetic coordinates for the ground location given in the
        (latitude, longitude) format where the values are in degrees and
        decimals.

    Attributes
    ----------
    name : str
        Name of the ground location used for identification and indexing.
    coor : tuple
        Geodetic coordinates for the location given in degrees and decimals.
    radius : float
        Earth's radius at the given coordinates.
    encounters : Dict
        Dictionary containing encounter information where the keys are the
        encounter names and the values are the corresponding `EncounterSpec`
        objects defining an encounter.

    Methods
    -------
    _radius(obsCoor)
        Calculate the Earth's geocentric radius using WGS84.
    add_encounter(name, encType, ang, angType, maxAng, solar)
        Define encounter for the ground location.
    """

    def __init__(self, name: str, coor: Tuple[float, float]) -> None:
        """Initialize attributes."""

        self.name = name
        self.coor = coor
        self.radius = self._radius(coor[0])
        self.encounters = {}

    def __str__(self) -> str:
        """Define information string."""

        output = [
            f"name={self.name}",
            f"coor={self.coor}",
            f"radius={round(self.radius, 5)}",
            f"{len(self.encounters)} encounter(s)"
        ]

        return ", ".join(output)

    def _radius(self, latitude: float) -> float:
        """Calculate the Earth's geocentric radius using WGS84.

        Parameters
        ----------
        latitude : float
            Latitude of the ground location in degrees and decimals.

        Returns
        -------
        float
            Earth's geocentric radius in kilometres and decimals.

        Notes
        -----
        By using an Earth ellipsoid with the WGS84 parameters of
        :math:`a=6378.137` and :math:`b=6356.7523142`, the geocentric radius
        can be calculated using the following formulation:

        .. math:: r = \sqrt{\frac{(a^2\cos(\beta))^2 + (b^2\sin(\beta))^2}{(a\cos(\beta))^2 + (b\sin(\beta))^2}}

        where :math:`\beta` is the observer's latitude.[1]_

        References
        ----------
        .. [1] Timur. Earth Radius by Latitude (WGS 84). 2018. url:
           https://planetcalc.com/7721/.
        """

        # Get lattidue parameter.
        phi = np.radians(latitude)

        # Define WGS84 Parameters.
        semi_major = 6378.137
        semi_minor = 6356.752314245

        c_phi, s_phi = np.cos(phi), np.sin(phi)

        num = (semi_major ** 2 * c_phi) ** 2 + (semi_minor ** 2 * s_phi) ** 2
        denom = (semi_major * c_phi) ** 2 + (semi_minor * s_phi) ** 2
        radius = np.sqrt(num / denom)

        return radius

    def add_encounter(self, name: str, encType: Literal["I", "T"], ang: float,
                      angType: Literal["A", "N"], solar: Literal[-1, 0, 1]=0,
                      sca: float=0) -> None:
        """Define encounter for the ground location.

        This method allows for the specification of a ground/satellite
        encounter for the specific Earth-based location. This association will
        be used in window generation and encounter planning.

        Parameters
        ----------
        name : str
            Encounter identifier.
        encType : {"I", "T"}
            Specifies encounter category as either imaging or transmission.
        ang : float
            Angular constraint for the encounter in degrees.
        angType : {"A", "N"}
            String specifying the constraint angle as either the altitude or
            the off-nadir angle type.
        solar : {-1, 0, 1}, optional
            Defines sunlight allowance where -1 allows for windows at night, 0
            allows for windows at day or night, and 1 allows for windows at day.
        sca : float, optional
            Float specifying the minimum angle between a satellite's position
            vector and the Sun's position vector in a ground-location-centric
            reference system. Refer to notes for more information.

        Examples
        --------
        >>> ground_pos = GroundPosition(name="Toronto", coor=(43.6532, -79.3832))
        >>> ground_pos.add_encounter(name="CYYZ IMG, encType="I", ang=30,
        ...                          angType="N", solar=0, SCA=0)
        """

        encounter_object = EncounterSpec(name, encType, ang, angType, solar, sca)
        self.encounters[name] = encounter_object
