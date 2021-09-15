"""Localize ground position information.

This module contains the `GoundPosition` class which stores ground location
data and is an instrumental part of encounter planning by managing the ground
location and satellite encounter relationship.
"""


from celest.core.decorators import set_module
from typing import Literal, Tuple
import numpy as np


_image_encounter = {
    "type": "I",
    "angType": "N",
    "lighting": 1,
    "sca": 0
}


_data_link_encounter = {
    "type": "T",
    "angType": "A",
    "lighting": 0,
    "sca": 30
}


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
    encType : {"image", "data_link"}
        Name of the encounter type as either an imaging encounter or a data
        link encounter. See notes for encounter specifications.
    ang : float
        Limit encounter angle in degrees.

    Attributes
    ----------
    name : str
        Name of the ground location used for identification and indexing.
    coor : tuple
        Geodetic coordinates for the location given in degrees and decimals.
    radius : float
        Earth's radius at the given coordinates.
    ang : float
        Angular constraint for the encounter in degrees.
    type : {"I", "T"}
        Specifies encounter category as either imaging or transmission.
    ang_type : {"A", "N"}
        String specifying the constraint angle as either the altitude or
        off-nadir angle type (see notes).
    lighting : {-1, 0, 1}
        Defines sunlight constraint where -1 gets windows at night, 0 gets
        windows at day or night, and 1 gets windows at day.
    solar_constraint_angle : float
        Float specifying the minimum angle between a satellite's position
        vector and the Sun's position vector in a ground-location-centric
        reference system. Refer to notes for more information.

    Methods
    -------
    _radius(obsCoor)
        Calculate the Earth's geocentric radius using WGS84.
    
    Notes
    -----
    The altitude angle is the position of an object measured in increasing
    degrees from the horizon. The off-nadir angle is the acute angle between
    the satellites nadir and the line joining the satellite with the ground
    location.

    When the "image" encounter type is selected, the valid encounter region
    is where the satellite's off-nadir angle is less then the constraint angle.
    The imaging encounters are constrained by daylight and no solar constraint
    angle.

    When the "data_link" encounter type is selected, the valid encounter region
    is where the satellites altitude angle is greater than the constraint
    angle. The data link encunters will be calculated at any time of day and
    will have a 30 degree constraint angle.

    In some instances of ground to satellite communication, hardware damage
    can be incurred when the ground station is within a certain angle of the
    Sun. The solar constraint angle allows encounters to be calculated out of
    direct alignment of the Sun by invalidating encounter regions where the
    Sun is behind or close to the satellite as seen from the ground station
    assuming the ground station is actively tracking the satellite.
    """

    def __init__(self, name: str, coor: Tuple[float, float], encType:
                 Literal["image", "data_link"], ang: float) -> None:
        """Initialize attributes."""

        self.name = name
        self.coor = coor
        self.radius = self._radius(coor[0])

        self.ang = ang
        
        if encType == "image":

            self.type = _image_encounter["type"]
            self.ang_type = _image_encounter["angType"]
            self.lighting = _image_encounter["lighting"]
            self.solar_constraint_angle = _image_encounter["sca"]

        elif encType == "data_link":

            self.type = _data_link_encounter["type"]
            self.ang_type = _data_link_encounter["angType"]
            self.lighting = _data_link_encounter["lighting"]
            self.solar_constraint_angle = _data_link_encounter["sca"]

        else:

            raise Exception("The encType parameter must be either \"image\" or \"data_link\"")

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
