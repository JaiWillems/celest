"""Localize encounter information.

The _encounterspec module contains the `EncounterSpec` object to localize
encounter information specific to a `GroundPosition` Object.
"""


from celest.core.decorators import set_module
from typing import Literal
import numpy as np
import pandas as pd


@set_module('celest.encounter')
class EncounterSpec(object):
    """Localize encounter information.

    The `EncounterSpec` object localizes encounter information specific to a
    `GroundPosition` Object.

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
        off-nadir angle type. Refer to notes for greater detail in the angle
        definitions.
    SCA : float, optional
        Float specifying the minimum angle between a satellite's position
        vector and the Sun's position vector in a ground-location-centric
        reference system. Refer to notes for more information.
    solar : {-1, 0, 1}, optional
        Defines sunlight constraint where -1 gets windows at night, 0 gets
        windows at day or night, and 1 gets windows at day.

    Attributes
    ----------
    name : str
        Encounter identifier.
    type : {"I", "T"}
        Specifies encounter category as either imaging or transmission.
    ang : float
        Angluar constraint for the encounter in degrees.
    ang_type : {"A", "N"}
        String specifying the constraint angle as either the altitude, or the
        off-nadir angle type.
    solar : {-1, 0, 1}, optional
        Defines sunlight allowance where -1 allows for windows at night, 0
        allows for windows at day or night, and 1 allows for windows at day.
    solar_constraint_ang : float, optional
        Float specifying the minimum angle between a satellite's position
        vector and the Sun's position vector in a ground-location-centric
        reference system. Refer to notes for more information.
    windows : np.array
        Array of shape (n, 3) of window start, end, and elapsed seconds data.
    length : int
        Length, n, of data attributes.
    encounter_indices : np.array
        Array containing arrays of encounter indices.
    
    Notes
    -----
    The altitude angle is defined as the angle of a celestial object measured
    in increasing degrees from the horizon. A valid encounter region with the
    altitude angle are the regions where `alt_angle > ang`. The off-nadir
    angle is defined as the angle of an Earth based location measured in
    increasing degrees from a satellite's nadir. A valid encounter region with
    the off-nadir angle are the regions where `nadir_ang < ang`.

    In some instances of ground to satellite communication, hardware damage
    can be incurred when the ground station is within a certain angle of the
    sun. The solar constraint angle allows encounters to be calculated out of
    direct alignment of the sun by invalidating encounter regions where the
    sun is behind or close to the satellite as seen from the ground station
    assuming the ground station is actively tracking the satellite.
    """
    
    def __init__(self, name: str, encType: Literal["I", "T"], ang: float,
                 angType: Literal["A", "N"], solar: Literal[-1, 0, 1]=0,
                 SCA: float=0) -> None:
        """Initialize attributes."""

        self.name = name
        self.type = encType
        self.ang = ang
        self.ang_type = angType
        self.solar = solar
        self.solar_constraint_ang = SCA
        self.windows = None
        self.length = None
        self.encounter_indices = None
    
    def __str__(self) -> str:
        """Define `EncounterSpec` informaiton string."""

        data = np.array([self.name, self.type, self.ang, self.ang_type, self.solar])
        index = np.array(["Name:", "Encounter Type:", "Constraint Angle:",
                          "Angle Type:", "Solar Allowance:", "Solar Constraint Angle:"])

        df = pd.DataFrame(data=data, index=index)

        return df.to_string()
