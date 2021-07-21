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
        off-nadir angle type.
    maxAng : bool
        Defines the contraint angle as a maximum constraint if True or as
        minimum constraint if False. Note that the off-nadir angle is
        measured to increase away from nadir.
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
    max_ang : bool
        Defines the contraint angle as a maximum constraint if True or as
        minimum constraint if False. Note that the off-nadir angle is
        measured to increase away from nadir.
    solar : {-1, 0, 1}, optional
        Defines sunlight constraint where -1 gets windows at night, 0 gets
        windows at day or night, and 1 gets windows at day.
    windows : np.array
        Array of shape (n, 3) of window start, end, and elapsed seconds data.
    length : int
        Length, n, of data attributes.
    encounter_indices : np.array
        Array containing arrays of encounter indices.
    """
    
    def __init__(self, name: str, encType: Literal["I", "T"], ang: float,
                 angType: Literal["A", "N"], maxAng: bool, solar:
                 Literal[-1, 0, 1]=0) -> None:
        """Initialize attributes."""

        self.name = name
        self.type = encType
        self.ang = ang
        self.ang_type = angType
        self.max_ang = maxAng
        self.solar = solar
        self.windows = None
        self.length = None
        self.encounter_indices = None
    
    def __str__(self) -> str:
        """Define `EncounterSpec` informaiton string."""

        data = np.array([self.name, self.type, self.ang, self.ang_type,
                         self.max_ang, self.solar])
        index = np.array(["Name:", "Encounter Type:", "Constraint Angle:",
                          "Angle Type:", "Is Max:", "Solar Constraint:"])

        df = pd.DataFrame(data=data, index=index)

        return df.to_string()