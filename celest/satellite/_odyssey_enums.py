"""Enumerations defined for `Satellite.generate_pointing_profiles`."""


from enum import Enum


class FlightStates(Enum):
    """Define satellite states as enumerations."""

    NORMAL_POINTING = 0
    PRE_ENCOUNTER = 1
    ENCOUNTER = 2
    POST_ENCOUNTER = 3
