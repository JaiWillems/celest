from enum import Enum


class FlightStates(Enum):
    NORMAL_POINTING = 0
    PRE_ENCOUNTER = 1
    ENCOUNTER = 2
    POST_ENCOUNTER = 3
