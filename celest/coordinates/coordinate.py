

from celest.coordinates.ground_location import GroundLocation
from celest.coordinates.frames.azel import AzEl
from celest.coordinates.frames.gcrs import GCRS
from celest.coordinates.frames.itrs import ITRS
from celest.coordinates.frames.wgs84 import WGS84
from celest.coordinates.transforms import (
    _gcrs_to_itrs,
    _itrs_to_gcrs,
    _itrs_to_wgs84,
    _wgs84_to_itrs,
    _itrs_to_azel
)
from typing import Literal


class Coordinate:
    """Coordinate(position)

    Coordinate frame manipulations and transformations.

    Parameters
    ----------
    position : {ITRS, GCRS, WGS84}
        Input position.

    Methods
    -------
    convert_to(frame, location=None)
        Convert the current frame to `frame`.

    Notes
    -----
    Conversions between the AzEl, GCRS, ITRS, and WGS84 frames are supported
    except for conversions from the AzEl frame.

    Examples
    --------
    Let `gcrs` be an initialized `GCRS` frame. We can then convert to the `ITRS`
    frame:

    >>> c = Coordinate(gcrs)
    >>> itrs = c._convert_to(ITRS)

    By specifying a ground location, we can convert to the `AzEl` frame:

    >>> location = GroundLocation(43.6532, -79.3832, 76, u.deg, u.m)
    >>> azel = c._convert_to(AzEl, location)
    """

    def __init__(self, position: Literal[ITRS, GCRS, WGS84]):

        acceptable_frames = [ITRS, GCRS, WGS84]
        if position.__class__ not in acceptable_frames:
            raise ValueError(f"Input frame must be among {', '.join([i.__name__ for i in acceptable_frames])}")

        self._base_position = position

    def convert_to(self, frame: Literal[ITRS, GCRS, WGS84, AzEl],
                   location: GroundLocation = None):
        if isinstance(self._base_position, frame):
            return self._base_position
        elif isinstance(self._base_position, GCRS):
            return self._convert_from_gcrs(frame, location)
        elif isinstance(self._base_position, ITRS):
            return self._convert_from_itrs(frame, location)
        elif isinstance(self._base_position, WGS84):
            return self._convert_from_wgs84(frame, location)
        else:
            return NotImplemented

    def _convert_from_gcrs(self, frame: Literal[ITRS, WGS84, AzEl],
                           location: GroundLocation = None):
        itrs_position = _gcrs_to_itrs(self._base_position)
        if frame == ITRS:
            return itrs_position
        elif frame == WGS84:
            return _itrs_to_wgs84(itrs_position)
        elif frame == AzEl:
            return _itrs_to_azel(itrs_position, location)
        else:
            return NotImplemented

    def _convert_from_itrs(self, frame: Literal[GCRS, WGS84, AzEl],
                           location: GroundLocation = None):
        if frame == GCRS:
            return _itrs_to_gcrs(self._base_position)
        elif frame == WGS84:
            return _itrs_to_wgs84(self._base_position)
        elif frame == AzEl:
            return _itrs_to_azel(self._base_position, location)
        else:
            return NotImplemented

    def _convert_from_wgs84(self, frame: Literal[ITRS, GCRS, AzEl],
                            location: GroundLocation = None):
        itrs_position = _wgs84_to_itrs(self._base_position)
        if frame == ITRS:
            return itrs_position
        elif frame == GCRS:
            return _itrs_to_gcrs(itrs_position)
        elif frame == AzEl:
            return _itrs_to_azel(itrs_position, location)
        else:
            return NotImplemented

