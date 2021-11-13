"""Generate window sets for specified encounters.

This module contains functions to generate window sets for various satellite
and ground-location encounters.
"""


from celest.core.decorators import set_module
from celest.core.interpolation import _interpolate
from celest.encounter._window_handling import Window, Windows
from celest.encounter._window_utils import _window_encounter_ind
from typing import Any, Literal
import numpy as np


_image_encounter = {
    "type": "I",
    "angType": 1,  # Off-nadir angle.
    "lighting": 1,
    "sca": 0
}


_data_link_encounter = {
    "type": "T",
    "angType": 0,  # Altitude angle.
    "lighting": 0,
    "sca": 30
}


@set_module('celest.encounter.windows')
def generate(satellite: Any, location: Any, enc: Literal["image", "data link"],
             ang: float, factor: int=5) -> Windows:
    """Return encounter windows.

    Parameters
    ----------
    satellite : Satellite
        Satellite taking part in ground interactions.
    location : GroundPosition
        Ground location of imaging site or ground station.
    enc : {"image", "data link"}
        Type of encounter as being an imaging or data link encounter.
    ang : float
        Encounter contraint angle in degrees.
    factor : int, optional
        Data interpolation factor for more precise encounter information, by
        default 5.

    Returns
    -------
    Windows
        The encounter opportunities between the satellite and ground location
        of the type specified.
    
    Notes
    -----
    Imaging encounters are use the off-nadir angle, measured in increasing
    degrees from the satellite's nadir to the ground location. When the
    off-nadir angle is used, the input `ang` provides a maximum constraint.
    Data link encounters use the altitude angle, measured in increasing
    degrees from the ground location's horizon to the satellite. When the
    altitude angle is used, the input `ang` provides a minimum constraint.

    Examples
    --------
    >>> toronto = GroundPosition(latitude=43.65, longitude=-79.38)
    >>> toronto_dl = windows.generate(satellite, toronto, "data link", 30)
    """

    windows = Windows()

    if enc == "image":
        ang_type = _image_encounter["angType"]
        lighting = _image_encounter["lighting"]
        sca = _image_encounter["sca"]
    else:
        ang_type = _data_link_encounter["angType"]
        lighting = _data_link_encounter["lighting"]
        sca = _data_link_encounter["sca"]

    if factor > 1:

        if ang_type:
            off_nadir = satellite.position.off_nadir(location)
            ind = np.where(off_nadir < ang)[0]

        elif not ang_type:
            alt, _ = satellite.position.horizontal(location)
            ind = np.where(alt > ang)[0]

        ind = np.split(ind, np.where(np.diff(ind) != 1)[0] + 1)
        ind = np.array(ind, dtype=object)

        julian_interp = _interpolate(satellite.time.julian(), factor, 2, ind)
        eci_interp = _interpolate(satellite.position.gcrs(), factor, 2, ind)

        satellite.position._GCRS = eci_interp
        satellite.position._ITRS = None
        satellite.position._GEO = None
        satellite.position.length = eci_interp.shape[0]
        satellite.position.time._julian = julian_interp
        satellite.position.time._length = julian_interp.size
        satellite.time = satellite.position.time

    enc_ind = _window_encounter_ind(satellite, location, ang, ang_type, sca, lighting)

    window_ind = np.split(enc_ind, np.where(np.diff(enc_ind) != 1)[0] + 1)
    window_ind = np.array(window_ind, dtype=object)

    times = satellite.time.julian()

    if window_ind.size != 0:

        n = len(window_ind)

        for j in range(n):

            start = times[window_ind[j][0]]
            end = times[window_ind[j][-1]]
            window = Window(satellite, location, start, end, enc, ang, lighting, sca)

            windows._add_window(window)
            
    return windows
