"""Window generation dependencies.

This module contains implementation specific functions used in identifying
viable encounters.
"""

from celest.encounter._encounter_math_utils import _analytical_encounter_ind
from celest.satellite.coordinate import Coordinate
from celest.satellite.time import Time
from typing import Any, Literal
from jplephem.spk import SPK
import numpy as np
import pkg_resources


def _sun_itrs(julian: np.ndarray) -> np.ndarray:
    """Return Sun itrs position.

    Parameters
    ----------
    julian : np.ndarray
        Array of shape (n,) containing Julian time data.

    Returns
    -------
    np.ndarray
        Array of shape (n, 3) containing Sun ECEF positions.
    """

    ephem = pkg_resources.resource_filename(__name__, '../data/de421.bsp')
    kernal = SPK.open(ephem)

    ssb2sunb = kernal[0, 10].compute(julian)
    ssb2eb = kernal[0, 3].compute(julian)
    eb2e = kernal[3, 399].compute(julian)
    e2sun = (ssb2sunb - ssb2eb - eb2e).T

    SPK.close(kernal)

    sun_itrs = Coordinate(e2sun, "gcrs", Time(julian)).itrs()

    return sun_itrs

def _get_ang(u: np.ndarray, v: np.ndarray) -> np.ndarray:
    """Calculate degree angle bewteen two vectors.

    Parameters
    ----------
    u, v : np.ndarray
        Arrays of shape (n,3) with rows of XYZ cartesian data.

    Returns
    -------
    np.ndarray
        Array of shape (n,) containing the degree angle between the two
        arrays.
    """

    num = np.sum(u * v, axis=1)
    denom = np.linalg.norm(u, axis=1) * np.linalg.norm(v, axis=1)
    ang = np.degrees(np.arccos(num / denom))

    return ang

def _window_encounter_ind(satellite: Any, location: Any, ang: float, form:
                          Literal[0, 1], sca: float=0, lighting:
                          Literal[-1, 0, 1]=0) -> np.ndarray:
    """Return indices for possible encounters.

    Parameters
    ----------
    satellite : Satellite
        Satellite involved with a ground location encounter.
    location : GroundPosition
        Ground location involved with a satellite encounter.
    ang : float
        Encounter constraint angle.
    form : {0, 1}
        Angle type used where 0 is for the altitude angle and 1 is for the
        off-nadir angle.
    sca : float, optional
        Solar constraint angle, by default 0
    lighting : {-1, 0, 1}, optional
        lighting conditions for encounters where -1 is for night time
        encounters, 0 is for all-time encoutners, and 1 is for day time
        encounters, by default 0

    Returns
    -------
    np.ndarray
        Array containing indices defining viable satellite encounter positions.
    """

    sat_itrs = satellite.position.itrs()
    time_data = satellite.time.julian()

    ground_GEO = [location.lat, location.lon]
    ground_GEO = np.repeat(np.array([ground_GEO]), time_data.size, axis=0)
    gnd_itrs = Coordinate(ground_GEO, "geo", Time(time_data)).itrs()

    enc_params = [sat_itrs, gnd_itrs, ang, form]
    analytical_ind = _analytical_encounter_ind(*enc_params)

    sun_itrs = _sun_itrs(time_data)
    sca_angs = _get_ang(sat_itrs - gnd_itrs, sun_itrs - gnd_itrs)
    sca_ind = np.where(sca_angs > sca)[0]

    enc_ind = np.intersect1d(analytical_ind, sca_ind)

    if lighting != 0:
        sun_angs = _get_ang(gnd_itrs, sun_itrs)

        if lighting == 1:
            # Isolate indices occuring at day.
            day_ind = np.where(sun_angs < 90)
            enc_ind = np.intersect1d(enc_ind, day_ind)
        else:
            # Isolate indices occuring at night.
            night_ind = np.where(sun_angs >= 90)
            enc_ind = np.intersect1d(enc_ind, night_ind)

    return enc_ind
