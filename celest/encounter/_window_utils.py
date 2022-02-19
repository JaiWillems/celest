

from celest.encounter._encounter_math_utils import _analytical_encounter_ind
from celest.satellite.coordinate import Coordinate
from typing import Any, Literal
from jplephem.spk import SPK
import numpy as np
import pkg_resources


def _sun_itrs(julian: np.ndarray) -> np.ndarray:
    """Return Sun itrs position.

    Parameters
    ----------
    julian : np.ndarray
        1-D array containing Julian times in the J2000 Epoch.

    Returns
    -------
    np.ndarray
        2-D array containing columns of x, y, z itrs position data of the sun.
    """

    ephem = pkg_resources.resource_filename(__name__, '../data/de421.bsp')
    kernal = SPK.open(ephem)

    ssb2sunb = kernal[0, 10].compute(julian)
    ssb2eb = kernal[0, 3].compute(julian)
    eb2e = kernal[3, 399].compute(julian)
    e2sun = (ssb2sunb - ssb2eb - eb2e).T

    SPK.close(kernal)

    sun_itrs = Coordinate(e2sun, "gcrs", julian).itrs()

    return sun_itrs


def _get_ang(u: np.ndarray, v: np.ndarray) -> np.ndarray:
    """Calculate degree angle bewteen two vectors.

    Parameters
    ----------
    u, v : np.ndarray
        2-D arrays containing row-vectors.

    Returns
    -------
    np.ndarray
        1-D array containing the degree angle between rows of `u` and `v`.
    """

    num = np.sum(u * v, axis=1)
    denom = np.linalg.norm(u, axis=1) * np.linalg.norm(v, axis=1)
    ang = np.degrees(np.arccos(num / denom))

    return ang


def _window_encounter_ind(satellite: Any, location: Any, ang: float, form:
                          Literal[0, 1], sca: float=0, lighting:
                          Literal[-1, 0, 1]=0) -> np.ndarray:
    """Return valid encounter indices.

    Parameters
    ----------
    satellite : Satellite
        Satellite involved with a ground location encounter.
    location : GroundPosition
        Ground location involved with a satellite encounter.
    ang : float
        Encounter constraint angle.
    form : {0, 1}
        Specifies the constraing angle as the altitude angle if 0 or the
        off-nadir angle if 1.
    sca : float, optional
        Solar constraint angle, by default 0.
    lighting : {-1, 0, 1}, optional
        Encounter lighting conditions where -1 is for night encounters, 1 is
        for day encounters, and 0 is for no constraint; default is 0.

    Returns
    -------
    np.ndarray
        1-D array containing viable satellite encounter positions.
    """

    sat_itrs = satellite.itrs()
    time_data = satellite._julian

    ground_GEO = np.array([[location.lat, location.lon, 0]])
    gnd_itrs = satellite._geo_to_itrs(ground_GEO)

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
