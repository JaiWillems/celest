

from celest.coordinates.coordinate import Coordinate
from celest.coordinates.frames.azel import AzEl
from celest.coordinates.frames.gcrs import GCRS
from celest.coordinates.ground_location import GroundLocation
from celest.encounter.window_handling import VisibleTimeWindow, WindowCollection
from celest.satellite import Satellite
from celest.units.quantity import Quantity
from celest import units as u
from enum import Enum
from jplephem.spk import SPK
import numpy as np
import pkg_resources


class Lighting(Enum):
    """Enum for different lighting conditions.

    Attributes
    ----------
    NIGHTTIME : int
        Night time lighting condition.
    ANYTTIME : int
        No lighting condition.
    DAYTIME : int
        Day time lighting condition.
    """

    NIGHTTIME = -1
    ANYTIME = 0
    DAYTIME = 1


def generate_vtws(satellite: Satellite, location: GroundLocation,
                  vis_threshold: float, lighting:
                  Lighting=Lighting.ANYTIME) -> WindowCollection:
    """Return visible time windows for a satellite-ground encounter.

    This function determines the time windows where the satellite has
    an elevation angle greater than `vis_threshold` and where lighting
    conditions are met.

    Parameters
    ----------
    satellite : Satellite
    location : GroundLocation
        Location of ground station.
    vis_threshold : float
        Visibility threshold in degrees.
        
        The visibility threshold is the minimum elevation angle of the satellite
        as seen from `location` where the satellite will be in visual range of
        `location`.
    lighting : Lighting, optional
        Lighting condition specifier, defaults to anytime lighting conditions.

        Lighting conditions can be specified using the `Lighting` enum. The
        different options are `Lighting.NIGHTTIME`, `Lighting.ANYTIME`, and
        `Lighting.DAYTIME`.

    Returns
    -------
    WindowCollection
        The visible time windows.
    """

    if satellite.velocity is None:
        raise Warning("Attitude data cannot be initialized because satellite "
                      "velocity data is unavailable.")
    if vis_threshold < 0 or 90 < vis_threshold:
        raise ValueError("vis_threshold must be between 0 and 90 degrees.")

    satellite_azel = Coordinate(satellite.position).convert_to(AzEl, location)
    julian = satellite.position.time.to(u.jd2000)

    window_indices = np.where(satellite_azel.elevation.to(u.deg) > vis_threshold)[0]

    lighting_indices = _lighting_constraint_indices(julian, location, lighting)
    window_indices = np.intersect1d(window_indices, lighting_indices)

    window_indices = np.split(window_indices,
                              np.where(np.diff(window_indices) != 1)[0] + 1)
    window_indices = np.array(window_indices, dtype=object)

    attitude = satellite.attitude(location)

    windows = WindowCollection()
    for window_index_set in window_indices:
        rise_time = julian[window_index_set[0] - 1]
        set_time = julian[window_index_set[-1]]
        windows.add_window(VisibleTimeWindow(rise_time, set_time, attitude))

    return windows


def _lighting_constraint_indices(julian: np.ndarray, location: GroundLocation,
                                 lighting: Lighting) -> np.ndarray:
    if lighting == Lighting.NIGHTTIME:
        return _get_night_constraint_indices(julian, location)
    elif lighting == Lighting.ANYTIME:
        return np.indices(julian.shape)
    elif lighting == Lighting.DAYTIME:
        return _get_day_constraint_indices(julian, location)
    else:
        raise ValueError("Invalid lighting constraint.")


def _get_night_constraint_indices(julian: np.ndarray, location:
                                  GroundLocation) -> np.ndarray:
    sun_elevation = _get_sun_elevation(julian, location)
    return np.where(sun_elevation.to(u.deg) < 0)[0]


def _get_sun_elevation(julian: np.ndarray, location: GroundLocation) -> Quantity:
    sun_coordinates = Coordinate(_get_sun_gcrs(julian))
    return sun_coordinates.convert_to(AzEl, location).elevation


def _get_sun_gcrs(julian: np.ndarray) -> GCRS:
    ephem = pkg_resources.resource_filename(__name__, '../data/de421.bsp')
    kernal = SPK.open(ephem)

    ssb2sunb_pos = kernal[0, 10].compute(julian)
    ssb2eb_pos = kernal[0, 3].compute(julian)
    eb2e_pos = kernal[3, 399].compute(julian)
    e2sun_pos = (ssb2sunb_pos - ssb2eb_pos - eb2e_pos).T

    SPK.close(kernal)

    return GCRS(julian, e2sun_pos[:, 0], e2sun_pos[:, 1], e2sun_pos[:, 2], u.km)


def _get_day_constraint_indices(julian: np.ndarray, location:
                                GroundLocation) -> np.ndarray:
    sun_elevation = _get_sun_elevation(julian, location)
    return np.where(sun_elevation.to(u.deg) > 0)[0]
