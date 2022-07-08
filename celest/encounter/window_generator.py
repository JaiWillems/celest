

from celest.coordinates.coordinate import Coordinate
from celest.coordinates.frames.azel import AzEl
from celest.coordinates.frames.gcrs import GCRS
from celest.coordinates.ground_location import GroundLocation
from celest.encounter.window_handling import VisibleTimeWindow, WindowHandler
from celest.satellite import Satellite
from celest import units as u
from enum import Enum
from jplephem.spk import SPK
import numpy as np
import pkg_resources


class Lighting(Enum):
    NIGHTTIME = -1
    ANYTIME = 0
    DAYTIME = 1


def generate_vtw(satellite: Satellite, location: GroundLocation,
                 vis_threshold: float, lighting:
                 Lighting=Lighting.ANYTIME) -> WindowHandler:
    """Return visible window times for a satellite-ground set.

    This function generates the visible time windows where the satellite has
    an elevation angle greater than `vis_threshold` and where lighting
    conditions are met.

    Parameters
    ----------
    satellite : Satellite
        Satellite object.
    location : GroundLocation
        Location of ground station.
    vis_threshold : float
        Visibility threshold in degrees.
    lighting : float, optional
        Lighting condition specifier, defaults to anytime lighting conditions.

        Lighting conditions can be specified using the `Lighting` enum. The
        different options are `Lighting.NIGHTTIME`, `Lighting.ANYTIME`, and
        `Lighting.DAYTIME`.

    Returns
    -------
    WindowHandler
        The visible time windows.
    """

    if satellite.velocity is None:
        raise Warning("Attitude data cannot be initialized because satellite "
                      "velocity data is unavailable.")
    if vis_threshold < 0 or 90 < vis_threshold:
        raise ValueError("vis_threshold must be between 0 and 90 degrees.")

    satellite_azel = Coordinate(satellite.position).convert_to(AzEl, location)
    satellite_elevation = satellite_azel.elevation.to(u.deg).data
    julian = satellite.position.time.data

    window_indices = np.where(satellite_elevation > vis_threshold)[0]

    lighting_indices = _lighting_constraint_indices(julian, location, lighting)
    window_indices = np.intersect1d(window_indices, lighting_indices)

    window_indices = np.split(window_indices,
                              np.where(np.diff(window_indices) != 1)[0] + 1)
    window_indices = np.array(window_indices, dtype=object)

    attitude = satellite.attitude(location)

    windows = WindowHandler()
    for window_index_set in window_indices:
        rise_time = julian[window_index_set[0] - 1]
        set_time = julian[window_index_set[-1]]
        windows.add_window(VisibleTimeWindow(rise_time, set_time, attitude))

    return windows


def _lighting_constraint_indices(julian, location, lighting):
    if lighting == Lighting.NIGHTTIME:
        return _get_night_constraint_indices(julian, location)
    elif lighting == Lighting.ANYTIME:
        return np.indices(julian.shape)
    elif lighting == Lighting.DAYTIME:
        return _get_day_constraint_indices(julian, location)


def _get_night_constraint_indices(julian, location):
    sun_elevation = _get_sun_elevation(julian, location)
    return np.where(sun_elevation.to(u.deg).data < 0)[0]


def _get_sun_elevation(julian, location):
    sun_coordinates = Coordinate(_sun_coordinates(julian))
    return sun_coordinates.convert_to(AzEl, location).elevation


def _sun_coordinates(julian) -> GCRS:
    ephem = pkg_resources.resource_filename(__name__, '../data/de421.bsp')
    kernal = SPK.open(ephem)

    ssb2sunb_pos = kernal[0, 10].compute(julian)
    ssb2eb_pos = kernal[0, 3].compute(julian)
    eb2e_pos = kernal[3, 399].compute(julian)
    e2sun_pos = (ssb2sunb_pos - ssb2eb_pos - eb2e_pos).T

    SPK.close(kernal)

    return GCRS(julian, e2sun_pos[:, 0], e2sun_pos[:, 1], e2sun_pos[:, 2], u.km)


def _get_day_constraint_indices(julian, location):
    sun_elevation = _get_sun_elevation(julian, location)
    return np.where(sun_elevation.to(u.deg).data > 0)[0]
