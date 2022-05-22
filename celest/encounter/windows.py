

from celest.encounter._window_handling import VTW, VTWHandler
from celest.satellite.coordinate import Coordinate
from jplephem.spk import SPK
import numpy as np
import pkg_resources


def generate_vtw(satellite, location, vis_threshold, lighting=0, tol=1e-5) -> VTWHandler:
    """Return visible window times for a satellite-ground combination.

    This function generates the visible window times where the satellite has
    an elevation angle greater than `vis_threshold` and where lighting
    conditions are met.

    Parameters
    ----------
    satellite : Satellite
        Satellite object.
    location : Location
        Location of ground station.
    vis_threshold : float
        Visibility threshold in degrees.
    lighting : float, optional
        Lighting condition specifier.

        `-1` indicates nighttime encounters, `0` indicates anytime encounters,
        and `1` indicates daytime encounters.
    tol : float, optional
        Window start and end time tolerance in days.

    Returns
    -------
    VTWHandler
        The visible time windows.
    """

    if vis_threshold < 0 or vis_threshold >= 90:
        raise ValueError("Visibility threshold shall be in the range [0, 90).")
    if lighting not in (-1, 0, 1):
        raise ValueError("Valid lighting conditions include -1, 0, and 1.")

    window_locations = _get_base_vtws(satellite, location, vis_threshold)
    window_locations *= _lighting_constraint_factor(satellite, location, lighting)

    julian = satellite._julian
    roll, pitch, yaw = satellite.attitude(location, stroke=True)

    window_indices = _get_window_indices(window_locations, julian)
    if window_indices.size == 0:
        return None

    windows = VTWHandler()
    for window_indices in window_indices:

        rt1, rt2 = julian[window_indices[0] - 1], julian[window_indices[0]]
        rise_time = _root_find(window_locations, rt1, rt2, tol)

        st1, st2 = julian[window_indices[-1]], julian[window_indices[-1] + 1]
        set_time = _root_find(window_locations, st1, st2, tol)

        window = VTW(rise_time, set_time, roll, pitch, yaw)
        windows._add_window(window)

    return windows


def _get_base_vtws(satellite, location, vis_threshold):

    elevation = satellite._elevation(location)
    return _constraint_highlighter(elevation, "ge", vis_threshold)


def _constraint_highlighter(decision_var, constraint, ang):
    """Return raw window highlighter.

    This function applies the constraint between the decision variable and
    constraint angle to generate a Stroke that is one where the constraint is
    satisfied and zero elsewhere.

    Parameters
    ----------
    decision_var : Stroke
        Stroke object containing the decision variable.
    constraint : {"gt", "ge", "lt", "le", "eq"}
        Constraint type.

        The constraint is applied with the decision variable on the left and
        the constraint angle on the right. For example, `constraint="gt"`
        implies `decision_var > ang`.
    ang : float
        Constraint angle in units equivalent to `decision_var`.

    Returns
    -------
    Stroke
        Stroke object containing the constraint highlighter.
    """

    if constraint == "gt":
        return decision_var > ang
    elif constraint == "ge":
        return decision_var >= ang
    elif constraint == "lt":
        return decision_var < ang
    elif constraint == "le":
        return decision_var <= ang
    elif constraint == "eq":
        return decision_var == ang
    else:
        raise ValueError("Invalid constraint type.")


def _lighting_constraint_factor(satellite, location, lighting):

    if lighting == -1:
        return _get_night_constraint_factor(satellite._julian, location)
    elif lighting == 0:
        return 1
    elif lighting == 1:
        return _get_day_constraint_factor(satellite._julian, location)


def _get_day_constraint_factor(julian, location):

    return _get_lighting_constraint_factor(julian, location, 1)


def _get_night_constraint_factor(julian, location):

    return _get_lighting_constraint_factor(julian, location, -1)


def _get_lighting_constraint_factor(julian, location, lighting):

    sun_elevation, _ = _sun_coor(julian).horizontal(location, stroke=True)
    comparison = "gt" if lighting == 1 else "le"
    return _constraint_highlighter(sun_elevation, comparison, 0)


def _sun_coor(julian) -> Coordinate:
    """Return the sun's coordinates.

    Parameters
    ----------
    julian : np.ndarray
        1-D array containing Julian times in the J2000 epoch.

    Returns
    -------
    Coordinate
        Coordinate object containing the sun's position.
    """

    ephem = pkg_resources.resource_filename(__name__, '../data/de421.bsp')
    kernal = SPK.open(ephem)

    ssb2sunb_pos, ssb2sunb_vel = kernal[0, 10].compute_and_differentiate(julian)
    ssb2eb_pos, ssb2eb_vel = kernal[0, 3].compute_and_differentiate(julian)
    eb2e_pos, eb2e_vel = kernal[3, 399].compute_and_differentiate(julian)

    e2sun_pos = (ssb2sunb_pos - ssb2eb_pos - eb2e_pos).T
    e2sun_vel = (ssb2sunb_vel - ssb2eb_vel - eb2e_vel).T

    SPK.close(kernal)

    return Coordinate(e2sun_pos, e2sun_vel, "gcrs", julian)


def _get_window_indices(window_locations, julian):

    ind = np.arange(0, julian.size, 1) * window_locations(julian).astype(int)
    ind = ind[ind != 0]
    ind = np.split(ind, np.where(np.diff(ind) != 1)[0] + 1)

    return np.array(ind, dtype=object)


def _root_find(func, lower_bound, upper_bound, tolerance) -> float:
    """Return root of f(t) between lg and rg.

    Parameters
    ----------
    func : Stroke
        Stroke object representing the function to be evaluated.
    tl : float
        Lower bound of the search interval.
    upper_bound : float
        Upper bound of the search interval.
    tolerance : float
        Search tolerance.

    Returns
    -------
    float
        Root of `f(t)` between `lg` and `rg`.
    """

    left_value, right_value = lower_bound, upper_bound
    left_function_value = int(func(left_value))

    while abs(right_value - left_value) > tolerance:

        center_value = (left_value + right_value) / 2
        center_function_value = int(func(center_value))

        if abs(left_function_value - center_function_value):
            right_value = center_value
        else:
            left_value = center_value
            left_function_value = center_function_value

    return (left_value if left_function_value else right_value)
