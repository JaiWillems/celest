

from celest.encounter._window_handling import Window, Windows
from celest.satellite.coordinate import Coordinate
from jplephem.spk import SPK
import numpy as np
import pkg_resources


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


def _sun_coor(julian):
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

    ssb2sunb = kernal[0, 10].compute(julian)
    ssb2eb = kernal[0, 3].compute(julian)
    eb2e = kernal[3, 399].compute(julian)
    e2sun = (ssb2sunb - ssb2eb - eb2e).T

    SPK.close(kernal)

    sun_coor = Coordinate(e2sun, "gcrs", julian)

    return sun_coor


def _root_find(f, tl, tr, tol):
    """Return root of f(t) between lg and rg.

    Parameters
    ----------
    f : Stroke
        Stroke object representing the function to be evaluated.
    tl : float
        Lower bound of the search interval.
    tr : float
        Upper bound of the search interval.
    tol : float
        Search tolerance.

    Returns
    -------
    float
        Root of `f(t)` between `lg` and `rg`.
    """

    l, r = tl, tr
    fl = int(f(l))

    while abs(r - l) > tol:

        c = (l + r) / 2
        fc = int(f(c))

        if abs(fl - fc) == 1:
            r = c
        else:
            l = c
            fl = fc

    return (l if fl == 1 else r)


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

    ua = None if u.ndim == 1 else 1
    va = None if v.ndim == 1 else 1

    n = np.sum(u * v, axis=(ua or va))
    d = np.linalg.norm(u, axis=ua) * np.linalg.norm(v, axis=va)
    ang = np.degrees(np.arccos(n / d))

    return ang


def generate(satellite, location, enc, ang, lighting=0, tol=1e-5):
    """Return windows for a given encounter.

    Parameters
    ----------
    satellite : Satellite
        Satellite object.
    location : Location
        Location of ground station.
    enc : {"image", "data_link"}
        Type of encounter.
    ang : float
        Constraint angle in degrees.
    lighting : float, optional
        Lighting condition specifier.

        `-1` indicates nighttime encounters, `0` indicates anytime encounters,
        and `1` indicates daytime encounters.
    tol : float, optional
        Window start and end time tolerance in days.

    Returns
    -------
    Windows
        The satellite window details.
    """

    # Determine windows based on the type of encounter.
    if enc == "image":
        dv = satellite._altitude(location)
        raw_windows_1 = _constraint_highlighter(dv, "ge", 0)

        dv = satellite.off_nadir(location, stroke=True)
        raw_windows_2 = _constraint_highlighter(dv, "lt", ang)

        raw_windows = raw_windows_1 * raw_windows_2
    elif enc == "data_link":
        dv = satellite._altitude(location)
        raw_windows = _constraint_highlighter(dv, "gt", ang)
    else:
        raise ValueError("Invalid encounter type.")

    if lighting != 0:
        sun_coor = _sun_coor(satellite._julian)
        sun_alt, _ = sun_coor.horizontal(location, stroke=True)
        comp = "gt" if lighting == 1 else "le"
        raw_windows = raw_windows * _constraint_highlighter(sun_alt, comp, 0)

    julian = satellite._julian
    ind = np.arange(0, julian.size, 1) * raw_windows(julian).astype(int)
    ind = ind[ind != 0]
    ind = np.split(ind, np.where(np.diff(ind) != 1)[0] + 1)
    ind = np.array(ind, dtype=object)

    windows = Windows()

    # Populate Windows object.
    if ind.size == 0:
        return windows

    for i in ind:

        st1, st2 = julian[i[0] - 1], julian[i[0]]
        start = _root_find(raw_windows, st1, st2, tol)

        et1, et2 = julian[i[-1]], julian[i[-1] + 1]
        end = _root_find(raw_windows, et1, et2, tol)

        window = Window(satellite, location, start, end, enc, ang, lighting)
        windows._add_window(window)

    return windows
