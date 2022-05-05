

from typing import Tuple


def _look_ang_time(pitch_stroke, look_ang, rise_time, set_time, tol=1e-5) -> float:
    """Return the time where the look angle is acheived.

    This method implements a bisection rooting algorithm to find the time where
    the look angle is acheived.

    Parameters
    ----------
    pitch_stroke : Stroke
        Stroke containing the pitch angle in degrees during the VTW.
    look_ang : float
        Look angle in degrees.
    rise_time : float
        Rise time of the VTW in Julian days.
    set_time : float
        Set time of the VTW in Julian days.
    tol : float, optional
        Tolerance for the bisection algorithm.

    Returns
    -------
    float
        Time where the look angle is acheived.
    """

    l = rise_time
    fl = pitch_stroke(l) - look_ang

    r = (rise_time + set_time) / 2
    fr = pitch_stroke(r) - look_ang

    if fl * fr > 0:
        return None

    while r - l > tol:

        c = (l + r) / 2
        fc = pitch_stroke(c) - look_ang

        if fc * fr > 0:
            r, fr = c, fc
        else:
            l, fl = c, fc

    return (l + r) / 2


def _get_OW_start(vtw, look_ang) -> float:
    """Return observation window start time.

    Parameters
    ----------
    vtw : VTW
        Visible time window under consideration.
    look_ang : float
        Encounter look angle.

    Returns
    -------
    float
        Observation window start time.
    """

    if look_ang is not None:
        start = _look_ang_time(vtw.pitch, look_ang, vtw.rise_time, vtw.set_time)
    else:
        start = (vtw.rise_time + vtw.set_time) / 2

    return start


def _insert_conflict(request_list, start, duration) -> bool:
    """Check if the insertion of a request conflicts with scheduled tasks.

    Parameters
    ----------
    request_list : list
        List of requests.
    start : float
        Start time of the inserting request in Julian days.
    duration : float
        Duration of the inserting request in seconds.

    Returns
    -------
    bool
        Whether the insertion conflicts with scheduled tasks.
    """

    for request in request_list:

        if request.is_scheduled:

            sch_start = request.scheduled_start
            sch_duration = request.scheduled_duration

            if sch_start + sch_duration / 86400 > start:
                return True

            if start + duration < sch_start:
                return True

    return False


def _insert_window(request_list, i) -> bool:
    """Attempt window insertion for request at index i.

    Parameters
    ----------
    request_list : list
        List of requests.
    i : int
        Request index.

    Returns
    -------
    bool
        True if window insertion was successful, False otherwise.
    """

    request = request_list[i]
    vtws = request.vtws

    for j, vtw in enumerate(vtws):

        start = _get_OW_start(vtw, request.look_ang)
        duration = request.duration

        if start + duration / 86400 > vtw.set_time:
            continue
        if start + duration / 86400 > request.deadline:
            continue
        if _insert_conflict(request_list, start, duration):
            continue

        request.is_scheduled = True
        request.scheduled_idx = j
        request.scheduled_start = start
        request.scheduled_duration = duration

        return True

    return False


def _overlap(vtw1, vtw2) -> bool:
    """Check if two visible time windows overlap.

    Parameters
    ----------
    vtw1, vtw2 : VTW

    Returns
    -------
    bool
        Returns True if an overlap exists, otherwise False.
    """

    s1, s2 = vtw1.rise_time, vtw2.rise_time
    e1, e2 = vtw1.set_time, vtw2.set_time

    if (e1 < s2) | (e2 < s1):
        return False
    else:
        return True


def _over(request_list, curr_vtw):
    """Determine the set of overlapping visible time windows.

    Parameters
    ----------
    request_list : list
        List of requests.
    curr_vtw : VTW

    Returns
    -------
    list
        List of visible time windows overlapping with `curr_vtw`.
    """

    overlap_vtws = []

    for request in request_list:
        for vtw in request.vtws:

            if _overlap(curr_vtw, vtw):
                overlap_vtws.append(vtw)

    return overlap_vtws


def _time_span(vtw1, vtw2) -> float:
    """Return the overlap between two visible time windows.

    Parameters
    ----------
    vtw1, vtw2 : VTW

    Returns
    -------
    float
        The overlap between the two visible time windows in Julian days.
    """

    s1, s2 = vtw1.rise_time, vtw2.rise_time
    e1, e2 = vtw1.set_time, vtw2.set_time

    if (e1 < s2) | (e2 < s1):
        return 0
    elif (s1 < s2) & (e2 < e1):
        return e2 - s2
    elif (s2 < s1) & (e1 < e2):
        return e1 - s1
    elif s2 < e1:
        return e1 - s2
    elif s1 < e2:
        return e2 - s1


def _conflict_degree(request_list, vtw) -> float:
    """Return the conflict degree of `vtw`.

    Parameters
    ----------
    request_list : list
        List of requests.
    vtw : VTW

    Returns
    -------
    float
        The conflict degree of `vtw`.
    """

    over = _over(request_list, vtw)
    time_sum = 0

    for tw in over:
        time_sum += _time_span(vtw, tw)

    cd = time_sum / len(over)

    return cd


def _min_conflict_degree(request_list, i) -> list:
    """Return the vtw conflict degrees for the ith request.

    Parameters
    ----------
    request_list : list
        List of requests.
    i : int
        Index of the request.

    Returns
    -------
    list
        List containing a tuple for each visible time window which holds the
        conflict degree and the index for the requests `vtws` attribute.

        The list is ordered in increasing order of conflict degree.
    """

    cd_arr = []

    for j, vtw in enumerate(request_list[i].vtws):

        cd = _conflict_degree(request_list, vtw)
        cd_arr.append((cd, j))

    cd_arr = sorted(cd_arr, key=lambda x: x[0])

    return cd_arr
