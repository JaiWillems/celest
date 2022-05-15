

from celest.schedule.adaptation_utils import (
    _get_OW_start,
    _insert_conflict
)


def _initial_solution(request_list):
    """Generate a greedy initial solution.

    Parameters
    ----------
    request_list : list
        List of requests.

    Returns
    -------
    list
        List of requests with items scheduled.
    """

    request_list = sorted(request_list, key=lambda x: x.priority, reverse=True)

    for i, request in enumerate(request_list):

        for j, vtw in enumerate(request.vtws):

            start = _get_OW_start(vtw, request.look_ang)
            duration = request.duration

            if start is None:
                continue
            if start < vtw.rise_time:
                continue
            if start + duration / 86400 > vtw.set_time:
                continue
            if start + duration / 86400 > request_list[i].deadline:
                continue
            if _insert_conflict(request_list, start, duration):
                continue

            request_list[i].is_scheduled = True
            request_list[i].scheduled_idx = j
            request_list[i].scheduled_start = start
            request_list[i].scheduled_duration = duration

            break

    return request_list


def _is_complete(request_list) -> bool:
    """Check if the schedule is complete.

    Parameters
    ----------
    request_list : list
        List of requests.

    Returns
    -------
    bool
        True if the schedule is complete, False otherwise.
    """

    for request in request_list:

        if not request.is_scheduled:
            return False

    return True


def _cost(request_list) -> float:
    """Evaluate the cost of the schedule.

    Parameters
    ----------
    request_list : list
        List of requests.

    Returns
    -------
    float
        Cost of the schedule.
    """

    cost = 0

    for request in request_list:
        cost -= request.priority * request.is_scheduled

    return cost
