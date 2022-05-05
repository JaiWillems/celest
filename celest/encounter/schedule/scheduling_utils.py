

from celest.encounter.schedule.adaptation_utils import (
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


def cost(request_list) -> float:
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
