

from celest.encounter.schedule.adaptation_utils import (
    _insert_window,
    _min_conflict_degree,
    _get_OW_start,
    _insert_conflict
)


def greedy_insertion(request_list, q):
    """Insert requests with highest priority.

    Parameters
    ----------
    request_list : list
        List of requests.
    q : int
        The number of requests to insert.

    Returns
    -------
    list
        List of requests with items inserted.
    """

    request_list = sorted(request_list, key=lambda x: x.priority, reverse=True)

    for i, request in enumerate(request_list):

        if not request.is_scheduled:
            if _insert_window(request_list, i):
                q -= 1
        if q == 0:
            break

    return request_list


def minimum_opportunity_insertion(request_list, q):
    """Insert requests with the fewest number of vtws.

    Parameters
    ----------
    request_list : list
        List of requests.
    q : int
        The number of requests to insert.

    Returns
    -------
    list
        List of requests with items inserted.
    """

    request_list = sorted(request_list, key=lambda x: len(x.vtws))

    for i, request in enumerate(request_list):

        if not request.is_scheduled:
            if _insert_window(request_list, i):
                q -= 1
        if q == 0:
            break

    return request_list


def minimum_conflict_insertion(request_list, q):

    for i, request in enumerate(request_list):

        if not request.is_scheduled:
            min_vtw, j = _min_conflict_degree(request_list, i)

            start = _get_OW_start(min_vtw, request.look_ang)
            duration = request.look_ang

            if start + duration / 86400 > min_vtw.set_time:
                continue
            if start + duration / 86400 > request.deadline:
                continue
            if _insert_conflict(request_list, start, duration):
                continue

            request.is_scheduled = True
            request.schedule_idx = j
            request.scheduled_start = start
            request.scheduled_duration = duration

            q -= 1

        if q == 0:
            break


_INSERTION_FUNCTIONS = [
    greedy_insertion,
    minimum_opportunity_insertion,
    minimum_conflict_insertion
]
