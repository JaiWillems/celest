

from celest.encounter.schedule.adaptation_utils import _conflict_degree
import random


def random_removal(request_list, q) -> list:
    """Randomly remove q requestrs from the request list.

    Parameters
    ----------
    request_list : list
        List of requests.
    q : int
        The number of requests to remove.

    Return
    ------
    list
        The request list with items removed.
    """

    while q > 0:

        i = random.randint(0, len(request_list) - 1)

        if request_list[i].is_scheduled:

            request_list[i].is_scheduled = False
            request_list[i].scheduled_idx = None
            request_list[i].scheduled_start = None
            request_list[i].scheduled_duration = None

            q -= 1

    return request_list


def priority_removal(request_list, q) -> list:
    """Remove q requests with the lowest priority from the request list.

    Parameters
    ----------
    request_list : list
        List of requests.
    q : int
        The number of requests to remove.

    Return
    ------
    list
        The request list with items removed.
    """

    request_list = sorted(request_list, key=lambda x: x.priority, reverse=True)

    i = 0
    while (q > 0) and (i < len(request_list)):

        if request_list[i].is_scheduled:

            request_list[i].is_scheduled = False
            request_list[i].scheduled_idx = None
            request_list[i].scheduled_start = None
            request_list[i].scheduled_duration = None

            q -= 1

        i += 1

    return request_list


def opportunity_removal(request_list, q) -> list:
    """Remove q requests with the most number of vtws from the request list.

    Parameters
    ----------
    request_list : list
        List of requests.
    q : int
        The number of requests to remove.

    Return
    ------
    list
        The request list with items removed.
    """

    request_list = sorted(request_list, key=lambda x: len(x.vtws), reverse=True)

    i = 0
    while (q > 0) and (i < len(request_list)):

        if request_list[i].is_scheduled:

            request_list[i].is_scheduled = False
            request_list[i].scheduled_idx = None
            request_list[i].scheduled_start = None
            request_list[i].scheduled_duration = None

            q -= 1

        i += 1

    return request_list


def conflict_removal(request_list, q) -> list:
    """Remove q requests with the highest conflict degree.

    Parameters
    ----------
    request_list : list
        List of requests.
    q : int
        The number of requests to remove.

    Returns
    -------
    list
        The request list with items removed.
    """

    cd_arr = []
    for i, request in enumerate(request_list):

        if request.is_scheduled:

            idx = request.scheduled_idx
            cd = _conflict_degree(request_list, request.vtws[idx])
            cd_arr.append((i, cd))

    for i, _ in sorted(cd_arr, key=lambda x: x[1], reverse=True):

        request_list[i].is_scheduled = False
        request_list[i].scheduled_idx = None
        request_list[i].scheduled_start = None
        request_list[i].scheduled_duration = None

        q -= 1

        if q == 0:
            break

    return request_list


_REMOVAL_FUNCTIONS = [
    random_removal,
    priority_removal,
    opportunity_removal,
    conflict_removal
]
