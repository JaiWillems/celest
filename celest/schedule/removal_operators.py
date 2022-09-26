

import random


def random_removal(request_handler, number_to_remove):
    """Remove n scheduled tasks at random using a uniform distribution.

    Parameters
    ----------
    request_handler : RequestHandler
    number_to_remove : int

    Returns
    -------
    RequestHandler
        `request_handler` with `number_to_remove` requests unscheduled.
    """

    if number_to_remove >= request_handler.number_of_scheduled_requests:
        request_handler.unschedule_all_requests()
        return request_handler

    unchecked_indices = [i for i in range(0, len(request_handler) - 1)]
    while number_to_remove > 0:
        request_index = random.choice(unchecked_indices)
        if request_handler.is_request_scheduled(request_index):
            request_handler.unschedule_request(request_index)
            number_to_remove -= 1
        unchecked_indices.remove(request_index)

    return request_handler


def priority_removal(request_handler, number_to_remove):
    """Remove n scheduled tasks in order of decreasing priority.

    Parameters
    ----------
    request_handler : RequestHandler
    number_to_remove : int

    Returns
    -------
    RequestHandler
        `request_handler` with `number_to_remove` requests unscheduled.
    """

    request_handler.sort_by_decreasing_priority()
    return _remove_first_n_scheduled_requests(request_handler, number_to_remove)


def _remove_first_n_scheduled_requests(request_handler, number_to_remove):
    if number_to_remove >= request_handler.number_of_scheduled_requests:
        request_handler.unschedule_all_requests()
        return request_handler

    index = 0
    while number_to_remove > 0:
        if request_handler.is_request_scheduled(index):
            request_handler.unschedule_request(index)
            number_to_remove -= 1
        index += 1

    return request_handler


def opportunity_removal(request_handler, number_to_remove):
    """Remove n scheduled tasks in order of decreasing opportunity.

    Parameters
    ----------
    request_handler : RequestHandler
    number_to_remove : int

    Returns
    -------
    RequestHandler
        `request_handler` with `number_to_remove` requests unscheduled.

    Notes
    -----
    The opportunity of a request is measured as the number of opportunities for
    the request to be fulfilled. Opportunity removal removes requests with the
    highest number of opportunities since they can more easily be inserted into
    a solution without conflict.
    """

    request_handler.sort_by_decreasing_opportunity()
    return _remove_first_n_scheduled_requests(request_handler, number_to_remove)


def conflict_removal(request_handler, number_to_remove):
    """Remove n scheduled tasks in order of decreasing conflict degree.

    Parameters
    ----------
    request_handler : RequestHandler
    number_to_remove : int

    Returns
    -------
    RequestHandler
        `request_handler` with `number_to_remove` requests unscheduled.

    Notes
    -----
    The conflict degree of a request is a measure of how much the opportunity
    windows of a request overlap with those of other requests. Conflict removal
    removes requests with the highest conflict degree as these are most likely
    to limit the number of scheduled requests.
    """

    request_handler.sort_by_decreasing_conflict_degree()
    return _remove_first_n_scheduled_requests(request_handler, number_to_remove)


_REMOVAL_FUNCTIONS = [
    random_removal,
    priority_removal,
    opportunity_removal,
    conflict_removal
]
