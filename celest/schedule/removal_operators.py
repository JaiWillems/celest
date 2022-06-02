

import random


def random_removal(request_handler, number_to_remove):

    if number_to_remove > request_handler.number_of_scheduled_requests:
        request_handler.unschedule_all_requests()
        return request_handler

    while (number_to_remove > 0):
        request_index = random.randint(0, len(request_handler.requests) - 1)
        if request_handler.is_request_scheduled(request_index):
            request_handler.unschedule_request(request_index)
            number_to_remove -= 1

    return request_handler


def priority_removal(request_handler, number_to_remove):

    if number_to_remove >= request_handler.number_of_scheduled_requests:
        request_handler.unschedule_all_requests()
        return request_handler

    request_handler.sort_by_decreasing_priority()
    return remove_first_n_scheduled_requests(request_handler, number_to_remove)


def remove_first_n_scheduled_requests(request_handler, number_to_remove):

    index = 0
    while (number_to_remove > 0):
        if request_handler.is_request_scheduled(index):
            request_handler.unschedule_request(index)
            number_to_remove -= 1
        index += 1

    return request_handler


def opportunity_removal(request_handler, number_to_remove):

    if number_to_remove >= request_handler.number_of_scheduled_requests:
        request_handler.unschedule_all_requests()
        return request_handler

    request_handler.sort_by_decreasing_opportunity()
    return remove_first_n_scheduled_requests(request_handler, number_to_remove)


def conflict_removal(request_handler, number_to_remove):

    if number_to_remove >= request_handler.number_of_scheduled_requests:
        request_handler.unschedule_all_requests()
        return request_handler

    request_handler.sort_by_decreasing_conflict_degree()
    return remove_first_n_scheduled_requests(request_handler, number_to_remove)


_REMOVAL_FUNCTIONS = [
    random_removal,
    priority_removal,
    opportunity_removal,
    conflict_removal
]
