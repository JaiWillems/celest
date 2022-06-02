

from celest.schedule.insertion_operators import insert_first_n_requests
from celest.schedule.request_handler import RequestIndices


def initialize_solution(request_handler):

    request_handler.sort_by_decreasing_priority()
    request_handler.sort_vtws_by_increasing_rise_time()
    insert_first_n_requests(request_handler, request_handler.number_of_requests)


def is_complete(request_handler):

    for request in request_handler:
        if not request[RequestIndices.is_scheduled]:
            return False

    return True


def cost(request_handler):

    solution_cost = 0
    for request in request_handler:
        if request[RequestIndices.is_scheduled]:
            solution_cost -= request[RequestIndices.priority]

    return solution_cost
