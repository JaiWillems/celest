

from celest.schedule.insertion_operators import insert_first_n_requests
from celest.schedule.request_handler import RequestIndices


def initialize_solution(request_handler):
    """Initial solution via a greedy heuristic.

    This function takes a `RequestHandler` and schedules, in place, a maximal
    number of requests by decreasing priority to get a good solution cost.

    Parameters
    ----------
    request_handler : RequestHandler

    Returns
    -------
    RequestHandler
        `request_handler` with initialized with the initial solution.
    """

    request_handler.sort_by_decreasing_priority()
    request_handler.sort_vtws_by_increasing_rise_time()
    insert_first_n_requests(request_handler, request_handler.number_of_requests)


def is_complete(request_handler):
    """Check if a solution is complete.

    A solution is complete if all requests are scheduled. A complete solution
    is the most optimal solution for a given set of requests.

    Parameters
    ----------
    request_handler : RequestHandler

    Returns
    -------
    boolean
        A boolean value showing if the solution is complete or not.
    """

    for request in request_handler:
        if not request[RequestIndices.is_scheduled]:
            return False

    return True


def cost(request_handler):
    """Return the cost of scheduled requests.

    The cost of scheduled requests is the negative sum of the scheduled request
    priorities.

    Parameters
    ----------
    request_handler : RequestHandler

    Returns
    -------
    int
        The cost of the scheduled requests in `request_handler`.
    """

    solution_cost = 0
    for request in request_handler:
        if request[RequestIndices.is_scheduled]:
            solution_cost -= request[RequestIndices.priority]

    return solution_cost
