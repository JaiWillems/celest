

from celest.schedule.alns import ALNS
from celest.schedule.insertion_operators import _INSERTION_FUNCTIONS
from celest.schedule.removal_operators import _REMOVAL_FUNCTIONS
from celest.schedule.scheduling_utils import (
    _initial_solution,
    _is_complete,
    _cost
)
from celest.encounter.windows import generate_vtw
from celest.encounter._window_handling import OW, OWHandler
import math


class Request:
    """Request(location, deadline, duration, priority=1, quality=1, look_ang=None)

    Parameters
    ----------
    location : GroundPosition
        Location of the request.
    deadline : float
        Deadline of the request in Julian days.
    duration : float
        Duration of the request in seconds.
    priority : float, optional
        Priority of the request. The default is 1.
    quality : float, optional
        Quality of the request. The default is 1.
    look_ang : float, optional
        Look angle of the request. The default is None.

    Attributes
    ----------
    location : GroundPosition
    deadline : float
    duration : float
    priority : float
    quality : float
    look_ang : float
    vtws : np.ndarray
        The visible time windows associated with the request.
    is_scheduled : bool
        Whether the request is scheduled or not.
    schedule_idx : int
        The index of the request in the visible time window array.
    scheduled_start : float
        The scheduled start time of the request in Julian days.
    scheduled_duration : float
        The scheduled duration of the request in seconds.
    """

    def __init__(self, location, deadline, duration, priority=1, quality=1,
                 look_ang=None) -> None:

        self.location = location
        self.deadline = deadline
        self.duration = duration
        self.priority = priority
        self.quality = quality
        self.look_ang = look_ang

    def _add_vtw_info(self, vtws) -> None:

        self.vtws = vtws
        self.is_scheduled = False
        self.scheduled_idx = None
        self.scheduled_start = None
        self.scheduled_duration = None


class Schedule(ALNS):
    """Schedule(satellite, vis_threshold)

    This class allows a user to specify a series of requests and creates a
    near optimal schedule using an adaptive large neighborhood search
    metaheuristic.

    Parameters
    ----------
    satellite : Satellite
        The satellite participating in the encounter.
    vis_threshold : float
        The elevation angle required by the satellite to be visible ot a ground
        target.
    """

    def __init__(self, satellite, vis_threshold) -> None:

        self.satellite = satellite
        self.vis_threshold = vis_threshold
        self.requests = []

    def add_request(self, location, deadline, duration, priority=1,
                    look_ang=None) -> None:
        """Add an encounter request to the satellite.

        Parameters
        ----------
        location : GroundPosition
            Ground location participating in the encounter.
        deadline : float
            The request deadline in Julian days.
        duration : float
            Minimum encounter duration.
        priority : int, optional
            Encounter priority on a range from 1 to 10, by default 1.

            If all requests have the same priority, the algorithm will
            maximize the number of requests scheduled.
        look_ang : float, optional
            A specific encounter look angle, by default None.
        """

        rq = Request(location, deadline, duration, priority, look_ang)
        vtws = generate_vtw(self.satellite, location, self.vis_threshold, 1)
        rq._add_vtw_info(vtws.to_numpy())

        self.requests.append(rq)

    def generate(self, max_iter, annealing_coeff, react_factor) -> OWHandler:
        """Return a near optimal schedule.

        This method will implement the adaptive large neighborhood search
        metaheuristic to find a near optimal and feasible solution.

        Parameters
        ----------
        max_iter : int
            Maximum number of iterations.
        annealing_coeff : float
            Annealing coefficient for simulated annealing acceptance criterion.
        react_factor : float
            Controls how sensitive weights are to updates. Must be in the
            range [0, 1].

        Returns
        -------
        OWHandler
            A container holding the scheduled observation time windows.
        """

        super().__init__(_initial_solution(self.requests))

        self.add_cost_func(_cost)
        self.add_destroy_funcs(_REMOVAL_FUNCTIONS)
        self.add_repair_funcs(_INSERTION_FUNCTIONS)
        self.add_is_complete_func(_is_complete)

        na = len(self.requests)
        nr = math.ceil(na / 4)

        xb = self.solve(max_iter, annealing_coeff, react_factor, nr, na)

        soln = OWHandler()
        for request in xb:

            if request.is_scheduled:

                idx = request.scheduled_idx
                roll = request.vtws[idx].roll(request.scheduled_start)
                pitch = request.vtws[idx].pitch(request.scheduled_start)
                yaw = request.vtws[idx].yaw(request.scheduled_start)

                ow = OW(request.scheduled_start, request.scheduled_duration,
                        request.location, request.deadline, roll, pitch, yaw)
                soln._add_window(ow)

        return soln
