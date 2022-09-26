

from celest.coordinates.frames.attitude import Attitude
from celest.coordinates.ground_location import GroundLocation
from celest.satellite import Satellite
from celest.schedule.alns import ALNS
from celest.schedule.insertion_operators import _INSERTION_FUNCTIONS
from celest.schedule.removal_operators import _REMOVAL_FUNCTIONS
from celest.schedule.request_handler import RequestHandler, RequestIndices
from celest.schedule.scheduling_utils import (
    cost,
    initialize_solution,
    is_complete
)
from celest.encounter.window_generator import generate_vtws, Lighting
from celest.encounter.window_handling import ObservationWindow, WindowCollection
from celest.units import Quantity
from celest import units as u
import math
import numpy as np


class Scheduler(ALNS):
    """Scheduler(satellite, vis_threshold)

    Creates a scheduling object for the given satellite and visibility
    threshold.

    This class implements an adaptive large neighborhood search metaheuristic
    for scheduling satellite observations.

    Parameters
    ----------
    satellite : Satellite
        The satellite executing the scheduled observations.
    vis_threshold : float
        Visibility threshold in degrees.

        The visibility threshold is the minimum elevation angle of the satellite
        as seen from `location` where the satellite will be in visual range of
        `location`.

    Methods
    -------
    add_request(location, deadline, duration, priority, quality, look_ang, lighting)
        Adds an encounter request to the schedule.
    generate(max_iter, annealing_coeff, react_factor)
        Generates a schedule for the given satellite.
    """

    def __init__(self, satellite: Satellite, vis_threshold: float) -> None:
        """Encounter scheduler.

        This class implements an adaptive large neighborhood search
        metaheuristic for scheduling satellite observations.

        Parameters
        ----------
        satellite : Satellite
            The satellite executing the scheduled observations.
        vis_threshold : float
            Visibility threshold in degrees.

            The visibility threshold is the minimum elevation angle of the
            satellite as seen from `location` where the satellite will be in
            visual range of `location`.
        """

        if not isinstance(satellite, Satellite):
            raise TypeError("satellite must be of type Satellite.")
        if not isinstance(vis_threshold, (float, int)):
            raise TypeError("vis_threshold must be a float or integer.")

        self.satellite = satellite
        self.visibility_threshold = vis_threshold
        self.request_handler = RequestHandler()

    def add_request(self, location: GroundLocation, deadline: float, duration:
                    float, priority: int, quality: int, look_ang: float=None,
                    lighting: Lighting=Lighting.DAYTIME) -> None:
        """Add a request to be scheduled.
        
        Parameters
        ----------
        location : GroundLocation
            The location associated with the encounter.
        deadline : float
            The deadline for the request to be scheduled by in the jd2000 epoch.
        duration : float
            The desired duration of the encounter in seconds.
        priority: int
            The priority of the request on a scale of 1 to 10 where higher
            priority requests are scheduled first.
        quality : int
            The quality of the request on a scale of 1 to 10. A higher quality
            request will occur with a smaller look_angle. For imaging requests,
            the imaging will be taken closer to nadir for a higher quality.
        look_ang : float, optional
            The desired look angle of the encounter in degrees, default is None.
        lighting : Lighting, optional
            The lighting condition of the encounter, default is DAYTIME.
        """

        vtws = generate_vtws(
            self.satellite,
            location,
            self.visibility_threshold,
            lighting
        )
        self.request_handler.add_request(
            location,
            Quantity(deadline, u.jd2000),
            Quantity(duration, u.s),
            priority,
            quality,
            look_ang,
            vtws
        )

    def generate(self, max_iter: int, annealing_coeff: float,
                 react_factor: float) -> WindowCollection:
        """Generate the schedule for the given satellite and requests.

        Parameters
        ----------
        max_iter : int
            The maximum number of iterations to perform in search for an
            improved solution.
        annealing_coeff : float
            The annealing coefficient between 0 and 1 for the simulated
            annealing process.
        react_factor : float
            Decay parameter within the range [0, 1].

            The decay parameter dictates the effect of previous iterations on
            the current iteration. A `react_factor` of 1 will cause the
            algorithm to ignore previous iterations and a `react_factor` of 0
            will cause the algorithm to only consider the previous iterations.

        Returns
        -------
        WindowCollection
            A container holding the scheduled windows.
        """

        if not self.request_handler.number_of_requests:
            raise Exception("No requests to schedule.")

        initialize_solution(self.request_handler)
        super().__init__(self.request_handler)

        self.add_cost_function(cost)
        self.add_destroy_functions(_REMOVAL_FUNCTIONS)
        self.add_repair_functions(_INSERTION_FUNCTIONS)
        self.add_completeness_function(is_complete)

        insert_per_iteration = self.request_handler.number_of_requests
        remove_per_iteration = math.ceil(insert_per_iteration / 4)

        best_solution = self.solve(
            max_iter,
            annealing_coeff,
            react_factor,
            remove_per_iteration,
            insert_per_iteration
        )

        return self._generate_window_handler_from_request_list(best_solution)

    def _generate_window_handler_from_request_list(self, request_list:
                                                   RequestHandler) -> WindowCollection:
        window_handler = WindowCollection()
        for request in request_list:
            if request[RequestIndices.is_scheduled]:

                idx = request[RequestIndices.vtw_index]
                attitude = request[RequestIndices.vtw_list][idx].attitude
                attitude_idx = np.where(
                    attitude.time.to(u.jd2000).data ==
                    request[RequestIndices.scheduled_start_time]
                )[0]
                time = attitude.time.data[attitude_idx]
                roll = attitude.roll.data[attitude_idx]
                pitch = attitude.pitch.data[attitude_idx]
                yaw = attitude.yaw.data[attitude_idx]
                unit = attitude.roll.unit
                location = attitude.location
                new_attitude = Attitude(time, roll, pitch, yaw, unit, location)

                ow = ObservationWindow(
                    request[RequestIndices.scheduled_start_time],
                    request[RequestIndices.scheduled_duration],
                    request[RequestIndices.location],
                    new_attitude
                )
                window_handler.add_window(ow)

        return window_handler
