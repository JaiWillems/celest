from celest.coordinates import Attitude
from celest.satellite import Satellite
from celest.schedule.alns import ALNS
from celest.schedule.insertion_operators import _INSERTION_FUNCTIONS
from celest.schedule.removal_operators import _REMOVAL_FUNCTIONS
from celest.schedule.request_handler import RequestHandler
from celest.schedule.scheduling_utils import (
    cost,
    initialize_solution,
    is_complete
)
from celest.encounter.window_generator import generate_vtw, Lighting
from celest.encounter.window_handling import ObservationWindow, WindowHandler
import math
import numpy as np


class Schedule(ALNS):

    def __init__(self, satellite, visibility_threshold):

        if not isinstance(satellite, Satellite):
            raise TypeError("satellite must be of type Satellite.")
        if not isinstance(visibility_threshold, (float, int)):
            raise TypeError("visibility_threshold must be a float or integer.")

        self.satellite = satellite
        self.visibility_threshold = visibility_threshold
        self.request_handler = RequestHandler()

    def add_request(self, location, deadline, duration, priority, quality, look_ang):

        # TODO: Make the inputs into Quantity objects.
        # TODO: Add a lighting parameter.
        vtws = generate_vtw(self.satellite, location, self.visibility_threshold, Lighting.DAYTIME)
        self.request_handler.add_request(location, deadline, duration, priority, quality, look_ang, vtws)

    def generate(self, max_iter, annealing_coeff, react_factor):

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

        best_solution = self.solve(max_iter, annealing_coeff, react_factor,
                                   remove_per_iteration, insert_per_iteration)

        return self._generate_OWHandler_from_request_list(best_solution)

    def _generate_OWHandler_from_request_list(self, request_list):

        ow_handler = WindowHandler()
        for request in request_list:

            # TODO: Add in the RequestIndices rather than hard coded values.
            if request[0]:

                idx = request[1]
                attitude = request[10][idx].attitude
                attitude_idx = np.where(attitude.time.data == request[2])[0]
                time = attitude.time.data[attitude_idx]
                roll = attitude.roll.data[attitude_idx]
                pitch = attitude.pitch.data[attitude_idx]
                yaw = attitude.yaw.data[attitude_idx]
                unit = attitude.roll.unit
                location = attitude.location
                new_attitude = Attitude(time, roll, pitch, yaw, unit, location)

                ow = ObservationWindow(request[2], request[3], request[5],
                                       request[4], new_attitude)
                ow_handler.add_window(ow)

        return ow_handler
