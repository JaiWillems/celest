

from celest.coordinates.frames.gcrs import GCRS
from celest.coordinates.ground_location import GroundLocation
from celest.encounter.window_handling import WindowCollection
from celest.satellite import Satellite
from celest.schedule.insertion_operators import _INSERTION_FUNCTIONS
from celest.schedule.removal_operators import _REMOVAL_FUNCTIONS
from celest.schedule.request_handler import RequestIndices
from celest.schedule.scheduler import Scheduler
from celest.schedule.scheduling_utils import (
    cost,
    initialize_solution,
    is_complete
)
from celest import units as u
from unittest import TestCase
import numpy as np


class TestSchedule(TestCase):

    def setUp(self):
        fname = 'Tests/test_data/coordinate_validation_long.txt'
        data = np.loadtxt(
            fname,
            usecols=[0, 11, 12, 13, 14, 15, 16],
            skiprows=1,
            max_rows=5000
        )
        julian = data[:, 0] + 2430000
        gcrs_pos = data[:, 1:4]
        gcrs_vel = data[:, 4:]

        gcrs_position = GCRS(
            julian,
            gcrs_pos[:, 0],
            gcrs_pos[:, 1],
            gcrs_pos[:, 2],
            u.km
        )
        gcrs_velocity = GCRS(
            julian,
            gcrs_vel[:, 0],
            gcrs_vel[:, 1],
            gcrs_vel[:, 2],
            u.m / u.s
        )

        self.satellite = Satellite(gcrs_position, gcrs_velocity)
        self.toronto = GroundLocation(43.6532, -79.3832, 0.076, u.deg, u.km)

    def test_schedule_raises_exception_if_null_satellite_is_passed(self):
        self.assertRaises(TypeError, Scheduler, None, 0)

    def test_schedule_raises_exception_if_null_time_is_passed(self):
        self.assertRaises(TypeError, Scheduler, self.satellite, None)

    def test_add_request(self):
        deadline = 2460467
        duration = 30
        priority = 1
        quality = 1
        look_ang = None

        schedule = Scheduler(self.satellite, 10)
        schedule.add_request(
            self.toronto,
            deadline,
            duration,
            priority,
            quality,
            look_ang
        )

        self.assertEqual(
            len(schedule.request_handler.requests),
            1
        )
        self.assertEqual(
            schedule.request_handler[0][RequestIndices.deadline].data,
            deadline
        )
        self.assertEqual(
            schedule.request_handler[0][RequestIndices.duration].data,
            duration
        )
        self.assertEqual(
            schedule.request_handler[0][RequestIndices.priority],
            priority
        )
        self.assertEqual(
            schedule.request_handler[0][RequestIndices.quality],
            quality
        )
        self.assertEqual(
            schedule.request_handler[0][RequestIndices.look_angle],
            look_ang
        )

    def test_generate_raises_exception_when_no_requests_are_added(self):
        schedule = Scheduler(self.satellite, 10)
        self.assertRaises(Exception, schedule.generate, 100, 0.8, 0.5)

    def test_generate_returns_WindowHandler_when_requests_are_added(self):
        schedule = Scheduler(self.satellite, 10)
        schedule.add_request(self.toronto, 2460467, 30, 1, 1, None)
        self.assertIsInstance(
            schedule.generate(100, 0.8, 0.5),
            WindowCollection
        )

    def test_generate_returns_improved_solutions(self):
        toronto = GroundLocation(43.65, -79.38, 0.076, u.deg, u.km)
        north_bay = GroundLocation(46.31, -79.46, 0.193, u.deg, u.km)
        sudbury = GroundLocation(46.49, -80.99, 0.348, u.deg, u.km)
        ottawa = GroundLocation(45.42, -75.70, 0.070, u.deg, u.km)
        kingston = GroundLocation(44.23, -76.49, 0.093, u.deg, u.km)
        niagara_falls = GroundLocation(43.09, -79.08, 0.099, u.deg, u.km)
        london = GroundLocation(42.98, -81.24, 0.251, u.deg, u.km)
        mississauga = GroundLocation(43.59, -79.64, 0.156, u.deg, u.km)
        timmins = GroundLocation(48.48, -81.33, 0.295, u.deg, u.km)
        tobermory = GroundLocation(45.25, -81.66, 0.271, u.deg, u.km)

        schedule = Scheduler(self.satellite, 10)

        schedule.add_request(toronto, 2460467, 30, 1, 1, None)
        schedule.add_request(north_bay, 2460467, 30, 1, 1, None)
        schedule.add_request(sudbury, 2460467, 30, 4, 1, None)
        schedule.add_request(ottawa, 2460467, 30, 2, 1, None)
        schedule.add_request(kingston, 2460467, 30, 7, 1, None)
        schedule.add_request(niagara_falls, 2460467, 30, 3, 1, None)
        schedule.add_request(london, 2460467, 30, 4, 1, None)
        schedule.add_request(mississauga, 2460467, 30, 5, 1, None)
        schedule.add_request(timmins, 2460467, 30, 1, 1, None)
        schedule.add_request(tobermory, 2460467, 30, 7, 1, None)

        initialize_solution(schedule.request_handler)
        initial_solution_cost = cost(schedule.request_handler)

        schedule.initial_solution = schedule.request_handler

        schedule.add_cost_function(cost)
        schedule.add_destroy_functions(_REMOVAL_FUNCTIONS)
        schedule.add_repair_functions(_INSERTION_FUNCTIONS)
        schedule.add_completeness_function(is_complete)

        insert_per_iteration = schedule.request_handler.number_of_requests
        remove_per_iteration = np.ceil(insert_per_iteration / 4)

        best_solution = schedule.solve(
            100,
            0.8,
            0.5,
            remove_per_iteration,
            insert_per_iteration
        )
        best_solution_cost = cost(best_solution)

        self.assertLess(best_solution_cost, initial_solution_cost)
