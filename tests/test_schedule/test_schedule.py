

from celest.encounter.groundposition import GroundPosition
from celest.encounter._window_handling import OWHandler
from celest.satellite.satellite import Satellite
from celest.schedule.insertion_operators import _INSERTION_FUNCTIONS
from celest.schedule.request_handler import RequestIndices
from celest.schedule.removal_operators import _REMOVAL_FUNCTIONS
from celest.schedule.schedule import Schedule
from celest.schedule.scheduling_utils import cost, is_complete, initialize_solution
from unittest import TestCase
import numpy as np


class TestSchedule(TestCase):

    def setUp(self):

        fname = 'Tests/test_data/coordinate_validation_long.txt'
        data = np.loadtxt(fname, usecols=[0, 11, 12, 13, 14, 15, 16], skiprows=1, max_rows=5000)
        julian, gcrs_pos, gcrs_vel = data[:, 0], data[:, 1:4], data[:, 4:]

        self.satellite = Satellite(gcrs_pos, gcrs_vel, "gcrs", julian, 2430000)
        self.toronto = GroundPosition(43.6532, -79.3832, 0.076)

    def test_schedule_raises_exception_if_null_satellite_is_passed(self):

        self.assertRaises(TypeError, Schedule, None, 0)

    def test_schedule_raises_exception_if_null_time_is_passed(self):

        self.assertRaises(TypeError, Schedule, self.satellite, None)

    def test_add_request(self):

        deadline = 2460467
        duration = 30
        priority = 1
        quality = 1
        look_ang = None

        schedule = Schedule(self.satellite, 10)
        schedule.add_request(self.toronto, deadline, duration, priority, quality, look_ang)

        self.assertEqual(len(schedule.request_handler.requests), 1)
        self.assertEqual(schedule.request_handler[0][RequestIndices.deadline], deadline)
        self.assertEqual(schedule.request_handler[0][RequestIndices.duration], duration)
        self.assertEqual(schedule.request_handler[0][RequestIndices.priority], priority)
        self.assertEqual(schedule.request_handler[0][RequestIndices.quality], quality)
        self.assertEqual(schedule.request_handler[0][RequestIndices.look_angle], look_ang)

    def test_generate_raises_exception_when_no_requests_are_added(self):

        schedule = Schedule(self.satellite, 10)
        self.assertRaises(Exception, schedule.generate, 100, 0.8, 0.5)

    def test_generate_returns_OWHandler_when_requests_are_added(self):

        schedule = Schedule(self.satellite, 10)
        schedule.add_request(self.toronto, 2460467, 30, 1, 1, None)
        self.assertIsInstance(schedule.generate(100, 0.8, 0.5), OWHandler)

    def test_generate_returns_improved_solutions(self):

        toronto = GroundPosition(43.65, -79.38, 0.076)
        north_bay = GroundPosition(46.31, -79.46, 0.193)
        sudbury = GroundPosition(46.49, -80.99, 0.348)
        ottawa = GroundPosition(45.42, -75.70, 0.070)
        kingston = GroundPosition(44.23, -76.49, 0.093)
        niagara_falls = GroundPosition(43.09, -79.08, 0.099)
        london = GroundPosition(42.98, -81.24, 0.251)
        mississauga = GroundPosition(43.59, -79.64, 0.156)
        timmins = GroundPosition(48.48, -81.33, 0.295)
        tobermory = GroundPosition(45.25, -81.66, 0.271)

        schedule = Schedule(self.satellite, 10)

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

        best_solution = schedule.solve(100, 0.8, 0.5, remove_per_iteration, insert_per_iteration)
        best_solution_cost = cost(best_solution)

        self.assertLess(best_solution_cost, initial_solution_cost)
