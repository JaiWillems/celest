

from celest.schedule.insertion_operators import (
    does_conflict_exists,
    greedy_insertion,
    image_quality,
    image_quality_is_met,
    _look_angle_time,
    minimum_conflict_insertion,
    minimum_opportunity_insertion
)
from celest.schedule.request_handler import RequestHandler
from polare import Stroke
from tests.test_schedule.request_handler_testing_utils import initialize_request_list
from unittest import TestCase
import numpy as np


class TestInsertionOperators(TestCase):

    def setUp(self):

        self.request_handler = RequestHandler()

        request_list = initialize_request_list()
        for request in request_list:
            self.request_handler.add_request(*request)

    def test_greedy_insert_one_item(self):

        greedy_insertion(self.request_handler, 1)
        self.assertEqual(self.request_handler.number_of_scheduled_requests, 1)
        self.assertTrue(self.request_handler[0][0])
        self.assertEqual(self.request_handler[0][7], 7)

    def test_greedy_insert_two_items(self):

        greedy_insertion(self.request_handler, 2)
        self.assertEqual(self.request_handler.number_of_scheduled_requests, 2)
        self.assertTrue(self.request_handler[0][0])
        self.assertEqual(self.request_handler[0][7], 7)
        self.assertTrue(self.request_handler[1][0])
        self.assertEqual(self.request_handler[1][7], 7)

    def test_greedy_insert_all_items(self):

        greedy_insertion(self.request_handler, 10)
        self.assertEqual(self.request_handler.number_of_scheduled_requests, 7)

    def test_look_angle_time_with_root(self):

        time = np.linspace(0, 1, 100)
        look_angle_stroke = Stroke(time, time)
        look_angle_time = _look_angle_time(look_angle_stroke, 0.5, 0, 1, )

        self.assertAlmostEqual(look_angle_time, 0.5, delta=1e-6)

    def test_look_angle_time_with_no_root(self):

        time = np.linspace(0, 1, 100)
        look_angle_stroke = Stroke(time, time)

        self.assertIsNone(_look_angle_time(look_angle_stroke, 2, 0, 1))

    def test_image_quality_with_nadir_time(self):

        vtw = self.request_handler[0][10][0]
        nadir_time = (vtw.rise_time + vtw.set_time) / 2
        self.assertEqual(image_quality(vtw, nadir_time), 10)

    def test_image_quality_with_start_time(self):

        vtw = self.request_handler[0][10][0]
        start_time = vtw.rise_time
        self.assertEqual(image_quality(vtw, start_time), 1)

    def test_image_quality_with_set_time(self):

        vtw = self.request_handler[0][10][0]
        start_time = vtw.set_time
        self.assertEqual(image_quality(vtw, start_time), 1)

    def test_image_quality_is_met_with_nadir_time(self):

        vtw = self.request_handler[0][10][0]
        nadir_time = (vtw.rise_time + vtw.set_time) / 2
        self.assertTrue(image_quality_is_met(vtw, nadir_time, 9))
        self.assertTrue(image_quality_is_met(vtw, nadir_time, 10))

    def test_image_quality_is_met_with_start_time(self):

        vtw = self.request_handler[0][10][0]
        start_time = vtw.rise_time
        self.assertTrue(image_quality_is_met(vtw, start_time, 1))
        self.assertFalse(image_quality_is_met(vtw, start_time, 2))

    def test_conflict_exists_with_no_conflicts(self):

        self._schedule_first_request()

        rise_time = self.request_handler[0][10][0].rise_time
        set_time = self.request_handler[0][10][0].set_time

        duration_in_days = set_time - rise_time
        start = rise_time + 2 * duration_in_days

        self.assertFalse(does_conflict_exists(self.request_handler, start, duration_in_days))

    def _schedule_first_request(self):

        rise_time = self.request_handler[0][10][0].rise_time
        set_time = self.request_handler[0][10][0].set_time

        self.request_handler.schedule_request(0, 0, rise_time, 86400 * (set_time - rise_time))

    def test_conflict_exists_with_window_around_scheduled_window(self):

        self._schedule_first_request()

        rise_time = self.request_handler[0][10][0].rise_time
        set_time = self.request_handler[0][10][0].set_time

        duration_in_days = 2 * (set_time - rise_time)
        start = rise_time - duration_in_days / 8

        self.assertTrue(does_conflict_exists(self.request_handler, start, duration_in_days))

    def test_conflict_exists_with_start_time_in_scheduled_window(self):

        self._schedule_first_request()

        rise_time = self.request_handler[0][10][0].rise_time
        set_time = self.request_handler[0][10][0].set_time

        duration_in_days = set_time - rise_time
        start = rise_time + duration_in_days / 2

        self.assertTrue(does_conflict_exists(self.request_handler, start, duration_in_days))

    def test_conflict_exists_with_end_time_in_scheduled_window(self):

        self._schedule_first_request()

        rise_time = self.request_handler[0][10][0].rise_time
        set_time = self.request_handler[0][10][0].set_time

        duration_in_days = set_time - rise_time
        start = set_time - duration_in_days / 2

        self.assertTrue(does_conflict_exists(self.request_handler, start, duration_in_days))

    def test_conflict_exists_with_window_in_scheduled_window(self):

        self._schedule_first_request()

        rise_time = self.request_handler[0][10][0].rise_time
        set_time = self.request_handler[0][10][0].set_time

        duration_in_days = (set_time - rise_time) / 2
        start = rise_time + duration_in_days / 2

        self.assertTrue(does_conflict_exists(self.request_handler, start, duration_in_days))

    def test_minimum_opportunity_insert_one_item(self):

        minimum_opportunity_insertion(self.request_handler, 1)
        self.assertEqual(self.request_handler.number_of_scheduled_requests, 1)
        self.assertTrue(self.request_handler[0][0])
        self.assertEqual(len(self.request_handler[0][10]), 20)

    def test_minimum_opportunity_insert_two_items(self):

        minimum_opportunity_insertion(self.request_handler, 2)
        self.assertEqual(self.request_handler.number_of_scheduled_requests, 2)
        self.assertTrue(self.request_handler[0][0])
        self.assertEqual(len(self.request_handler[0][10]), 20)
        self.assertTrue(self.request_handler[1][0])
        self.assertEqual(len(self.request_handler[1][10]), 19)

    def test_minimum_opportunity_insert_all_items(self):

        minimum_opportunity_insertion(self.request_handler, 10)
        self.assertEqual(self.request_handler.number_of_scheduled_requests, 8)

    def test_minimum_conflict_insert_one_item(self):

        minimum_conflict_insertion(self.request_handler, 1)
        self.assertEqual(self.request_handler.number_of_scheduled_requests, 1)
        self.assertTrue(self.request_handler[0][0])

    def test_minimum_conflict_insert_two_items(self):

        minimum_conflict_insertion(self.request_handler, 2)
        self.assertEqual(self.request_handler.number_of_scheduled_requests, 2)
        self.assertTrue(self.request_handler[0][0])
        self.assertTrue(self.request_handler[1][0])

    def test_minimum_conflict_insert_all_items(self):

        minimum_conflict_insertion(self.request_handler, 10)
        self.assertEqual(self.request_handler.number_of_scheduled_requests, 4)
