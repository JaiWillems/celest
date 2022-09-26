

from celest.schedule.insertion_operators import (
    does_conflict_exist,
    greedy_insertion,
    image_quality,
    is_image_quality_met,
    _look_angle_time,
    minimum_conflict_insertion,
    minimum_opportunity_insertion
)
from celest.schedule.request_handler import RequestHandler, RequestIndices
from celest.units.quantity import Quantity
from celest import units as u
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
        self.assertTrue(self.request_handler[0][RequestIndices.is_scheduled])
        self.assertEqual(self.request_handler[0][RequestIndices.priority], 7)

    def test_greedy_insert_two_items(self):
        greedy_insertion(self.request_handler, 2)
        self.assertEqual(self.request_handler.number_of_scheduled_requests, 2)
        self.assertTrue(self.request_handler[0][RequestIndices.is_scheduled])
        self.assertEqual(self.request_handler[0][RequestIndices.priority], 7)
        self.assertTrue(self.request_handler[1][RequestIndices.is_scheduled])
        self.assertEqual(self.request_handler[1][RequestIndices.priority], 7)

    def test_greedy_insert_all_items(self):
        greedy_insertion(self.request_handler, 10)
        self.assertEqual(self.request_handler.number_of_scheduled_requests, 7)

    def test_look_angle_time_with_root(self):
        look_angle = Quantity(np.linspace(0, 1, 100), u.deg)
        desired_look_angle = Quantity(0.5, u.deg)
        julian = Quantity(np.linspace(0, 1, 100), u.jd2000)
        start = Quantity(0, u.jd2000)
        end = Quantity(1, u.jd2000)
        look_angle_time = _look_angle_time(
            look_angle,
            desired_look_angle,
            julian,
            start,
            end
        )

        self.assertAlmostEqual(look_angle_time.data, 0.5, delta=1e-2)

    def test_look_angle_time_with_no_root(self):
        look_angle = Quantity(np.linspace(0, 1, 100), u.deg)
        desired_look_angle = Quantity(2, u.deg)
        julian = Quantity(np.linspace(0, 1, 100), u.jd2000)
        start = Quantity(0, u.jd2000)
        end = Quantity(1, u.jd2000)
        look_angle_time = _look_angle_time(
            look_angle,
            desired_look_angle,
            julian,
            start,
            end
        )

        self.assertIsNone(look_angle_time)

    def test_image_quality_with_nadir_time(self):
        vtw = self.request_handler[0][RequestIndices.vtw_list][0]
        nadir_time = (vtw.rise_time + vtw.set_time) / 2
        self.assertEqual(image_quality(vtw, nadir_time), 10)

    def test_image_quality_with_start_time(self):
        vtw = self.request_handler[0][RequestIndices.vtw_list][0]
        start_time = vtw.rise_time
        self.assertEqual(image_quality(vtw, start_time), 1)

    def test_image_quality_with_set_time(self):
        vtw = self.request_handler[0][RequestIndices.vtw_list][0]
        start_time = vtw.set_time
        self.assertEqual(image_quality(vtw, start_time), 1)

    def test_image_quality_is_met_with_nadir_time(self):
        vtw = self.request_handler[0][RequestIndices.vtw_list][0]
        nadir_time = (vtw.rise_time + vtw.set_time) / 2
        self.assertTrue(is_image_quality_met(vtw, nadir_time, 9))
        self.assertTrue(is_image_quality_met(vtw, nadir_time, 10))

    def test_image_quality_is_met_with_start_time(self):
        vtw = self.request_handler[0][RequestIndices.vtw_list][0]
        start_time = vtw.rise_time
        self.assertTrue(is_image_quality_met(vtw, start_time, 1))
        self.assertFalse(is_image_quality_met(vtw, start_time, 2))

    def test_conflict_exists_with_no_conflicts(self):
        self._schedule_first_request()

        vtw = self.request_handler[0][RequestIndices.vtw_list][0]
        rise_time = vtw.rise_time
        set_time = vtw.set_time

        duration_in_days = set_time - rise_time
        start = rise_time + 2 * duration_in_days

        self.assertFalse(does_conflict_exist(self.request_handler, start, duration_in_days))

    def _schedule_first_request(self):
        vtw = self.request_handler[0][RequestIndices.vtw_list][0]
        rise_time = vtw.rise_time
        set_time = vtw.set_time
        duration = Quantity(86400 * (set_time.to(u.jd2000) - rise_time.to(u.jd2000)), u.s)

        self.request_handler.schedule_request(0, 0, rise_time, duration)

    def test_conflict_exists_with_window_around_scheduled_window(self):
        self._schedule_first_request()

        vtw = self.request_handler[0][RequestIndices.vtw_list][0]
        rise_time = vtw.rise_time
        set_time = vtw.set_time

        duration_in_days = 2 * (set_time - rise_time)
        start = rise_time - duration_in_days / 8

        self.assertTrue(does_conflict_exist(self.request_handler, start, duration_in_days))

    def test_conflict_exists_with_start_time_in_scheduled_window(self):
        self._schedule_first_request()

        vtw = self.request_handler[0][RequestIndices.vtw_list][0]
        rise_time = vtw.rise_time
        set_time = vtw.set_time

        duration_in_days = set_time - rise_time
        start = rise_time + duration_in_days / 2

        self.assertTrue(does_conflict_exist(self.request_handler, start, duration_in_days))

    def test_conflict_exists_with_end_time_in_scheduled_window(self):
        self._schedule_first_request()

        vtw = self.request_handler[0][RequestIndices.vtw_list][0]
        rise_time = vtw.rise_time
        set_time = vtw.set_time

        duration_in_days = set_time - rise_time
        start = set_time - duration_in_days / 2

        self.assertTrue(does_conflict_exist(self.request_handler, start, duration_in_days))

    def test_conflict_exists_with_window_in_scheduled_window(self):
        self._schedule_first_request()

        vtw = self.request_handler[0][RequestIndices.vtw_list][0]
        rise_time = vtw.rise_time
        set_time = vtw.set_time

        duration_in_days = (set_time - rise_time) / 2
        start = rise_time + duration_in_days / 2

        self.assertTrue(does_conflict_exist(self.request_handler, start, duration_in_days))

    def test_minimum_opportunity_insert_one_item(self):
        minimum_opportunity_insertion(self.request_handler, 1)
        self.assertEqual(self.request_handler.number_of_scheduled_requests, 1)
        self.assertTrue(self.request_handler[0][RequestIndices.is_scheduled])
        self.assertEqual(len(self.request_handler[0][RequestIndices.vtw_list]), 18)

    def test_minimum_opportunity_insert_two_items(self):
        minimum_opportunity_insertion(self.request_handler, 2)
        self.assertEqual(self.request_handler.number_of_scheduled_requests, 2)
        self.assertTrue(self.request_handler[0][RequestIndices.is_scheduled])
        self.assertEqual(len(self.request_handler[0][RequestIndices.vtw_list]), 18)
        self.assertTrue(self.request_handler[1][RequestIndices.is_scheduled])
        self.assertEqual(len(self.request_handler[1][RequestIndices.vtw_list]), 18)

    def test_minimum_opportunity_insert_all_items(self):
        minimum_opportunity_insertion(self.request_handler, 10)
        self.assertEqual(self.request_handler.number_of_scheduled_requests, 8)

    def test_minimum_conflict_insert_one_item(self):
        minimum_conflict_insertion(self.request_handler, 1)
        self.assertEqual(self.request_handler.number_of_scheduled_requests, 1)
        self.assertTrue(self.request_handler[0][RequestIndices.is_scheduled])

    def test_minimum_conflict_insert_two_items(self):
        minimum_conflict_insertion(self.request_handler, 2)
        self.assertEqual(self.request_handler.number_of_scheduled_requests, 2)
        self.assertTrue(self.request_handler[0][RequestIndices.is_scheduled])
        self.assertTrue(self.request_handler[1][RequestIndices.is_scheduled])

    def test_minimum_conflict_insert_all_items(self):
        minimum_conflict_insertion(self.request_handler, 10)
        self.assertEqual(self.request_handler.number_of_scheduled_requests, 4)
