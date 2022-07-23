

from celest.schedule.request_handler import RequestHandler, RequestIndices
from tests.test_schedule.request_handler_testing_utils import initialize_request_list
from unittest import TestCase


class TestRequestHandler(TestCase):

    def setUp(self) -> None:
        self.request_handler = RequestHandler()
        for request in initialize_request_list():
            self.request_handler.add_request(*request)

    def test_getitem(self):
        self.assertListEqual(self.request_handler[0],
                             self.request_handler.requests[0])

    def test_request_handler_iteration(self):
        for request in self.request_handler:
            self.assertIsInstance(request, list)

    def test_add_request(self):
        request_list = initialize_request_list()

        self.assertListEqual(self.request_handler[0][4:-1], request_list[0][:-1])
        self.assertListEqual(self.request_handler[1][4:-1], request_list[1][:-1])
        self.assertListEqual(self.request_handler[2][4:-1], request_list[2][:-1])
        self.assertListEqual(self.request_handler[3][4:-1], request_list[3][:-1])
        self.assertListEqual(self.request_handler[4][4:-1], request_list[4][:-1])
        self.assertListEqual(self.request_handler[5][4:-1], request_list[5][:-1])
        self.assertListEqual(self.request_handler[6][4:-1], request_list[6][:-1])
        self.assertListEqual(self.request_handler[7][4:-1], request_list[7][:-1])
        self.assertListEqual(self.request_handler[8][4:-1], request_list[8][:-1])
        self.assertListEqual(self.request_handler[9][4:-1], request_list[9][:-1])

        self.assertEqual(self.request_handler.number_of_requests, len(request_list))

    def test_schedule_request(self):
        self.request_handler.schedule_request(0, 1, 2, 3)

        self.assertTrue(self.request_handler[0][RequestIndices.is_scheduled])
        self.assertEqual(self.request_handler[0][RequestIndices.vtw_index], 1)
        self.assertEqual(self.request_handler[0][RequestIndices.scheduled_start_time], 2)
        self.assertEqual(self.request_handler[0][RequestIndices.scheduled_duration], 3)

    def test_unschedule_all_requests(self):
        self.request_handler.schedule_request(0, 1, 2, 3)
        self.request_handler.schedule_request(1, 1, 2, 3)
        self.request_handler.unschedule_all_requests()

        self.assertFalse(self.request_handler[0][RequestIndices.is_scheduled])
        self.assertFalse(self.request_handler[1][RequestIndices.is_scheduled])
        self.assertFalse(self.request_handler[2][RequestIndices.is_scheduled])
        self.assertFalse(self.request_handler[3][RequestIndices.is_scheduled])
        self.assertFalse(self.request_handler[4][RequestIndices.is_scheduled])
        self.assertFalse(self.request_handler[5][RequestIndices.is_scheduled])
        self.assertFalse(self.request_handler[6][RequestIndices.is_scheduled])
        self.assertFalse(self.request_handler[7][RequestIndices.is_scheduled])
        self.assertFalse(self.request_handler[8][RequestIndices.is_scheduled])
        self.assertFalse(self.request_handler[9][RequestIndices.is_scheduled])

    def test_unschedule_request(self):
        self.request_handler.schedule_request(0, 1, 2, 3)
        self.request_handler.unschedule_request(0)

        self.assertFalse(self.request_handler[0][RequestIndices.is_scheduled])
        self.assertIsNone(self.request_handler[0][RequestIndices.vtw_index])
        self.assertIsNone(self.request_handler[0][RequestIndices.scheduled_start_time])
        self.assertIsNone(self.request_handler[0][RequestIndices.scheduled_duration])

    def test_is_request_scheduled(self):
        self.assertFalse(self.request_handler.is_request_scheduled(0))
        self.request_handler.schedule_request(0, 1, 2, 3)
        self.assertTrue(self.request_handler.is_request_scheduled(0))

    def test_is_request_scheduled(self):
        self.assertFalse(self.request_handler.is_request_scheduled(0))
        self.request_handler.schedule_request(0, 1, 2, 3)
        self.assertTrue(self.request_handler.is_request_scheduled(0))

    def test_vtw_index(self):
        self.assertIsNone(self.request_handler.vtw_index(0))
        self.request_handler.schedule_request(0, 1, 2, 3)
        self.assertEqual(self.request_handler.vtw_index(0), 1)

    def test_start_time(self):
        self.assertIsNone(self.request_handler.start_time(0))
        self.request_handler.schedule_request(0, 1, 2, 3)
        self.assertEqual(self.request_handler.start_time(0), 2)

    def test_scheduled_duration(self):
        self.assertIsNone(self.request_handler.scheduled_duration(0))
        self.request_handler.schedule_request(0, 1, 2, 3)
        self.assertEqual(self.request_handler.scheduled_duration(0), 3)

    def test_deadline(self):
        test_deadline = self.request_handler.deadline(0)
        actual_deadline = self.request_handler[0][RequestIndices.deadline]
        self.assertEqual(test_deadline, actual_deadline)

    def test_duration(self):
        test_duration = self.request_handler.duration(0)
        actual_duration = self.request_handler[0][RequestIndices.duration]
        self.assertEqual(test_duration, actual_duration)

    def test_priority(self):
        test_priority = self.request_handler.priority(0)
        actual_priority = self.request_handler[0][RequestIndices.priority]
        self.assertEqual(test_priority, actual_priority)

    def test_image_quality(self):
        test_image_quality = self.request_handler.image_quality(0)
        actual_image_quality = self.request_handler[0][RequestIndices.quality]
        self.assertEqual(test_image_quality, actual_image_quality)

    def test_look_angle(self):
        test_look_angle = self.request_handler.look_angle(0)
        actual_look_angle = self.request_handler[0][RequestIndices.look_angle]
        self.assertEqual(test_look_angle, actual_look_angle)

    def test_vtw_list(self):
        test_vtw_list = self.request_handler.vtw_list(0)
        actual_vtw_list = self.request_handler[0][RequestIndices.vtw_list]
        self.assertEqual(test_vtw_list, actual_vtw_list)

    def test_number_of_scheduled_items(self):
        self.request_handler.schedule_request(0, 1, 2, 3)
        self.request_handler.schedule_request(1, 1, 2, 3)
        self.assertEqual(self.request_handler.number_of_scheduled_requests, 2)

    def test_sort_by_decreasing_priority(self):
        self.request_handler.sort_by_decreasing_priority()

        self.assertEqual(self.request_handler[0][RequestIndices.priority], 7)
        self.assertEqual(self.request_handler[1][RequestIndices.priority], 7)
        self.assertEqual(self.request_handler[2][RequestIndices.priority], 5)
        self.assertEqual(self.request_handler[3][RequestIndices.priority], 4)
        self.assertEqual(self.request_handler[4][RequestIndices.priority], 4)
        self.assertEqual(self.request_handler[5][RequestIndices.priority], 3)
        self.assertEqual(self.request_handler[6][RequestIndices.priority], 2)
        self.assertEqual(self.request_handler[7][RequestIndices.priority], 1)
        self.assertEqual(self.request_handler[8][RequestIndices.priority], 1)
        self.assertEqual(self.request_handler[9][RequestIndices.priority], 1)

    def test_sort_by_decreasing_opportunity(self):
        self.request_handler.sort_by_decreasing_opportunity()

        self.assertEqual(len(self.request_handler[0][RequestIndices.vtw_list]), 20)
        self.assertEqual(len(self.request_handler[1][RequestIndices.vtw_list]), 19)
        self.assertEqual(len(self.request_handler[2][RequestIndices.vtw_list]), 19)
        self.assertEqual(len(self.request_handler[3][RequestIndices.vtw_list]), 19)
        self.assertEqual(len(self.request_handler[4][RequestIndices.vtw_list]), 19)
        self.assertEqual(len(self.request_handler[5][RequestIndices.vtw_list]), 19)
        self.assertEqual(len(self.request_handler[6][RequestIndices.vtw_list]), 18)
        self.assertEqual(len(self.request_handler[7][RequestIndices.vtw_list]), 18)
        self.assertEqual(len(self.request_handler[8][RequestIndices.vtw_list]), 18)
        self.assertEqual(len(self.request_handler[9][RequestIndices.vtw_list]), 18)

    def test_sort_vtws_by_increasing_conflict_degree(self):
        self.request_handler.sort_vtws_by_increasing_conflict_degree()
        for request in self.request_handler:
            previous_conflict_degree = None
            for vtw in request[RequestIndices.vtw_list]:
                current_conflict_degree = self.request_handler._conflict_degree(vtw)
                if previous_conflict_degree is not None:
                    self.assertLessEqual(previous_conflict_degree, current_conflict_degree)
                previous_conflict_degree = current_conflict_degree

    def test_sort_vtws_by_increasing_rise_time(self):
        self.request_handler.sort_vtws_by_increasing_rise_time()
        for request in self.request_handler:
            previous_rise_time = None
            for vtw in request[RequestIndices.vtw_list]:
                if previous_rise_time is not None:
                    self.assertLessEqual(previous_rise_time, vtw.rise_time)
                previous_rise_time = vtw.rise_time
