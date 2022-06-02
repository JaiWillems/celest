

from celest.schedule.request_handler import RequestHandler, RequestIndices
from celest.schedule.scheduling_utils import (
    cost,
    initialize_solution,
    is_complete
)
from tests.test_schedule.test_request_handler_utils import initialize_request_list
from unittest import TestCase


class TestSchedulingUtils(TestCase):

    def setUp(self) -> None:

        self.request_handler = RequestHandler()
        for request in initialize_request_list():
            self.request_handler.add_request(*request)

    def test_initialize_solution(self):

        initialize_solution(self.request_handler)
        self.assertEqual(self.request_handler.number_of_scheduled_requests, 7)

    def test_is_complete_none_scheduled(self):

        self.assertFalse(is_complete(self.request_handler))

    def test_is_complete_some_scheduled(self):

        self.request_handler.schedule_request(0, 0, 0, 0)
        self.request_handler.schedule_request(1, 0, 0, 0)
        self.assertFalse(is_complete(self.request_handler))

    def test_is_complete_all_scheduled(self):

        self.request_handler.schedule_request(0, 0, 0, 0)
        self.request_handler.schedule_request(1, 0, 0, 0)
        self.request_handler.schedule_request(2, 0, 0, 0)
        self.request_handler.schedule_request(3, 0, 0, 0)
        self.request_handler.schedule_request(4, 0, 0, 0)
        self.request_handler.schedule_request(5, 0, 0, 0)
        self.request_handler.schedule_request(6, 0, 0, 0)
        self.request_handler.schedule_request(7, 0, 0, 0)
        self.request_handler.schedule_request(8, 0, 0, 0)
        self.request_handler.schedule_request(9, 0, 0, 0)
        self.assertTrue(is_complete(self.request_handler))

    def test_cost_none_scheduled(self):

        self.assertEqual(cost(self.request_handler), 0)

    def test_cost_some_scheduled(self):

        self.request_handler.schedule_request(0, 0, 0, 0)
        self.request_handler.schedule_request(1, 0, 0, 0)
        actual_cost = - (self.request_handler[0][RequestIndices.priority] +
                         self.request_handler[1][RequestIndices.priority])
        self.assertEqual(cost(self.request_handler), actual_cost)
