
from celest.schedule.removal_operators import (
    conflict_removal,
    opportunity_removal,
    priority_removal,
    random_removal
)
from celest.schedule.request_handler import RequestHandler
from tests.test_schedule.request_handler_testing_utils import initialize_request_list
from unittest import TestCase


class TestRemovalOperators(TestCase):

    def setUp(self):
        self.request_handler = RequestHandler()

        for request in initialize_request_list():
            self.request_handler.add_request(*request)

    def test_random_removal_with_no_items(self):
        request_handler = random_removal(self.request_handler, 1)
        self.assertEqual(request_handler.number_of_scheduled_requests, 0)

    def test_random_removal_one_item(self):
        self._schedule_all_requests()
        request_handler = random_removal(self.request_handler, 1)
        self.assertEqual(request_handler.number_of_scheduled_requests, 9)

    def _schedule_all_requests(self):
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

    def test_random_removal_two_items(self):
        self._schedule_all_requests()
        request_handler = random_removal(self.request_handler, 2)
        self.assertEqual(request_handler.number_of_scheduled_requests, 8)

    def test_priority_removal_with_no_items(self):
        request_handler = priority_removal(self.request_handler, 1)
        self.assertEqual(request_handler.number_of_scheduled_requests, 0)

    def test_priority_removal_one_item(self):
        self._schedule_all_requests()
        request_handler = priority_removal(self.request_handler, 1)
        self.assertEqual(request_handler.number_of_scheduled_requests, 9)

    def test_priority_removal_two_items(self):
        self._schedule_all_requests()
        request_handler = priority_removal(self.request_handler, 2)
        self.assertEqual(request_handler.number_of_scheduled_requests, 8)

    def test_opportunity_removal_with_no_items(self):
        request_handler = opportunity_removal(self.request_handler, 1)
        self.assertEqual(request_handler.number_of_scheduled_requests, 0)

    def test_opportunity_removal_one_item(self):
        self._schedule_all_requests()
        request_handler = opportunity_removal(self.request_handler, 1)
        self.assertEqual(request_handler.number_of_scheduled_requests, 9)

    def test_opportunity_removal_two_items(self):
        self._schedule_all_requests()
        request_handler = opportunity_removal(self.request_handler, 2)
        self.assertEqual(request_handler.number_of_scheduled_requests, 8)

    def test_conflict_removal_with_no_items(self):
        request_handler = conflict_removal(self.request_handler, 1)
        self.assertEqual(request_handler.number_of_scheduled_requests, 0)

    def test_conflict_removal_one_item(self):
        self._schedule_all_requests()
        request_handler = conflict_removal(self.request_handler, 1)
        self.assertEqual(request_handler.number_of_scheduled_requests, 9)

    def test_conflict_removal_two_items(self):
        self._schedule_all_requests()
        request_handler = conflict_removal(self.request_handler, 2)
        self.assertEqual(request_handler.number_of_scheduled_requests, 8)
