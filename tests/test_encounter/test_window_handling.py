

from celest.coordinates.frames.attitude import Attitude
from celest.coordinates.ground_location import GroundLocation
from celest.encounter.window_handling import (
    VisibleTimeWindow,
    ObservationWindow,
    WindowHandler
)
from celest import units as u
from unittest import TestCase
import numpy as np
import os


class TestVisibleTimeWindow(TestCase):

    def setUp(self):
        self.rise_time = 0
        self.set_time = 1
        self.attitude = Attitude(
            np.random.rand(5),
            np.random.rand(5),
            np.random.rand(5),
            np.random.rand(5),
            u.deg,
            GroundLocation(0, 0, 0, u.deg, u.km)
        )
        self.vtw = VisibleTimeWindow(
            self.rise_time,
            self.set_time,
            self.attitude
        )

    def test_initialization(self):
        self.assertIsInstance(self.vtw, VisibleTimeWindow)

    def test_str(self):
        self.assertEqual(f"Rise time: {self.rise_time} {u.jd2000}, Set time: "
                         f"{self.set_time} {u.jd2000}, Attitude: "
                         f"{self.attitude}", str(self.vtw))

    def test_repr(self):
        self.assertEqual(f"VisibleTimeWindow({self.rise_time}, "
                         f"{self.set_time}, {self.attitude})", repr(self.vtw))

    def test_rise_time_property(self):
        self.assertEqual(self.rise_time, self.vtw.rise_time.data)

    def test_set_time_property(self):
        self.assertEqual(self.set_time, self.vtw.set_time.data)

    def test_attitude_property(self):
        self.assertEqual(self.attitude, self.vtw.attitude)


class TestObservationWindow(TestCase):

    def setUp(self):
        self.start_time = 0
        self.duration = 1
        self.location = GroundLocation(0, 0, 0, u.deg, u.km)
        self.deadline = 2
        self.attitude = Attitude(
            np.random.rand(5),
            np.random.rand(5),
            np.random.rand(5),
            np.random.rand(5),
            u.deg,
            self.location
        )
        self.obw = ObservationWindow(
            self.start_time,
            self.duration,
            self.deadline,
            self.location,
            self.attitude
        )

    def test_initialization(self):
        self.assertIsInstance(self.obw, ObservationWindow)

    def test_str(self):
        self.assertEqual(f"Start time: {self.start_time} {u.jd2000}, "
                         f"Duration: {self.duration} {u.s}, Deadline: "
                         f"{self.deadline} {u.jd2000}, Location: "
                         f"{self.location}, Attitude: {self.attitude}",
                         str(self.obw))

    def test_repr(self):
        self.assertEqual(f"ObservationWindow({self.start_time}, "
                         f"{self.duration}, {self.deadline}, {self.location}, "
                         f"{self.attitude})", repr(self.obw))

    def test_start_time_property(self):
        self.assertEqual(self.start_time, self.obw.start_time.data)

    def test_duration_property(self):
        self.assertEqual(self.duration, self.obw.duration.data)

    def test_location_property(self):
        self.assertEqual(self.location, self.obw.location)

    def test_deadline_property(self):
        self.assertEqual(self.deadline, self.obw.deadline.data)

    def test_attitude_property(self):
        self.assertEqual(self.attitude, self.obw.attitude)


class TestWindowHandler(TestCase):

    def setUp(self):
        self.window_handler = WindowHandler()
        self.ground_location = GroundLocation(0, 0, 0, u.deg, u.km)
        self.vtw_attitude = Attitude(
            np.random.rand(5),
            np.random.rand(5),
            np.random.rand(5),
            np.random.rand(5),
            u.deg,
            self.ground_location
        )
        self.ow_attitude = Attitude(
            np.random.rand(1),
            np.random.rand(1),
            np.random.rand(1),
            np.random.rand(1),
            u.deg,
            self.ground_location
        )
        self.test_vtw = VisibleTimeWindow(0, 1, self.vtw_attitude)
        self.test_obw = ObservationWindow(0, 1, 2, self.ground_location, self.ow_attitude)

    def test_initialization(self):
        self.assertIsInstance(self.window_handler, WindowHandler)
        self.assertListEqual(self.window_handler._window_data, [])

    def test_add_window_raises_error_for_non_window(self):
        self.assertRaises(TypeError, self.window_handler.add_window,
                          "not a window")

    def test_add_window_alters_window_data_with_vtw(self):
        self.window_handler.add_window(self.test_vtw)
        self.assertListEqual(self.window_handler._window_data, [self.test_vtw])

    def test_add_window_alters_window_data_with_ow(self):
        self.window_handler.add_window(self.test_obw)
        self.assertListEqual(self.window_handler._window_data, [self.test_obw])

    def test_add_window_enforces_same_window_type(self):
        self.window_handler.add_window(self.test_vtw)
        self.assertRaises(TypeError, self.window_handler.add_window,
                          self.test_obw)

    def test_str(self):
        self.window_handler.add_window(self.test_vtw)
        self.window_handler.add_window(self.test_vtw)
        self.assertEqual(str([self.test_vtw, self.test_vtw]),
                         str(self.window_handler))

    def test_repr(self):
        self.window_handler.add_window(self.test_vtw)
        self.window_handler.add_window(self.test_vtw)
        self.assertEqual(f"WindowHandler({repr(self.test_vtw)}, {repr(self.test_vtw)})",
                         repr(self.window_handler))

    def test_len_returns_correct_length(self):
        self.assertEqual(len(self.window_handler), 0)
        self.window_handler.add_window(self.test_vtw)
        self.assertEqual(len(self.window_handler), 1)
        self.window_handler.add_window(self.test_vtw)
        self.assertEqual(len(self.window_handler), 2)

    def test_iteration(self):
        self.window_handler.add_window(self.test_vtw)
        self.window_handler.add_window(self.test_vtw)
        self.assertListEqual([window for window in self.window_handler],
                             self.window_handler._window_data)

    def test_indexing_with_exact_start_time(self):
        pass

    def test_indexing_with_inexact_start_time(self):
        pass

    def test_indexing_with_range(self):
        pass

    def test_save_text_file_raises_error_for_no_windows(self):
        self.assertRaises(Exception, self.window_handler.save_text_file)

    def test_save_text_file_for_vtws(self):
        self.window_handler.add_window(self.test_vtw)
        self.window_handler.add_window(self.test_vtw)

        file_name = "vtw_test_file"
        self.window_handler.save_text_file(file_name)
        self.assertTrue(os.path.exists(file_name + ".txt"))
        if os.path.exists(file_name + ".txt"):
            os.remove(file_name + ".txt")

    def test_save_text_file_for_ows(self):
        self.window_handler.add_window(self.test_obw)
        self.window_handler.add_window(self.test_obw)

        file_name = "ow_test_file"
        self.window_handler.save_text_file(file_name)
        self.assertTrue(os.path.exists(file_name + ".txt"))
        if os.path.exists(file_name + ".txt"):
            os.remove(file_name + ".txt")
