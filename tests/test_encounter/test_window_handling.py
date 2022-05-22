

from celest.encounter.groundposition import GroundPosition
from celest.encounter._window_handling import (
    VTW,
    OW,
    WindowHandler,
    VTWHandler,
    OWHandler
)
from unittest import TestCase
import numpy as np
import unittest
import os


class TestVTW(TestCase):

    def test_vtw_initialization(self):

        rise_time = 21282.4
        set_time = 21282.5
        self.window = VTW(rise_time, set_time, None, None, None)


class TestOW(TestCase):

    def test_ow_initialization(self):

        start_time = 21282.4
        duration = 30
        location = GroundPosition(43.6532, 79.3832, 0.076)
        deadline = 21284.3

        self.window = OW(start_time, duration, location, deadline, None, None, None)


class TestWindowHandler(TestCase):

    def setUp(self):

        self.windows = WindowHandler()

        self.W1 = VTW(1319, 1319.2, None, None, None)
        self.W2 = VTW(1322.3, 1322.31, None, None, None)
        self.W3 = VTW(1324, 1324.2, None, None, None)
        self.W4 = VTW(1324.3, 1325.3, None, None, None)
        self.W5 = VTW(1325, 1326.1, None, None, None)

        self.windows._add_window_to_window_handler(self.W1)
        self.windows._add_window_to_window_handler(self.W2)
        self.windows._add_window_to_window_handler(self.W3)
        self.windows._add_window_to_window_handler(self.W4)
        self.windows._add_window_to_window_handler(self.W5)

    def test_getitem(self):

        self.assertEqual(self.windows[1310], self.W1)
        self.assertEqual(self.windows[1319], self.W1)
        self.assertEqual(self.windows[1323], self.W2)
        self.assertEqual(self.windows[1324.4], self.W4)
        self.assertEqual(self.windows[1330], self.W5)

        test_windows = self.windows[1322.5, 1325]
        true_windows = np.array([self.W2, self.W5], dtype=object)
        self.assertTrue(np.array_equiv(test_windows, true_windows))

        self.assertEqual(self.windows[1323.5, 1324.1], self.W3)

        test_windows = self.windows[1322.5:1326]
        true_windows = np.array([self.W3, self.W4, self.W5], dtype=object)
        self.assertTrue(np.array_equiv(test_windows, true_windows))

        test_windows = self.windows[1322.3:1325]
        true_windows = np.array([self.W2, self.W3, self.W4, self.W5], dtype=object)
        self.assertTrue(np.array_equiv(test_windows, true_windows))

    def test_add_window(self):

        test_val = 0
        for w in self.windows:
            window_start = w.rise_time if isinstance(w, VTW) else w.start_time
            self.assertLessEqual(test_val, window_start)
            test_val = window_start

    def test_get_window(self):

        self.assertEqual(self.windows.get_window(1320), self.W1)
        self.assertEqual(self.windows.get_window(1322.5), self.W2)
        self.assertEqual(self.windows.get_window(1324.5), self.W4)
        self.assertEqual(self.windows.get_window(1325.5), self.W5)

    def test_get_windows_in_range(self):

        test_windows = self.windows.get_windows_in_range(1320, 1326)
        true_windows = np.array([self.W2, self.W3, self.W4, self.W5], dtype=object)
        self.assertTrue(np.array_equiv(test_windows, true_windows))

    def test_to_numpy(self):

        self.assertIsInstance(self.windows.to_numpy(), np.ndarray)


class TestVTWHandler(TestCase):

    def setUp(self) -> None:

        self.windows = VTWHandler()

        self.W1 = VTW(1319, 1319.2, None, None, None)
        self.W2 = VTW(1322.3, 1322.31, None, None, None)
        self.W3 = VTW(1324, 1324.2, None, None, None)
        self.W4 = VTW(1324.3, 1325.3, None, None, None)
        self.W5 = VTW(1325, 1326.1, None, None, None)

        self.windows._add_window(self.W1)
        self.windows._add_window(self.W2)
        self.windows._add_window(self.W3)
        self.windows._add_window(self.W4)
        self.windows._add_window(self.W5)

    def test_stats(self):

        self.assertIsNotNone(self.windows.stats())

    def test_save(self):

        self.windows.save('test_vtw.csv', ',')
        self.assertTrue(os.path.isfile('test_vtw.csv'))

        loaded_data = np.loadtxt('test_vtw.csv', delimiter=',', skiprows=1)
        for i, window in enumerate(self.windows):
            self.assertEqual(window.rise_time, loaded_data[i, 0])
            self.assertEqual(window.set_time, loaded_data[i, 1])
            self.assertEqual(window.duration, loaded_data[i, 2])

        os.remove('test_vtw.csv')


class TestOWHandler(TestCase):

    def setUp(self) -> None:

        self.windows = OWHandler()
        location = GroundPosition(43.6532, 79.3832, 0.076)

        self.W1 = OW(1319, 30, location, 1319.2, None, None, None)
        self.W2 = OW(1322.3, 20, location, 1322.31, None, None, None)
        self.W3 = OW(1324, 32, location, 1324.2, None, None, None)
        self.W4 = OW(1324.3, 28, location, 1325.3, None, None, None)
        self.W5 = OW(1325, 10, location, 1326.1, None, None, None)

        self.windows._add_window(self.W1)
        self.windows._add_window(self.W2)
        self.windows._add_window(self.W3)
        self.windows._add_window(self.W4)
        self.windows._add_window(self.W5)

    def test_save(self):

        self.windows.save('test_ow.csv', ',')
        self.assertTrue(os.path.isfile('test_ow.csv'))

        loaded_data = np.loadtxt('test_ow.csv', delimiter=',', skiprows=1)
        for i, window in enumerate(self.windows):

            self.assertEqual(window.location.latitude, loaded_data[i, 0])
            self.assertEqual(window.location.longitude, loaded_data[i, 1])
            self.assertEqual(window.start_time, loaded_data[i, 2])
            self.assertEqual(window.duration, loaded_data[i, 3])
            self.assertEqual(window.deadline, loaded_data[i, 4])

        os.remove('test_ow.csv')


if __name__ == "__main__":
    unittest.main()
