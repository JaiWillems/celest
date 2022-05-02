

from celest.encounter.groundposition import GroundPosition
from celest.encounter._window_handling import VTW, OW, WindowHandler, VTWHandler, OWHandler
from unittest import TestCase
import numpy as np
import unittest
import os


class TestVTW(TestCase):

    def setUp(self):

        rise_time = 21282.4
        set_time = 21282.5

        self.window = VTW(rise_time, set_time)

    def test_copy(self):

        wind_copy = self.window.copy()

        self.assertTrue(self.window.rise_time == wind_copy.rise_time)
        self.assertTrue(self.window.set_time == wind_copy.set_time)
        self.assertTrue(self.window.duration == wind_copy.duration)


class TestOW(TestCase):

    def setUp(self):

        start_time = 21282.4
        duration = 30
        location = GroundPosition(43.6532, 79.3832, 0.076)
        quality = 9
        deadline = 21284.3

        self.window = OW(start_time, duration, location, quality, deadline)

    def test_copy(self):

        wind_copy = self.window.copy()

        self.assertTrue(self.window.start_time == wind_copy.start_time)
        self.assertTrue(self.window.duration == wind_copy.duration)
        self.assertTrue(self.window.location.lat == wind_copy.location.lat)
        self.assertTrue(self.window.location.lon == wind_copy.location.lon)
        self.assertTrue(self.window.quality == wind_copy.quality)
        self.assertTrue(self.window.deadline == wind_copy.deadline)


class TestWindowHandler(TestCase):

    def setUp(self):

        self.windows = WindowHandler()

        self.W1 = VTW(1319, 1319.2)
        self.W2 = VTW(1322.3, 1322.31)
        self.W3 = VTW(1324, 1324.2)
        self.W4 = VTW(1324.3, 1325.3)
        self.W5 = VTW(1325, 1326.1)

        self.windows._add_window_base(self.W1)
        self.windows._add_window_base(self.W2)
        self.windows._add_window_base(self.W3)
        self.windows._add_window_base(self.W4)
        self.windows._add_window_base(self.W5)

    def test_getitem(self):

        # Test int/float indexing.
        self.assertEqual(self.windows[1310], self.W1)
        self.assertEqual(self.windows[1319], self.W1)
        self.assertEqual(self.windows[1323], self.W2)
        self.assertEqual(self.windows[1324.4], self.W4)
        self.assertEqual(self.windows[1330], self.W5)

        # Test tuple indexing with unique mapping.
        arr = [self.W2, self.W5]
        val_wind = np.array(arr, dtype=object)
        self.assertTrue(np.array_equiv(self.windows[1322.5, 1325], val_wind))

        # Test tuple indexing without unique mapping.
        self.assertEqual(self.windows[1323.5, 1324.1], self.W3)

        # Test slice indexing without boundaries cases.
        arr = [self.W3, self.W4, self.W5]
        val_wind = np.array(arr, dtype=object)
        self.assertTrue(np.array_equiv(self.windows[1322.5:1326], val_wind))

        # Test slice indexing with boundary cases.
        arr = [self.W2, self.W3, self.W4, self.W5]
        val_wind = np.array(arr, dtype=object)
        self.assertTrue(np.array_equiv(self.windows[1322.3:1325], val_wind))

    def test_add_window(self):

        test_val = 0
        for window in self.windows:
            start = window.rise_time if isinstance(window, VTW) else window.start_time
            self.assertLessEqual(test_val, start)
            test_val = start

    def test_get_window(self):

        self.assertEqual(self.windows.get_window(1320), self.W1)
        self.assertEqual(self.windows.get_window(1322.5), self.W2)
        self.assertEqual(self.windows.get_window(1324.5), self.W4)
        self.assertEqual(self.windows.get_window(1325.5), self.W5)

    def test_get_windows_in_range(self):

        arr = [self.W2, self.W3, self.W4, self.W5]
        val_wind = np.array(arr, dtype=object)
        self.assertTrue(np.array_equiv(self.windows.get_windows_in_range(1320, 1326), val_wind))

    def test_to_numpy(self):

        np_arr = self.windows.to_numpy()
        self.assertIsInstance(np_arr, np.ndarray)


class TestVTWHandler(TestCase):

    def setUp(self) -> None:

        self.windows = VTWHandler()

        self.W1 = VTW(1319, 1319.2)
        self.W2 = VTW(1322.3, 1322.31)
        self.W3 = VTW(1324, 1324.2)
        self.W4 = VTW(1324.3, 1325.3)
        self.W5 = VTW(1325, 1326.1)

        self.windows._add_window(self.W1)
        self.windows._add_window(self.W2)
        self.windows._add_window(self.W3)
        self.windows._add_window(self.W4)
        self.windows._add_window(self.W5)

    def test_stats(self):

        stats = self.windows.stats()
        self.assertIsNotNone(stats)

    def test_save(self):

        self.windows.save('test_vtw.csv', ',')
        self.assertTrue(os.path.isfile('test_vtw.csv'))

        data = np.loadtxt('test_vtw.csv', delimiter=',', skiprows=1)
        for i, window in enumerate(self.windows):

            self.assertEqual(window.rise_time, data[i, 0])
            self.assertEqual(window.set_time, data[i, 1])
            self.assertEqual(window.duration, data[i, 2])

        os.remove('test_vtw.csv')


class TestOWHandler(TestCase):

    def setUp(self) -> None:

        self.windows = OWHandler()

        location = GroundPosition(43.6532, 79.3832, 0.076)

        self.W1 = OW(1319, 30, location, 8, 1319.2)
        self.W2 = OW(1322.3, 20, location, 9, 1322.31)
        self.W3 = OW(1324, 32, location, 4, 1324.2)
        self.W4 = OW(1324.3, 28, location, 7, 1325.3)
        self.W5 = OW(1325, 10, location, 10, 1326.1)

        self.windows._add_window(self.W1)
        self.windows._add_window(self.W2)
        self.windows._add_window(self.W3)
        self.windows._add_window(self.W4)
        self.windows._add_window(self.W5)

    def test_save(self):

        self.windows.save('test_ow.csv', ',')
        self.assertTrue(os.path.isfile('test_ow.csv'))

        data = np.loadtxt('test_ow.csv', delimiter=',', skiprows=1)
        for i, window in enumerate(self.windows):

            self.assertEqual(window.location.lat, data[i, 0])
            self.assertEqual(window.location.lon, data[i, 1])
            self.assertEqual(window.start_time, data[i, 2])
            self.assertEqual(window.duration, data[i, 3])
            self.assertEqual(window.quality, data[i, 4])
            self.assertEqual(window.deadline, data[i, 5])

        os.remove('test_ow.csv')


if __name__ == "__main__":
    unittest.main()
