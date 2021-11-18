"""Testing module for the window data structures."""


from celest.encounter.groundposition import GroundPosition
from celest.encounter._window_handling import Window, Windows
from unittest import TestCase
import numpy as np
import unittest


class TestWindow(TestCase):

    def setUp(self):
        """Test fixure for test method execution."""

        location = GroundPosition(43.6532, 79.3832)
        self.window = Window("", location, 21282.4, 21282.5, "image", 30, 1, 7)

    def test_str(self):
        """Test `Window.__str__`."""

        print(self.window)

    def test_copy(self):
        """Test `Window.copy`."""

        wind_copy = self.window.copy()

        self.assertTrue(self.window.satellite == wind_copy.satellite)
        self.assertTrue(self.window.coor == wind_copy.coor)
        self.assertTrue(self.window.type == wind_copy.type)
        self.assertTrue(self.window.start == wind_copy.start)
        self.assertTrue(self.window.end == wind_copy.end)


class TestWindows(TestCase):

    def setUp(self):
        """Test fixure for test method execution."""

        location = GroundPosition(43.6532, 79.3832)
        self.windows = Windows()

        self.wind_one = Window("", location, 1319, 1319.2, "image", 30, 1, 0)
        self.wind_two = Window("", location, 1322.3, 1322.31, "data link", 10, -1, 30)
        self.wind_three = Window("", location, 1324, 1324.2, "image", 30, 1, 0)
        self.wind_four = Window("", location, 1324.3, 1325.3, "image", 30, 1, 0)
        self.wind_five = Window("", location, 1325, 1326.1, "image", 30, 1, 0)

        self.windows._add_window(self.wind_one)
        self.windows._add_window(self.wind_two)
        self.windows._add_window(self.wind_three)
        self.windows._add_window(self.wind_four)
        self.windows._add_window(self.wind_five)
    
    def test_getitem(self):
        """Test `Window.__getitem__`."""

        # Test int/float indexing.
        self.assertEqual(self.windows[1310], self.wind_one)
        self.assertEqual(self.windows[1319], self.wind_one)
        self.assertEqual(self.windows[1323], self.wind_two)
        self.assertEqual(self.windows[1324.4], self.wind_four)
        self.assertEqual(self.windows[1330], self.wind_five)

        # Test tuple indexing with unique mapping.
        arr = [self.wind_two, self.wind_five]
        val_wind = np.array(arr, dtype=object)
        self.assertTrue(np.array_equiv(self.windows[1322.5, 1325], val_wind))

        # Test tuple indexing without unique mapping.
        self.assertEqual(self.windows[1323.5, 1324.1], self.wind_three)

        # Test slice indexing without boundaries cases.
        arr = [self.wind_three, self.wind_four, self.wind_five]
        val_wind = np.array(arr, dtype=object)
        self.assertTrue(np.array_equiv(self.windows[1322.5:1326], val_wind))

        # Test slice indexing with boundary cases.
        arr = [self.wind_two, self.wind_three, self.wind_four, self.wind_five]
        val_wind = np.array(arr, dtype=object)
        self.assertTrue(np.array_equiv(self.windows[1322.3:1325], val_wind))

    def test_add_window(self):
        """Test `Window._add_window`."""

        test_val = 0
        for window in self.windows:
            start = window.start
            self.assertLessEqual(test_val, start)
            test_val = start

    def test_stats(self):
        """Test `Window.stats`."""

        df = self.windows.stats()
        self.assertIsNotNone(df)

    def test_to_numpy(self):
        """Test `Window.to_numpy`."""
        
        np_arr = self.windows.to_numpy()
        self.assertIsInstance(np_arr, np.ndarray)

    def test_save_encounters(self):
        pass


if __name__ == "__main__":
    unittest.main()