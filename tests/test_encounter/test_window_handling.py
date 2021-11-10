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

        window_copy = self.window.copy()

        self.assertTrue(self.window.satellite == window_copy.satellite)
        self.assertTrue(self.window.coor == window_copy.coor)
        self.assertTrue(self.window.type == window_copy.type)
        self.assertTrue(self.window.start == window_copy.start)
        self.assertTrue(self.window.end == window_copy.end)


class TestWindows(TestCase):

    def setUp(self):
        """Test fixure for test method execution."""


        location = GroundPosition(43.6532, 79.3832)
        self.windows = Windows()

        self.window_one = Window("", location, 1319, 1319.2, "image", 30, 1, 0)
        self.window_two = Window("", location, 1322.3, 1322.31, "data link", 10, -1, 30)
        self.window_three = Window("", location, 1320, 1320.2, "image", 30, 1, 0)
        self.window_three = Window("", location, 1320, 1320.3, "image", 30, 1, 0)
        self.window_three = Window("", location, 1320, 1320.1, "image", 30, 1, 0)

        self.windows._add_window(self.window_one)
        self.windows._add_window(self.window_two)
        self.windows._add_window(self.window_three)

    def test_add_window(self):
        """Test `Window._add_window`."""

        test_val = 0
        for window in self.windows:
            start = window.start
            self.assertLessEqual(test_val, start)
            test_val = start

    def test_encounter_stats(self):
        """Test `Window.encounter_stats`."""

        df = self.windows.encounter_stats()
        self.assertIsNotNone(df)

    def test_to_numpy(self):
        """Test `Window.to_numpy`."""
        
        np_arr = self.windows.to_numpy()
        self.assertIsInstance(np_arr, np.ndarray)

    def test_save_encounters(self):
        pass


if __name__ == "__main__":
    unittest.main()