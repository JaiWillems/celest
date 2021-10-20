

import numpy as np
import unittest
from unittest import TestCase
from celest.encounter.groundposition import GroundPosition
from celest.encounter._window_handling import Window, Windows
from celest.encounter.encounter import generate
from celest.encounter._window_utils import _window_encounter_ind
from celest.satellite.coordinate import Coordinate
from celest.satellite.satellite import Satellite
from celest.satellite.time import Time


class TestEncounter(TestCase):

    def setUp(self):
        """Test fixure for test method execution."""

        fname = "tests/test_data/coordinate_validation_set.txt"
        data = np.loadtxt(fname=fname, delimiter="\t", skiprows=1)
        times, ECEF = data[:, 0], data[:, 10:]

        timeData = Time(times, 2430000)
        coor = Coordinate(ECEF, "ecef", timeData)

        self.finch = Satellite(coor)

    def test_windows(self):
        """Test `Encounter.windows.generate.`"""

        location_1 = GroundPosition(43.6532, -79.3832)
        location_2 = GroundPosition(52.1579, -106.6702)

        # Test case 1.
        windows_1 = generate(self.finch, location_1, "image", 30, 0)
        indices = _window_encounter_ind(self.finch, location_1, 30, 1, 0, 1)
        indices = np.split(indices, np.where(np.diff(indices) != 1)[0] + 1)

        val_windows = Windows()
        time = self.finch.time.julian()
        for reg in indices:
                
            start = time[reg[0]]
            end = time[reg[-1]]
                
            window = Window(self.finch, location_1, start, end, "image", 30, 1, 0)
            val_windows._add_window(window)
        
        self.assertEqual(len(val_windows), len(windows_1))

        windows_val, windows_calc = val_windows.windows, windows_1.windows
        for val, calc in zip(windows_val, windows_calc):
            self.assertEqual(val.start, calc.start)
            self.assertEqual(val.end, calc.end)

        # Test case 2.
        windows_2 = generate(self.finch, location_1, "data_link", 10, 0)
        indices = _window_encounter_ind(self.finch, location_1, 10, 0, 30, 0)
        indices = np.split(indices, np.where(np.diff(indices) != 1)[0] + 1)

        val_windows = Windows()
        time = self.finch.time.julian()
        for reg in indices:
                
            start = time[reg[0]]
            end = time[reg[-1]]
                
            window = Window(self.finch, location_1, start, end, "data_link", 10, 0, 30)
            val_windows._add_window(window)
        
        self.assertEqual(len(val_windows), len(windows_2))

        windows_val, windows_calc = val_windows.windows, windows_2.windows
        for val, calc in zip(windows_val, windows_calc):
            self.assertEqual(val.start, calc.start)
            self.assertEqual(val.end, calc.end)

        windows_3 = generate(self.finch, location_2, "image", 30, 0)
        indices = _window_encounter_ind(self.finch, location_2, 30, 1, 0, 1)
        indices = np.split(indices, np.where(np.diff(indices) != 1)[0] + 1)

        val_windows = Windows()
        time = self.finch.time.julian()
        for reg in indices:
                
            start = time[reg[0]]
            end = time[reg[-1]]
                
            window = Window(self.finch, location_2, start, end, "image", 30, 1, 0)
            val_windows._add_window(window)
        
        self.assertEqual(len(val_windows), len(windows_3))

        windows_val, windows_calc = val_windows.windows, windows_3.windows
        for val, calc in zip(windows_val, windows_calc):
            self.assertEqual(val.start, calc.start)
            self.assertEqual(val.end, calc.end)
        

        windows_4 = generate(self.finch, location_2, "data_link", 10)
        indices = _window_encounter_ind(self.finch, location_2, 10, 0, 30, 0)
        indices = np.split(indices, np.where(np.diff(indices) != 1)[0] + 1)

        val_windows = Windows()
        time = self.finch.time.julian()
        for reg in indices:
                
            start = time[reg[0]]
            end = time[reg[-1]]
                
            window = Window(self.finch, location_1, start, end, "data_link", 30, 1, 0)
            val_windows._add_window(window)
        
        self.assertEqual(len(val_windows), len(windows_4))

        windows_val, windows_calc = val_windows.windows, windows_4.windows
        for val, calc in zip(windows_val, windows_calc):
            self.assertEqual(val.start, calc.start)
            self.assertEqual(val.end, calc.end)


if __name__ == "__main__":
    unittest.main()