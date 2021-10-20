

import numpy as np
import unittest
from unittest import TestCase
from celest.encounter.groundposition import GroundPosition
from celest.encounter._window_handling import Window, Windows
from celest.satellite.coordinate import Coordinate
from celest.satellite.satellite import Satellite
from celest.satellite.time import Time


class TestSatellite(TestCase):

    def setUp(self):
        """Test fixure for test method execution."""

        fname = "tests/test_data/coordinate_validation_set.txt"
        data = np.loadtxt(fname=fname, delimiter="\t", skiprows=1)

        self.times = data[:, 0]
        self.ECEF = data[:, 10:]

        self.timeData = Time(self.times, 2430000)
        self.coor = Coordinate(self.ECEF, "ecef", self.timeData)

        self.finch = Satellite(self.coor)

    def test_interpolate(self):
        """Test `Satellite.interpolate`."""

        location = GroundPosition(0, 0)

        window = Window(None, location, self.times[100] + 2430000, self.times[150] + 2430000, None, None, None, None)
        window_list = Windows()
        window_list._add_window(window)

        len_i = len(self.finch)

        self.finch.interpolate(windows=window_list, factor=5)

        len_f = len(self.finch)

        self.assertGreater(len_f, len_i)

    def test_save_data(self):
        """Test `Satellite.save_data`."""

        times = ("julian", "ut1", "gmst", "gast")
        positions = ("geo", "eci", "ecef")
        self.finch.save_data("test_data.csv", ",", times, positions)

        data = np.loadtxt("test_data.csv", delimiter=",", skiprows=1)
        load_julian = data[:, 1]
        load_ut1 = data[:, 2]
        load_gmst = data[:, 3]
        load_gast = data[:, 4]
        load_geo = data[:, 5:8]
        load_eci = data[:, 8:11]
        load_ecef = data[:, 11:]

        julian = self.finch.time.julian()
        ut1 = self.finch.time.UT1()
        gmst = self.finch.time.GMST()
        gast = self.finch.time.GAST()
        geo = self.finch.position.GEO()
        eci = self.finch.position.ECI()
        ecef = self.finch.position.ECEF()

        self.assertTrue(np.array_equal(load_julian, julian))
        self.assertTrue(np.array_equal(load_ut1, ut1))
        self.assertTrue(np.array_equal(load_gmst, gmst))
        self.assertTrue(np.array_equal(load_gast, gast))
        self.assertTrue(np.array_equal(load_geo, geo))
        self.assertTrue(np.array_equal(load_eci, eci))
        self.assertTrue(np.array_equal(load_ecef, ecef))

        import os
        os.remove("test_data.csv")


if __name__ == "__main__":
    unittest.main()
