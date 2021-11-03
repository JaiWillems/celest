"""Testing module for the `Satellite` class."""


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
        self.ITRS = data[:, 10:]

        self.timeData = Time(self.times, 2430000)
        self.coor = Coordinate(self.ITRS, "itrs", self.timeData)

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
        positions = ("geo", "gcrs", "itrs")
        self.finch.save_data("test_data.csv", ",", times, positions)

        data = np.loadtxt("test_data.csv", delimiter=",", skiprows=1)
        load_julian = data[:, 1]
        load_ut1 = data[:, 2]
        load_gmst = data[:, 3]
        load_gast = data[:, 4]
        load_geo = data[:, 5:8]
        load_gcrs = data[:, 8:11]
        load_itrs = data[:, 11:]

        julian = self.finch.time.julian()
        ut1 = self.finch.time.ut1()
        gmst = self.finch.time.gmst()
        gast = self.finch.time.gast()
        geo = self.finch.position.geo()
        gcrs = self.finch.position.gcrs()
        itrs = self.finch.position.itrs()

        self.assertTrue(np.array_equal(load_julian, julian))
        self.assertTrue(np.array_equal(load_ut1, ut1))
        self.assertTrue(np.array_equal(load_gmst, gmst))
        self.assertTrue(np.array_equal(load_gast, gast))
        self.assertTrue(np.array_equal(load_geo, geo))
        self.assertTrue(np.array_equal(load_gcrs, gcrs))
        self.assertTrue(np.array_equal(load_itrs, itrs))

        import os
        os.remove("test_data.csv")


if __name__ == "__main__":
    unittest.main()
