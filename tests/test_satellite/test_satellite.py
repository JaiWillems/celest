

from celest.encounter.groundposition import GroundPosition
from celest.satellite.satellite import Satellite
from unittest import TestCase
import numpy as np
import unittest


class TestSatellite(TestCase):

    def setUp(self):

        fname = "tests/test_data/coordinate_validation_long.txt"
        cols = (0, 11, 12, 13, 14, 15, 16)
        skiprows = 1
        max_rows = 5000
        data = np.loadtxt(fname=fname, usecols=cols, skiprows=skiprows,
                          max_rows=max_rows)

        self.times = data[:, 0]
        self.GCRS = data[:, 1:4]
        self.GCRS_vel = data[:, 4:7]

        self.offset = 2430000
        self.satellite = Satellite(self.GCRS, self.GCRS_vel, "gcrs",
                                   self.times, self.offset)

        self.julian = self.satellite.julian()
        self.ut1 = self.satellite.ut1()
        self.gmst = self.satellite.gmst()
        self.gast = self.satellite.gast()
        self.geo = self.satellite.geo()
        self.gcrs = self.satellite.gcrs()
        self.itrs = self.satellite.itrs()

    def test_attitude(self):

        location = GroundPosition(52.1579, -106.6702, 0.482)
        roll, pitch, yaw = self.satellite.attitude(location)

        self.assertIsNotNone(roll)
        self.assertIsNotNone(pitch)
        self.assertIsNotNone(yaw)

    def test_save(self):

        times = ("julian", "ut1", "gmst", "gast")
        positions = ("geo", "gcrs", "itrs")
        self.satellite.save(times, positions, path="test_data.csv", sep=",")

        data = np.loadtxt("test_data.csv", delimiter=",", skiprows=1)
        load_julian = data[:, 1]
        load_ut1 = data[:, 2]
        load_gmst = data[:, 3]
        load_gast = data[:, 4]
        load_lat = data[:, 5]
        load_lon = data[:, 6]
        load_alt = data[:, 7]
        load_gcrs = data[:, 8:14]
        load_itrs = data[:, 14:20]

        self.assertTrue(np.array_equal(load_julian, self.julian))
        self.assertTrue(np.array_equal(load_ut1, self.ut1))
        self.assertTrue(np.array_equal(load_gmst, self.gmst))
        self.assertTrue(np.array_equal(load_gast, self.gast))
        self.assertTrue(np.array_equal(load_lat, self.geo[0]))
        self.assertTrue(np.array_equal(load_lon, self.geo[1]))
        self.assertTrue(np.array_equal(load_alt, self.geo[2]))
        self.assertTrue(np.array_equal(load_gcrs[:, 0], self.gcrs[0]))
        self.assertTrue(np.array_equal(load_gcrs[:, 1], self.gcrs[1]))
        self.assertTrue(np.array_equal(load_gcrs[:, 2], self.gcrs[2]))
        self.assertTrue(np.array_equal(load_gcrs[:, 3], self.gcrs[3]))
        self.assertTrue(np.array_equal(load_gcrs[:, 4], self.gcrs[4]))
        self.assertTrue(np.array_equal(load_gcrs[:, 5], self.gcrs[5]))
        self.assertTrue(np.array_equal(load_itrs[:, 0], self.itrs[0]))
        self.assertTrue(np.array_equal(load_itrs[:, 1], self.itrs[1]))
        self.assertTrue(np.array_equal(load_itrs[:, 2], self.itrs[2]))
        self.assertTrue(np.array_equal(load_itrs[:, 3], self.itrs[3]))
        self.assertTrue(np.array_equal(load_itrs[:, 4], self.itrs[4]))
        self.assertTrue(np.array_equal(load_itrs[:, 5], self.itrs[5]))

        import os
        os.remove("test_data.csv")


if __name__ == "__main__":
    unittest.main()
