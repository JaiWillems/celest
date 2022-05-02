

from celest.satellite.satellite import Satellite
from unittest import TestCase
import numpy as np
import unittest


class TestSatellite(TestCase):

    def setUp(self):

        fname = "tests/test_data/coordinate_validation_set.txt"
        data = np.loadtxt(fname=fname, delimiter="\t", skiprows=1)

        self.times = data[:, 0]
        self.ITRS = data[:, 10:]

        self.offset = 2430000
        self.finch = Satellite(self.ITRS, "itrs", self.times, self.offset)

        self.julian = self.finch.julian()
        self.ut1 = self.finch.ut1()
        self.gmst = self.finch.gmst()
        self.gast = self.finch.gast()
        self.lat, self.lon, self.alt = self.finch.geo()
        self.gcrs_x, self.gcrs_y, self.gcrs_z = self.finch.gcrs()
        self.itrs_x, self.itrs_y, self.itrs_z = self.finch.itrs()

    def test_save(self):
        """Test `Satellite.save`."""

        times = ("julian", "ut1", "gmst", "gast")
        positions = ("geo", "gcrs", "itrs")
        self.finch.save(times, positions, path="test_data.csv", sep=",")

        data = np.loadtxt("test_data.csv", delimiter=",", skiprows=1)
        load_julian = data[:, 1]
        load_ut1 = data[:, 2]
        load_gmst = data[:, 3]
        load_gast = data[:, 4]
        load_lat = data[:, 5]
        load_lon = data[:, 6]
        load_alt = data[:, 7]
        load_gcrs = data[:, 8:11]
        load_itrs = data[:, 11:]

        self.assertTrue(np.array_equal(load_julian, self.julian))
        self.assertTrue(np.array_equal(load_ut1, self.ut1))
        self.assertTrue(np.array_equal(load_gmst, self.gmst))
        self.assertTrue(np.array_equal(load_gast, self.gast))
        self.assertTrue(np.array_equal(load_lat, self.lat))
        self.assertTrue(np.array_equal(load_lon, self.lon))
        self.assertTrue(np.array_equal(load_alt, self.alt))
        self.assertTrue(np.array_equal(load_gcrs[:, 0], self.gcrs_x))
        self.assertTrue(np.array_equal(load_gcrs[:, 1], self.gcrs_y))
        self.assertTrue(np.array_equal(load_gcrs[:, 2], self.gcrs_z))
        self.assertTrue(np.array_equal(load_itrs[:, 0], self.itrs_x))
        self.assertTrue(np.array_equal(load_itrs[:, 1], self.itrs_y))
        self.assertTrue(np.array_equal(load_itrs[:, 2], self.itrs_z))

        import os
        os.remove("test_data.csv")


if __name__ == "__main__":
    unittest.main()
