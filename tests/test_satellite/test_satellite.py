

import numpy as np
import unittest
from unittest import TestCase
from celest.encounter import Encounter, GroundPosition
from celest.satellite import Coordinate, Satellite, Time


class TestSatellite(TestCase):

    def setUp(self):
        """Test fixure for test method execution."""

        fname = "tests/test_data/coordinate_validation_set.txt"
        data = np.loadtxt(fname=fname, delimiter="\t", skiprows=1)

        self.times = data[:, 0]
        self.ECEF = data[:, 10:]

        self.timeData = Time(self.times, 2430000)
        self.coor = Coordinate(self.ECEF, "ECEF", self.timeData)

        self.finch = Satellite(self.coor)

    def test_generate_pointing_profiles(self):
        """Test `Satellite.generate_pointing_profiles`."""

        groundPos = GroundPosition("Toronto", (43.6532, -79.3832), "image", 30)

        encounter = Encounter(self.finch)
        enc_ind = encounter._window_encounter_ind(groundPos=groundPos)

        rotations = self.finch.generate_pointing_profiles(groundPos, enc_ind, 10)

        self.assertIsNotNone(rotations)

    def test_save_data(self):
        """Test `Satellite.save_data`."""

        self.finch.save_data("test_data.csv", ",", ["GEO", "ECI", "ECEF"])

        data = np.loadtxt("test_data.csv", delimiter=",", skiprows=1)
        load_julian = data[:, 1]
        load_GEO = data[:, 2:5]
        load_ECI = data[:, 5:8]
        load_ECEF = data[:, 8:]

        julian = self.finch.times.julian()
        GEO = self.finch.position.GEO()
        ECI = self.finch.position.ECI()
        ECEF = self.finch.position.ECEF()

        self.assertTrue(np.array_equal(load_julian, julian))
        self.assertTrue(np.array_equal(load_GEO, GEO))
        self.assertTrue(np.array_equal(load_ECI, ECI))
        self.assertTrue(np.array_equal(load_ECEF, ECEF))

        import os
        os.remove("test_data.csv")


if __name__ == "__main__":
    unittest.main()
