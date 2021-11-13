"""Testing module for the encounter mathematical basis."""


from celest.encounter.groundposition import GroundPosition
from celest.encounter._encounter_math_utils import _analytical_encounter_ind
from celest.satellite.coordinate import Coordinate
from celest.satellite.time import Time
from unittest import TestCase
import numpy as np
import unittest


class TestEncounterMathUtils(TestCase):

    def setUp(self):
        """Test fixure for test method execution."""

        fname = "tests/test_data/coordinate_validation_set.txt"
        data = np.loadtxt(fname=fname, delimiter="\t", skiprows=1)

        times = data[:, 0]
        itrs = data[:, 10:]

        self.timeData = Time(times, 2430000)
        self.coor = Coordinate(itrs, "itrs", self.timeData)
        self.length = data.shape[0]

    def test_anaytical_encounter_ind(self):
        """Test `analytical_encounter_ind`."""

        # Set up shared parameters.
        lat, lon = 43.6532, -79.3832
        GEO_data = self.coor._geo_to_itrs(np.array([[lat, lon, 0]]))
        gnd_itrs = np.repeat(GEO_data, self.length, 0)
        sat_itrs = self.coor.itrs()

        sat_itrs = self.coor.itrs()

        # Test imaging encounter.
        location = GroundPosition(lat, lon)

        altitude, _ = self.coor.horizontal(location)
        off_nadir = self.coor.off_nadir(location)

        ind = np.where((0 <= altitude) & (off_nadir < 30))[0]
        calc_ind = _analytical_encounter_ind(sat_itrs, gnd_itrs, 30, 1)

        self.assertTrue(np.array_equal(ind, calc_ind))

        # Test transmission encounter.
        location = GroundPosition(lat, lon)

        altitude, _ = self.coor.horizontal(location)

        ind = np.where(30 <= altitude)[0]
        calc_ind = _analytical_encounter_ind(sat_itrs, gnd_itrs, 30, 0)

        self.assertTrue(np.array_equal(ind, calc_ind))


if __name__ == "__main__":
    unittest.main()
