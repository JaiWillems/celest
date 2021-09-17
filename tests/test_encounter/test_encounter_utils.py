from celest.encounter import GroundPosition, analytical_encounter_ind
from celest.satellite import Coordinate, Time
from unittest import TestCase
import numpy as np
import unittest


class TestEncounterUtils(TestCase):

    def setUp(self):
        """Test fixure for test method execution."""

        fname = "tests/test_data/coordinate_validation_set.txt"
        data = np.loadtxt(fname=fname, delimiter="\t", skiprows=1)

        times = data[:, 0]
        ECEF = data[:, 10:]

        self.timeData = Time(times, 2430000)
        self.coor = Coordinate(ECEF, "ECEF", self.timeData)
        self.length = data.shape[0]

    def test_anaytical_encounter_ind(self):
        """Test `analytical_encounter_ind`."""

        # Set up shared parameters.
        coor = (43.6532, -79.3832)
        GEO_data = self.coor._GEO_to_ECEF(np.array([[coor[0], coor[1], 0]]))
        gndECEF = np.repeat(GEO_data, self.length, 0)
        satECEF = self.coor.ECEF()

        satECEF = self.coor.ECEF()

        # Test imaging encounter.
        groundPos = GroundPosition("", coor, "image", 30)

        altitude, _ = self.coor.horizontal(groundPos=groundPos)
        off_nadir = self.coor.off_nadir(groundPos=groundPos)

        ind = np.where((0 <= altitude) & (off_nadir < 30))[0]
        calc_ind = analytical_encounter_ind(satECEF, gndECEF, 30, 1)

        self.assertTrue(np.array_equal(ind, calc_ind))

        # Test transmission encounter.
        groundPos = GroundPosition("", coor, "data_link", 30)

        altitude, _ = self.coor.horizontal(groundPos=groundPos)

        ind = np.where(30 <= altitude)[0]
        calc_ind = analytical_encounter_ind(satECEF, gndECEF, 30, 0)

        self.assertTrue(np.array_equal(ind, calc_ind))


if __name__ == "__main__":
    unittest.main()
