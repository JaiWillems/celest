

from celest.satellite._angle_representations import sexagesimal, _ISO6709_representation
from unittest import TestCase
import numpy as np


class TestAngleRepresentations(TestCase):

    def test_sexagesimal(self):
        """Test `Coordinate.sexagesimal`."""

        ang = np.array([43.6532, -79.3832, -33.2833, 149.1000])
        true_sexagesimal = ["+43\u00B039\u203211.52\u2033",
                            "-79\u00B022\u203259.52\u2033",
                            "-33\u00B016\u203259.88\u2033",
                            "+149\u00B006\u203200.00\u2033"]

        calc_ang = sexagesimal(ang)

        for i in range(calc_ang.shape[0]):
            with self.subTest(i=i):
                self.assertTrue(true_sexagesimal[i] == calc_ang[i])

    def test_ISO6709_representation(self):
        """Test `Coordinate._ISO6709_representation`."""

        latitude = np.array([43.6532, -33.2833])
        longitude = np.array([-79.3832, 149.1000])
        altitude = np.array([430.23, 532.98])

        true_sexagesimal = ["43\u00B039\u203211.52\u2033N 79\u00B022\u203259.52\u2033W 430.23km",
                            "33\u00B016\u203259.88\u2033S 149\u00B006\u203200.00\u2033E 532.98km"]

        calc_ang = _ISO6709_representation(latitude, longitude, altitude)

        for i in range(calc_ang.shape[0]):
            with self.subTest(i=i):
                self.assertTrue(true_sexagesimal[i] == calc_ang[i])
