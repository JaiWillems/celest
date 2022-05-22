

from celest.encounter.groundposition import GroundPosition
from unittest import TestCase
import numpy as np
import unittest


WGS84_MAJOR_AXIS = 6378.137
WGS84_MINOR_AXIS = 6356.752314245


class TestGroundPosition(TestCase):

    def setUp(self):

        self.latitude = 43.6532
        self.longitude = -79.3832
        self.elevation = 0.076

    def test_radius(self):
        """Test `GroundPosition._WGS84_radius`.

        Notes
        -----
        Validation against online calculator methodology.[1]_

        References
        ----------
        .. [1] Timur. Earth Radius by Latitude (WGS 84). 2018. url:
           https://planetcalc.com/7721/.
        """

        location = GroundPosition(self.latitude, self.longitude, self.elevation)

        cos_latitude = np.cos(np.radians(self.latitude))
        sin_latitude = np.sin(np.radians(self.latitude))

        numinator = (WGS84_MAJOR_AXIS ** 2 * cos_latitude) ** 2 + \
            (WGS84_MINOR_AXIS ** 2 * sin_latitude) ** 2
        denominator = (WGS84_MAJOR_AXIS * cos_latitude) ** 2 + \
            (WGS84_MINOR_AXIS * sin_latitude) ** 2
        radius = np.sqrt(numinator / denominator) + self.elevation

        self.assertAlmostEqual(radius, location.radius, delta=0.001)


if __name__ == "__main__":
    unittest.main()
