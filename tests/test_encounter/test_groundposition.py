

from celest.encounter.groundposition import GroundPosition
from unittest import TestCase
import numpy as np
import unittest


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

        ground_pos = GroundPosition(self.latitude, self.longitude, self.elevation)

        a = 6378.1370
        b = 6356.7523142
        latitude = np.radians(ground_pos.latitude)
        clat, slat = np.cos(latitude), np.sin(latitude)
        num = (a ** 2 * clat) ** 2 + (b ** 2 * slat) ** 2
        denom = (a * clat) ** 2 + (b * slat) ** 2
        radius = np.sqrt(num / denom) + self.elevation

        self.assertAlmostEqual(radius, ground_pos.radius, delta=0.001)


if __name__ == "__main__":
    unittest.main()
