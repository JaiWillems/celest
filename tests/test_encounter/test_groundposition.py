"""Test module for the `GroundPosition` class."""


from celest.encounter.groundposition import GroundPosition
from unittest import TestCase
import numpy as np
import unittest


class TestGroundPosition(TestCase):

    def setUp(self):
        """Test fixure for test method execution."""

        self.name = "Toronto"
        self.lat = 43.6532
        self.lon = -79.3832

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

        ground_pos = GroundPosition(self.lat, self.lon)

        a = 6378.1370
        b = 6356.7523142
        lat = np.radians(ground_pos.lat)
        clat, slat = np.cos(lat), np.sin(lat)
        num = (a ** 2 * clat) ** 2 + (b ** 2 * slat) ** 2
        denom = (a * clat) ** 2 + (b * slat) ** 2
        radius = np.sqrt(num / denom)

        self.assertAlmostEqual(radius, ground_pos.radius, delta=0.001)


if __name__ == "__main__":
    unittest.main()