

import numpy as np
import unittest
from unittest import TestCase
from celest.encounter import GroundPosition


class TestGroundPosition(TestCase):

    def setUp(self):
        """Test fixure for test method execution."""

        self.name = "Toronto"
        self.coor = (43.6532, -79.3832)
    
    def test_encounter_attributes(self):
        """Test `GroundPosition` instantiation.
        
        Notes
        -----
        Test that the appropriate parameters are set basef on the input
        encounter type.
        """

        ang = 20
        ground_pos = GroundPosition(self.name, self.coor, "image", ang)

        self.assertEqual(ground_pos.ang, ang)
        self.assertTrue(ground_pos.type == "I")
        self.assertTrue(ground_pos.ang_type == "N")
        self.assertEqual(ground_pos.lighting, 1)
        self.assertEqual(ground_pos.solar_constraint_angle, 0)

        ang = 20
        ground_pos = GroundPosition(self.name, self.coor, "data_link", ang)

        self.assertEqual(ground_pos.ang, ang)
        self.assertTrue(ground_pos.type == "T")
        self.assertTrue(ground_pos.ang_type == "A")
        self.assertEqual(ground_pos.lighting, 0)
        self.assertEqual(ground_pos.solar_constraint_angle, 30)

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

        ground_pos = GroundPosition(self.name, self.coor, "image", 0)

        a = 6378.1370
        b = 6356.7523142
        lat = np.radians(ground_pos.coor[0])
        clat, slat = np.cos(lat), np.sin(lat)
        num = (a ** 2 * clat) ** 2 + (b ** 2 * slat) ** 2
        denom = (a * clat) ** 2 + (b * slat) ** 2
        radius = np.sqrt(num / denom)

        self.assertAlmostEqual(radius, ground_pos.radius, delta=0.001)


if __name__ == "__main__":
    unittest.main()