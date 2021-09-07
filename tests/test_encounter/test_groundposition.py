

import numpy as np
import unittest
from unittest import TestCase
from celest.encounter import GroundPosition


class TestGroundPosition(TestCase):

    def setUp(self):

        name = "Toronto"
        coor = (43.6532, -79.3832)

        self.ground_pos = GroundPosition(name=name, coor=coor)

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

        a = 6378.1370
        b = 6356.7523142
        lat = np.radians(self.ground_pos.coor[0])
        clat, slat = np.cos(lat), np.sin(lat)
        num = (a ** 2 * clat) ** 2 + (b ** 2 * slat) ** 2
        denom = (a * clat) ** 2 + (b * slat) ** 2
        radius = np.sqrt(num / denom)

        self.assertAlmostEqual(radius, self.ground_pos.radius, delta=0.001)

    def test_add_encounter(self):
        """Test `GroundPosition.add_encounter`."""
        
        name = "CYYZ IMG"
        encType = "I"
        ang = 30
        angType = "N"
        solar = 0
        sca = 0

        self.ground_pos.add_encounter(name=name, encType=encType, ang=ang,
                                      angType=angType, solar=solar, sca=sca)
        
        enc = self.ground_pos.encounters[name]

        self.assertTrue(enc.name == name)
        self.assertTrue(enc.type == encType)
        self.assertTrue(enc.ang == ang)
        self.assertTrue(enc.ang_type == angType)
        self.assertTrue(enc.solar == solar)
        self.assertTrue(enc.solar_constraint_ang == sca)



if __name__ == "__main__":
    unittest.main()