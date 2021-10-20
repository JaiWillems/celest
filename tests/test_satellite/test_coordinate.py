

from celest.encounter.groundposition import GroundPosition
from celest.satellite.coordinate import Coordinate
from celest.satellite.time import Time
from unittest import TestCase
import numpy as np
import unittest


class TestCoordinate(TestCase):

    def setUp(self):
        """Test fixure for test method execution."""

        fname = "tests/test_data/coordinate_validation_set.txt"
        data = np.loadtxt(fname=fname, delimiter="\t", skiprows=1)

        self.times = data[:, 0]
        self.geo = data[:, 1:3]
        self.alt = data[:, 3]
        self.local_altaz = data[:, 4:7]
        self.ECI = data[:, 7:10]
        self.ECEF = data[:, 10:]
        self.length = data.shape[0]

        geo = np.concatenate((self.geo, self.alt.reshape((-1, 1))), axis=1)

        self.timeData = Time(self.times, 2430000)
        self.coor_ECI = Coordinate(self.ECI, "eci", self.timeData)
        self.coor_ECEF = Coordinate(self.ECEF, "ecef", self.timeData)
        self.coor_GEO = Coordinate(geo, "geo", self.timeData)

    def test_set_base_position(self):
        """Test `Coordinate._set_base_position`."""

        coor_1 = Coordinate(self.geo, "geo", self.timeData)

        self.assertEqual(coor_1._GEO.shape[1], 3)

        self.assertIsNotNone(coor_1.time)
        self.assertIsNotNone(coor_1._GEO)
        self.assertIsNone(coor_1._ECI)
        self.assertIsNone(coor_1._ECEF)
        self.assertEqual(self.length, coor_1.length)

        coor_2 = Coordinate(self.ECI, "eci", self.timeData)

        self.assertIsNotNone(coor_2.time)
        self.assertIsNone(coor_2._GEO)
        self.assertIsNotNone(coor_2._ECI)
        self.assertIsNone(coor_2._ECEF)
        self.assertEqual(self.length, coor_2.length)

        coor_3 = Coordinate(self.ECEF, "ecef", self.timeData)

        self.assertIsNotNone(coor_3.time)
        self.assertIsNone(coor_3._GEO)
        self.assertIsNone(coor_3._ECI)
        self.assertIsNotNone(coor_3._ECEF)
        self.assertEqual(self.length, coor_3.length)

    def test_GEO_to_ECEF(self):
        """Test `Coordinate._GEO_to_ECEF`.

        Notes
        -----
        Test cases are taken from a GMAT data set.
        """

        geo = np.concatenate((self.geo, self.alt.reshape((-1, 1))), axis=1)
        calc_ECEF = self.coor_GEO._GEO_to_ECEF(geo)

        for i in range(self.length):
            with self.subTest(i=i):
                self.assertAlmostEqual(calc_ECEF[i, 0], self.ECEF[i, 0], delta=0.001)

        for i in range(self.length):
            with self.subTest(i=i):
                self.assertAlmostEqual(calc_ECEF[i, 1], self.ECEF[i, 1], delta=0.001)

        for i in range(self.length):
            with self.subTest(i=i):
                self.assertAlmostEqual(calc_ECEF[i, 2], self.ECEF[i, 2], delta=0.001)

    def test_ECEF_to_GEO(self):
        """Test `Coordinate._ECEF_to_GEO`.

        Notes
        -----
        Test cases are taken from a GMAT data set.
        """

        calc_GEO = self.coor_ECEF._ECEF_to_GEO(self.ECEF)

        for i in range(self.length):
            with self.subTest(i=i):
                self.assertAlmostEqual(calc_GEO[i, 0], self.geo[i, 0], delta=0.18)

        for i in range(self.length):
            with self.subTest(i=i):
                self.assertAlmostEqual(calc_GEO[i, 1], self.geo[i, 1], delta=0.00001)

        for i in range(self.length):
            with self.subTest(i=i):
                self.assertAlmostEqual(calc_GEO[i, 2], self.alt[i], delta=0.06)

    def test_GEO(self):
        """Test `Coordinate.GEO`."""

        calc_geo_1 = self.coor_ECI.GEO()
        calc_geo_2 = self.coor_ECEF.GEO()

        for i in range(self.length):
            with self.subTest(i=i):

                epsilon = 1
                calc_lon, test_lon = calc_geo_1[i, 1], self.geo[i, 1]
                condition_one = abs(180 - calc_lon) < epsilon and abs(180 + test_lon) < epsilon
                condition_two = abs(180 + calc_lon) < epsilon and abs(180 - test_lon) < epsilon
                if condition_one:
                    test_lon = (test_lon + 360) % 360
                elif condition_two:
                    calc_lon = (calc_lon + 360) % 360

                self.assertAlmostEqual(calc_geo_1[i, 0], self.geo[i, 0], delta=0.5)
                self.assertAlmostEqual(calc_lon, test_lon, delta=0.9)
                self.assertAlmostEqual(calc_geo_1[i, 2], self.alt[i], delta=0.1)

        for i in range(self.length):
            with self.subTest(i=i):
                self.assertAlmostEqual(calc_geo_2[i, 0], self.geo[i, 0], delta=0.18)
                self.assertAlmostEqual(calc_geo_2[i, 1], self.geo[i, 1], delta=0.00001)
                self.assertAlmostEqual(calc_geo_2[i, 2], self.alt[i], delta=0.06)

    def test_ERA(self):
        """Test `Coordinate.ERA`.

        Notes
        -----
        The `Coordinate.ERA` method was validated through
        `Coordinate._ECI_and_ECEF` validation. This test method was implemented
        using the calculated output of the `Coordinate.ERA` method after it was
        shown to be correct.
        """

        era = np.array([6.2360075])

        timeData = Time(np.array([2454545]))
        basePos = np.array([[6343.82, -2640.87, -11.26]])
        posData = Coordinate(basePos, "ecef", timeData)
        calc_era = posData.ERA()

        self.assertAlmostEqual(era[0], calc_era[0], delta=0.01)

    def test_ECI_and_ECEF(self):
        """Test `Coordinate._ECI_and_ECEF`.

        Notes
        -----
        Test cases are taken from a GMAT data set.
        """

        calc_ECEF = self.coor_ECI._ECI_and_ECEF(self.ECI, frame="eci")
        calc_ECI = self.coor_ECEF._ECI_and_ECEF(self.ECEF, frame="ecef")

        for i in range(self.length):
            with self.subTest(i=i):
                self.assertAlmostEqual(calc_ECEF[i, 0], self.ECEF[i, 0], delta=20)
                self.assertAlmostEqual(calc_ECEF[i, 1], self.ECEF[i, 1], delta=20)
                self.assertAlmostEqual(calc_ECEF[i, 2], self.ECEF[i, 2], delta=20)

        for i in range(self.length):
            with self.subTest(i=i):
                self.assertAlmostEqual(calc_ECI[i, 0], self.ECI[i, 0], delta=20)
                self.assertAlmostEqual(calc_ECI[i, 1], self.ECI[i, 1], delta=20)
                self.assertAlmostEqual(calc_ECI[i, 2], self.ECI[i, 2], delta=20)

    def test_ECI(self):
        """Test `Coordinate.ECI`."""

        calc_eci_1 = self.coor_ECEF.ECI()
        calc_eci_2 = self.coor_GEO.ECI()

        for i in range(self.length):
            with self.subTest(i=i):
                self.assertAlmostEqual(calc_eci_1[i, 0], self.ECI[i, 0], delta=17)
                self.assertAlmostEqual(calc_eci_1[i, 1], self.ECI[i, 1], delta=17)
                self.assertAlmostEqual(calc_eci_1[i, 2], self.ECI[i, 2], delta=17)

        for i in range(self.length):
            with self.subTest(i=i):
                self.assertAlmostEqual(calc_eci_2[i, 0], self.ECI[i, 0], delta=17)
                self.assertAlmostEqual(calc_eci_2[i, 1], self.ECI[i, 1], delta=17)
                self.assertAlmostEqual(calc_eci_2[i, 2], self.ECI[i, 2], delta=17)

    def test_ECEF(self):
        """Test `Coordinate.ECEF`."""

        calc_ecef_1 = self.coor_ECI.ECEF()
        calc_ecef_2 = self.coor_GEO.ECEF()

        for i in range(self.length):
            with self.subTest(i=i):
                self.assertAlmostEqual(calc_ecef_1[i, 0], self.ECEF[i, 0], delta=17)
                self.assertAlmostEqual(calc_ecef_1[i, 1], self.ECEF[i, 1], delta=17)
                self.assertAlmostEqual(calc_ecef_1[i, 2], self.ECEF[i, 2], delta=17)

        for i in range(self.length):
            with self.subTest(i=i):
                self.assertAlmostEqual(calc_ecef_2[i, 0], self.ECEF[i, 0], delta=5)
                self.assertAlmostEqual(calc_ecef_2[i, 1], self.ECEF[i, 1], delta=5)
                self.assertAlmostEqual(calc_ecef_2[i, 2], self.ECEF[i, 2], delta=5)

    def test_get_ang(self):
        """Test `Coordinate._get_ang`."""

        vec_one = np.array([[56, 92, 76], [9238, 8479, 9387], [2, 98, 23]])
        vec_two = np.array([[36, 29, 38], [2703, 947, 8739], [9827, 921, 1]])
        ang = np.array([16.28, 37, 83.65])

        calc_ang = self.coor_GEO._get_ang(vec_one, vec_two)

        for i in range(calc_ang.size):
            with self.subTest(i=i):
                self.assertAlmostEqual(ang[i], calc_ang[i], delta=0.01)

    def test_horizontal(self):
        """Test `Coordinate.horizontal`.

        Notes
        -----
        Test cases are generated using the `Astropy` Python package.
        """

        from astropy.coordinates import SkyCoord, ITRS, GCRS, EarthLocation, AltAz
        from astropy import units as u
        from astropy import time

        # Set up observer location.
        obsCoor = (52.1579, -106.6702)
        loc = EarthLocation.from_geodetic(obsCoor[1]*u.deg, obsCoor[0]*u.deg)

        # Prepare time and position information.
        timeData = time.Time(self.times + 2430000, format="jd")
        x, y, z = self.ECI[:, 0], self.ECI[:, 1], self.ECI[:, 2]

        # Define coordinate frames.
        gcrs = GCRS(obstime=timeData)
        itrs = ITRS(obstime=timeData)
        altaz = AltAz(obstime=timeData, location=loc)

        # Convert between frames.
        gcrsCoor = SkyCoord(x=x, y=y, z=z, unit='km', frame=gcrs, representation_type='cartesian')
        itrsCoor = gcrsCoor.transform_to(itrs)
        altazCoor = itrsCoor.transform_to(altaz)

        # Get validation data.
        itrsData = np.transpose(np.array([itrsCoor.x, itrsCoor.y, itrsCoor.z]))
        alt = altazCoor.alt.degree
        az = altazCoor.az.degree

        # Get Celest results.
        coor = Coordinate(itrsData, "ecef", self.timeData)
        groundPos = GroundPosition(obsCoor)
        calc_alt, calc_az = coor.horizontal(groundPos)

        for i in range(calc_alt.size-5000):
            with self.subTest(i=i):
                self.assertAlmostEqual(alt[i], calc_alt[i], delta=0.25)

        for i in range(calc_az.size-5000):
            with self.subTest(i=i):
                self.assertAlmostEqual(az[i], calc_az[i], delta=7.5)

    def test_off_nadir(self):
        """Test `Coordinate.off_nadir`.

        Notes
        -----
        Test cases generated from geometric model designed by Mingde Yin.

        The tolerance required to pass all test cases is artificially inflated
        due to the low number of significant figures and coarse parameter
        selection of the geometric model used to generated test cases.
        """

        off_nadir = np.array([66.88, 65.09, 63.90, 63.22, 62.46, 61.67, 58.42,
                              38.27, 23.73, 56.29])

        obsCoor = (52.1579, -106.6702)
        groundPos = GroundPosition(obsCoor)

        timeData = Time(self.times[210:220], 2430000)
        coor = Coordinate(self.ECEF[210:220], "ecef", timeData)
        calc_off_nadir = coor.off_nadir(groundPos)

        for i in range(10):
            with self.subTest(i=i):
                self.assertAlmostEqual(off_nadir[i], calc_off_nadir[i], delta=0.3)

    def test_wgs84_radius(self):
        """Test `Coordinate._WGS84_radius`.

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
        lat = np.radians(self.geo[:, 0])
        clat, slat = np.cos(lat), np.sin(lat)
        num = (a ** 2 * clat) ** 2 + (b ** 2 * slat) ** 2
        denom = (a * clat) ** 2 + (b * slat) ** 2
        radius = np.sqrt(num / denom)

        calc_radius = self.coor_GEO._WGS84_radius(self.geo[:, 0])

        for i in range(5):
            with self.subTest(i=i):
                self.assertAlmostEqual(radius[i], calc_radius[i], delta=0.001)

    def test_altitude(self):
        """Test `Coordinate.altitude`.

        Notes
        -----
        Test cases are taken from a GMAT data set.
        """

        calc_alt = self.coor_ECEF.altitude()

        for i in range(calc_alt.size):
            with self.subTest(i=i):
                self.assertAlmostEqual(self.alt[i], calc_alt[i], delta=0.06)

    def test_distance(self):
        """Test `Coordinate.distance`."""

        times = np.array([30462.5, 30462.50069444, 30462.50171802,
                          30462.50278711, 30462.50386162])
        position = np.array([[6343.81620221, -2640.87223125, -11.25541802],
                             [6295.64583763, -2718.09271472, 443.08232543],
                             [6173.04658005, -2808.91831102, 1108.5854422],
                             [5980.91111229, -2872.96257946, 1792.17964249],
                             [5724.3020284, -2904.09809986, 2460.25799377]])
        dist = np.array([9070.49268746, 8776.7179543, 8330.99543153,
                         7851.70082642, 7359.09189844])

        obsCoor = (52.1579, -106.6702)
        groundPos = GroundPosition(obsCoor)

        timeData = Time(times, 2430000)
        coor = Coordinate(position, "ecef", timeData)
        calc_dist = coor.distance(groundPos)

        for i in range(5):
            with self.subTest(i=i):
                self.assertAlmostEqual(calc_dist[i], dist[i], delta=0.1)


if __name__ == "__main__":
    unittest.main()
