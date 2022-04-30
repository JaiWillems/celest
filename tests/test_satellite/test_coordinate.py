

from celest.encounter.groundposition import GroundPosition
from celest.satellite.coordinate import Coordinate
from polare import Stroke
from unittest import TestCase
import numpy as np
import unittest


class TestCoordinate(TestCase):

    def setUp(self):

        fname = "tests/test_data/coordinate_validation_set.txt"
        data = np.loadtxt(fname=fname, delimiter="\t", skiprows=1)

        self.times = data[:, 0]
        self.geo = data[:, 1:3]
        self.alt = data[:, 3]
        self.local_altaz = data[:, 4:7]
        self.GCRS = data[:, 7:10]
        self.ITRS = data[:, 10:]
        self._length = data.shape[0]

        self.tri_geo = np.concatenate((self.geo, self.alt.reshape((-1, 1))), axis=1)

        self.offset = 2430000
        self.coor_gcrs = Coordinate(self.GCRS, "gcrs", self.times, self.offset)
        self.coor_itrs = Coordinate(self.ITRS, "itrs", self.times, self.offset)
        self.coor_geo = Coordinate(self.tri_geo, "geo", self.times, self.offset)

    def test_set_base_position(self):
        """Test `Coordinate._set_base_position`."""

        coor_1 = Coordinate(self.geo, "geo", self.times, self.offset)

        self.assertIsNotNone(coor_1._julian)
        self.assertIsNone(coor_1._GCRS)
        self.assertIsNotNone(coor_1._ITRS)
        self.assertEqual(self._length, coor_1._length)

        coor_2 = Coordinate(self.GCRS, "gcrs", self.times, self.offset)

        self.assertIsNotNone(coor_2._julian)
        self.assertIsNotNone(coor_2._GCRS)
        self.assertIsNone(coor_2._ITRS)
        self.assertEqual(self._length, coor_2._length)

        coor_3 = Coordinate(self.ITRS, "itrs", self.times, self.offset)

        self.assertIsNotNone(coor_3._julian)
        self.assertIsNone(coor_3._GCRS)
        self.assertIsNotNone(coor_3._ITRS)
        self.assertEqual(self._length, coor_3._length)

    def test_geo_to_itrs(self):
        """Test `Coordinate._geo_to_itrs`.

        Notes
        -----
        Test cases generated using a GMAT data set.
        """

        calc_itrs = self.coor_geo._geo_to_itrs(self.tri_geo)

        for i in range(self._length):
            with self.subTest(i=i):
                self.assertAlmostEqual(calc_itrs[i, 0], self.ITRS[i, 0],
                                       delta=0.001)
                self.assertAlmostEqual(calc_itrs[i, 1], self.ITRS[i, 1],
                                       delta=0.001)
                self.assertAlmostEqual(calc_itrs[i, 2], self.ITRS[i, 2],
                                       delta=0.001)

    def test_itrs_to_geo(self):
        """Test `Coordinate._itrs_to_geo`.

        Notes
        -----
        Test cases generated using a GMAT data set.
        """

        itrs = np.array([Stroke(self.times, self.ITRS[:, i]) for i in range(3)])
        calc_lat, calc_lon = self.coor_itrs._itrs_to_geo(itrs)

        for i in range(self._length):
            with self.subTest(i=i):
                self.assertAlmostEqual(calc_lat(self.times[i]), self.geo[i, 0],
                                       delta=0.18)
                self.assertAlmostEqual(calc_lon(self.times[i]), self.geo[i, 1],
                                       delta=0.00001)

    def test_geo(self):
        """Test `Coordinate.geo`."""

        lat1, lon1, alt1 = self.coor_gcrs.geo()
        lat2, lon2, alt2 = self.coor_itrs.geo()

        for i in range(self._length):
            with self.subTest(i=i):

                eps = 1
                calc_lon, test_lon = lon1[i], self.geo[i, 1]
                cond_one = abs(180 - calc_lon) < eps and abs(180 + test_lon) < eps
                cond_two = abs(180 + calc_lon) < eps and abs(180 - test_lon) < eps

                if cond_one:
                    test_lon = (test_lon + 360) % 360
                elif cond_two:
                    calc_lon = (calc_lon + 360) % 360

                self.assertAlmostEqual(lat1[i], self.geo[i, 0], delta=0.5)
                self.assertAlmostEqual(calc_lon, test_lon, delta=0.9)
                self.assertAlmostEqual(alt1[i], self.alt[i], delta=0.1)

        for i in range(self._length):
            with self.subTest(i=i):
                self.assertAlmostEqual(lat2[i], self.geo[i, 0], delta=0.18)
                self.assertAlmostEqual(lon2[i], self.geo[i, 1], delta=0.00001)
                self.assertAlmostEqual(alt2[i], self.alt[i], delta=0.06)

    def test_era(self):
        """Test `Coordinate.era`.

        Notes
        -----
        `Coordinate.era` was validated through `Coordinate._gcrs_and_itrs`
        validation. Test cases were generated from the `Coordinate.era` output
        after it was shown to be correct to ensure functionality is consistent.
        """

        era = np.degrees(np.array([6.2360075]))

        julian = [2454545]
        position = [[6343.82, -2640.87, -11.26]]
        coor = Coordinate(position, "itrs", julian)
        calc_era = coor.era()

        self.assertAlmostEqual(era[0], calc_era[0], delta=0.01)

    def test_gcrs_and_itrs(self):
        """Test `Coordinate._gcrs_and_itrs`.

        Notes
        -----
        Test cases generated using a GMAT data set.
        """

        calc_itrs = self.coor_gcrs._gcrs_and_itrs(self.GCRS, frame="gcrs")
        calc_gcrs = self.coor_itrs._gcrs_and_itrs(self.ITRS, frame="itrs")

        for i in range(self._length):
            with self.subTest(i=i):
                self.assertAlmostEqual(calc_itrs[i, 0], self.ITRS[i, 0],
                                       delta=0.35)
                self.assertAlmostEqual(calc_itrs[i, 1], self.ITRS[i, 1],
                                       delta=0.35)
                self.assertAlmostEqual(calc_itrs[i, 2], self.ITRS[i, 2],
                                       delta=0.35)

        for i in range(self._length):
            with self.subTest(i=i):
                self.assertAlmostEqual(calc_gcrs[i, 0], self.GCRS[i, 0],
                                       delta=0.35)
                self.assertAlmostEqual(calc_gcrs[i, 1], self.GCRS[i, 1],
                                       delta=0.35)
                self.assertAlmostEqual(calc_gcrs[i, 2], self.GCRS[i, 2],
                                       delta=0.35)

    def test_gcrs(self):
        """Test `Coordinate.gcrs`."""

        calc_gcrs_x1, calc_gcrs_y1, calc_gcrs_z1 = self.coor_itrs.gcrs()
        calc_gcrs_x2, calc_gcrs_y2, calc_gcrs_z2 = self.coor_geo.gcrs()

        for i in range(self._length):
            with self.subTest(i=i):
                self.assertAlmostEqual(calc_gcrs_x1[i], self.GCRS[i, 0],
                                       delta=0.35)
                self.assertAlmostEqual(calc_gcrs_y1[i], self.GCRS[i, 1],
                                       delta=0.35)
                self.assertAlmostEqual(calc_gcrs_z1[i], self.GCRS[i, 2],
                                       delta=0.35)

        for i in range(self._length):
            with self.subTest(i=i):
                self.assertAlmostEqual(calc_gcrs_x2[i], self.GCRS[i, 0],
                                       delta=0.35)
                self.assertAlmostEqual(calc_gcrs_y2[i], self.GCRS[i, 1],
                                       delta=0.35)
                self.assertAlmostEqual(calc_gcrs_z2[i], self.GCRS[i, 2],
                                       delta=0.35)

    def test_itrs(self):
        """Test `Coordinate.itrs`."""

        calc_itrs_x1, calc_itrs_y1, calc_itrs_z1 = self.coor_gcrs.itrs()
        calc_itrs_x2, calc_itrs_y2, calc_itrs_z2 = self.coor_geo.itrs()

        for i in range(self._length):
            with self.subTest(i=i):
                self.assertAlmostEqual(calc_itrs_x1[i], self.ITRS[i, 0],
                                       delta=0.35)
                self.assertAlmostEqual(calc_itrs_y1[i], self.ITRS[i, 1],
                                       delta=0.35)
                self.assertAlmostEqual(calc_itrs_z1[i], self.ITRS[i, 2],
                                       delta=0.35)

        for i in range(self._length):
            with self.subTest(i=i):
                self.assertAlmostEqual(calc_itrs_x2[i], self.ITRS[i, 0],
                                       delta=0.35)
                self.assertAlmostEqual(calc_itrs_y2[i], self.ITRS[i, 1],
                                       delta=0.35)
                self.assertAlmostEqual(calc_itrs_z2[i], self.ITRS[i, 2],
                                       delta=0.35)

    def test_get_ang(self):
        """Test `Coordinate._get_ang`."""

        vec_one = np.array([[56, 92, 76], [9238, 8479, 9387], [2, 98, 23]])
        vec_two = np.array([[36, 29, 38], [2703, 947, 8739], [9827, 921, 1]])
        ang = np.array([16.28, 37, 83.65])

        calc_ang = self.coor_geo._get_ang(vec_one, vec_two)

        for i in range(calc_ang.size):
            with self.subTest(i=i):
                self.assertAlmostEqual(ang[i], calc_ang[i], delta=0.01)

    def test_horizontal(self):
        """Test `Coordinate.horizontal`.

        Notes
        -----
        Test cases generated using the `Astropy` Python package.
        """

        from astropy.coordinates import SkyCoord, ITRS, GCRS, EarthLocation, AltAz
        from astropy import units as u
        from astropy import time

        # Set up observer location.
        lat, lon = 52.1579, -106.6702
        loc = EarthLocation.from_geodetic(lon*u.deg, lat*u.deg)

        # Prepare time and position information.
        times = time.Time(self.times + self.offset, format="jd")
        x, y, z = self.GCRS.T

        # Define coordinate frames.
        gcrs = GCRS(obstime=times)
        itrs = ITRS(obstime=times)
        altaz = AltAz(obstime=times, location=loc)

        # Convert between frames.
        gcrsCoor = SkyCoord(x=x, y=y, z=z, unit='km',
                            frame=gcrs, representation_type='cartesian')
        itrsCoor = gcrsCoor.transform_to(itrs)
        altazCoor = itrsCoor.transform_to(altaz)

        # Get validation data.
        itrsData = np.transpose(np.array([itrsCoor.x, itrsCoor.y, itrsCoor.z]))
        alt = altazCoor.alt.degree
        az = altazCoor.az.degree

        # Get Celest results.
        coor = Coordinate(itrsData, "itrs", self.times, self.offset)
        groundPos = GroundPosition(lat, lon)
        calc_alt, calc_az = coor.horizontal(groundPos)

        for i in range(calc_alt.size):
            with self.subTest(i=i):
                self.assertAlmostEqual(alt[i], calc_alt[i], delta=0.25)
                self.assertAlmostEqual(az[i], calc_az[i], delta=7.5)

    def test_off_nadir(self):
        """Test `Coordinate.off_nadir`.

        Notes
        -----
        Test cases generated using a geometric model designed by Mingde Yin.

        The tolerance required to pass all test cases is artificially inflated
        due to the low precision in the model's output and parameters used to
        generate the test cases.
        """

        off_nadir = np.array([66.88, 65.09, 63.90, 63.22, 62.46, 61.67, 58.42,
                              38.27, 23.73, 56.29])

        location = GroundPosition(52.1579, -106.6702)

        julian = self.times[210:220]
        coor = Coordinate(self.ITRS[210:220], "itrs", julian, self.offset)
        calc_off_nadir = coor.off_nadir(location)

        for i in range(10):
            with self.subTest(i=i):
                self.assertAlmostEqual(off_nadir[i], calc_off_nadir[i], delta=0.3)

    def test_WGS84_radius(self):
        """Test `Coordinate._WGS84_radius`.

        Notes
        -----
        Validation against online calculator methodology.[1]_

        References
        ----------
        .. [1] Timur. Earth Radius by Latitude (WGS 84). 2018. url:
           https://planetcalc.com/7721/.
        """

        a, b = 6378.1370, 6356.7523142
        lat = np.radians(self.geo[:, 0])
        clat, slat = np.cos(lat), np.sin(lat)
        num = (a ** 2 * clat) ** 2 + (b ** 2 * slat) ** 2
        denom = (a * clat) ** 2 + (b * slat) ** 2
        radius = np.sqrt(num / denom)

        lat = Stroke(self.times, self.geo[:, 0], "cubic")

        calc_radius = self.coor_geo._WGS84_radius(lat)

        for i in range(radius.size):
            with self.subTest(i=i):
                self.assertAlmostEqual(radius[i], calc_radius(self.times[i]),
                                       delta=0.001)

    def test_altitude(self):
        """Test `Coordinate.altitude`.

        Notes
        -----
        Test cases are taken from a GMAT data set.
        """

        calc_alt = self.coor_itrs.altitude()

        for i in range(calc_alt.size):
            with self.subTest(i=i):
                self.assertAlmostEqual(self.alt[i], calc_alt[i], delta=0.06)

    def test_distance(self):
        """Test `Coordinate.distance`."""

        julian = [30462.5, 30462.50069444, 30462.50171802, 30462.50278711,
                  30462.50386162]
        position = [[6343.81620221, -2640.87223125, -11.25541802],
                    [6295.64583763, -2718.09271472, 443.08232543],
                    [6173.04658005, -2808.91831102, 1108.5854422],
                    [5980.91111229, -2872.96257946, 1792.17964249],
                    [5724.3020284, -2904.09809986, 2460.25799377]]
        dist = [9070.49268746, 8776.7179543, 8330.99543153, 7851.70082642,
                7359.09189844]

        location = GroundPosition(52.1579, -106.6702)

        coor = Coordinate(position, "itrs", julian, self.offset)
        calc_dist = coor.distance(location)

        for i in range(5):
            with self.subTest(i=i):
                self.assertAlmostEqual(calc_dist[i], dist[i], delta=0.1)


if __name__ == "__main__":
    unittest.main()
