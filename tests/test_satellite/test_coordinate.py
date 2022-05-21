

from celest.encounter.groundposition import GroundPosition
from celest.satellite.coordinate import Coordinate
from polare import Stroke
from unittest import TestCase
import numpy as np
import unittest


WGS84_MAJOR_AXIS = 6378.137
WGS84_MINOR_AXIS = 6356.752314245


class TestCoordinate(TestCase):

    def setUp(self):

        fname = "tests/test_data/coordinate_validation_long.txt"
        cols = (0, 5, 6, 7, 11, 12, 13, 14, 15, 16, 17, 18, 19)
        skiprows = 1
        max_rows = 5000
        data = np.loadtxt(fname=fname, usecols=cols, skiprows=skiprows,
                          max_rows=max_rows)

        self.times = data[:, 0]
        self.geographical = data[:, 1:3]
        self.altitude = data[:, 3]
        self.GCRS = data[:, 4:7]
        self.GCRS_vel = data[:, 7:10]
        self.ITRS = data[:, 10:]
        self._length = data.shape[0]

        self.tri_geo = np.concatenate((self.geographical,
                                       self.altitude.reshape((-1, 1))), axis=1)

        self.offset = 2430000
        self.coor_gcrs = Coordinate(self.GCRS, self.GCRS_vel, "gcrs",
                                    self.times, self.offset)
        self.coor_itrs = Coordinate(self.ITRS, self.GCRS_vel, "itrs",
                                    self.times, self.offset)

    def test_set_base_data(self):

        self.assertIsNotNone(self.coor_gcrs._julian)
        self.assertIsNotNone(self.coor_gcrs._GCRS)
        self.assertIsNotNone(self.coor_gcrs._ITRS)
        self.assertEqual(self._length, self.coor_gcrs._length)

        self.assertIsNotNone(self.coor_itrs._julian)
        self.assertIsNotNone(self.coor_itrs._GCRS)
        self.assertIsNotNone(self.coor_itrs._ITRS)
        self.assertEqual(self._length, self.coor_itrs._length)

    def test_geo_to_itrs(self):

        test_itrs = Coordinate._geo_to_itrs(self, self.tri_geo)

        self.assertTrue(np.allclose(test_itrs[:, 0], self.ITRS[:, 0],
                                    atol=0.001))
        self.assertTrue(np.allclose(test_itrs[:, 1], self.ITRS[:, 1],
                                    atol=0.001))
        self.assertTrue(np.allclose(test_itrs[:, 2], self.ITRS[:, 2],
                                    atol=0.001))

    def test_itrs_to_geo(self):

        itrs = np.array([Stroke(self.times, self.ITRS[:, i]) for i in range(3)])
        test_lat, test_lon = Coordinate._itrs_to_geo(self, itrs)

        self.assertTrue(np.allclose(test_lat(self.times),
                        self.geographical[:, 0], atol=0.18))
        self.assertTrue(np.allclose(test_lon(self.times),
                        self.geographical[:, 1], atol=0.00001))

    def test_geo(self):

        lat1, lon1, alt1 = self.coor_gcrs.geo()
        lat2, lon2, alt2 = self.coor_itrs.geo()

        for i in range(self._length):

            eps = 1
            calc_lon, test_lon = lon1[i], self.geographical[i, 1]
            cond_one = abs(180 - calc_lon) < eps and abs(180 + test_lon) < eps
            cond_two = abs(180 + calc_lon) < eps and abs(180 - test_lon) < eps

            if cond_one:
                test_lon = (test_lon + 360) % 360
            elif cond_two:
                calc_lon = (calc_lon + 360) % 360

            self.assertAlmostEqual(lat1[i], self.geographical[i, 0], delta=0.5)
            self.assertAlmostEqual(calc_lon, test_lon, delta=0.9)
            self.assertAlmostEqual(alt1[i], self.altitude[i], delta=0.1)

        self.assertTrue(np.allclose(lat2, self.geographical[:, 0], atol=0.18))
        self.assertTrue(np.allclose(lon2, self.geographical[:, 1], atol=0.00001))
        self.assertTrue(np.allclose(alt2, self.altitude, atol=0.06))

    def test_era(self):

        julian = [2454545]
        position = [[6343.82, -2640.87, -11.26]]
        velocity = [[0, 0, 0]]
        coor = Coordinate(position, velocity, "itrs", julian)

        true_era = np.degrees(np.array([6.2360075]))
        test_era = coor.era()

        self.assertAlmostEqual(true_era[0], test_era[0], delta=0.01)

    def test_gcrs_and_itrs(self):

        calc_itrs = self.coor_gcrs._gcrs_and_itrs(self.GCRS, frame="gcrs")
        calc_gcrs = self.coor_itrs._gcrs_and_itrs(self.ITRS, frame="itrs")

        self.assertTrue(np.allclose(calc_itrs[:, 0], self.ITRS[:, 0], atol=0.35))
        self.assertTrue(np.allclose(calc_itrs[:, 1], self.ITRS[:, 1], atol=0.35))
        self.assertTrue(np.allclose(calc_itrs[:, 2], self.ITRS[:, 2], atol=0.35))

        self.assertTrue(np.allclose(calc_gcrs[:, 0], self.GCRS[:, 0], atol=0.35))
        self.assertTrue(np.allclose(calc_gcrs[:, 1], self.GCRS[:, 1], atol=0.35))
        self.assertTrue(np.allclose(calc_gcrs[:, 2], self.GCRS[:, 2], atol=0.35))

    def test_gcrs(self):

        calc_gcrs_x, calc_gcrs_y, calc_gcrs_z, _, _, _ = self.coor_itrs.gcrs()

        self.assertTrue(np.allclose(calc_gcrs_x, self.GCRS[:, 0], atol=0.35))
        self.assertTrue(np.allclose(calc_gcrs_y, self.GCRS[:, 1], atol=0.35))
        self.assertTrue(np.allclose(calc_gcrs_z, self.GCRS[:, 2], atol=0.35))

    def test_itrs(self):

        calc_itrs_x, calc_itrs_y, calc_itrs_z, _, _, _ = self.coor_gcrs.itrs()

        self.assertTrue(np.allclose(calc_itrs_x, self.ITRS[:, 0], atol=0.35))
        self.assertTrue(np.allclose(calc_itrs_y, self.ITRS[:, 1], atol=0.35))
        self.assertTrue(np.allclose(calc_itrs_z, self.ITRS[:, 2], atol=0.35))

    def test_lvlh(self):

        lvlh_data = self.coor_gcrs.lvlh()
        self.assertIsNotNone(lvlh_data)

    def test_get_ang(self):

        vector_one = np.array([[56, 92, 76], [9238, 8479, 9387], [2, 98, 23]])
        vector_two = np.array([[36, 29, 38], [2703, 947, 8739], [9827, 921, 1]])

        true_angle = np.array([16.28, 37, 83.65])
        test_angle = Coordinate._get_ang(self, vector_one, vector_two)

        self.assertTrue(np.allclose(test_angle, true_angle, atol=0.01))

    def test_horizontal(self):

        from astropy.coordinates import SkyCoord, GCRS, EarthLocation, AltAz
        from astropy import units as u
        from astropy import time

        latitude, longitude, height = 52.1579, -106.6702, 0.482
        location = EarthLocation.from_geodetic(longitude * u.deg,
                                               latitude * u.deg,
                                               height * u.km)

        times = time.Time(self.times + self.offset, format="jd")
        x, y, z = self.GCRS.T

        gcrs = SkyCoord(x=x, y=y, z=z, unit='km', frame=GCRS(obstime=times),
                        representation_type='cartesian')
        el_az = gcrs.transform_to(AltAz(obstime=times, location=location))

        true_elevation = el_az.alt.degree
        true_azimuth = el_az.az.degree

        c = Coordinate(self.GCRS, self.GCRS_vel, "gcrs", self.times, self.offset)
        location = GroundPosition(latitude, longitude, height)
        test_elevation, test_azimuth = c.horizontal(location)

        self.assertTrue(np.allclose(true_elevation, test_elevation, atol=0.32))
        self.assertTrue(np.allclose(true_azimuth, test_azimuth, atol=7.3))

    def test_off_nadir(self):
        """
        Notes
        -----
        Test cases generated using a geometric model designed by Mingde Yin.

        The tolerance required to pass all test cases is artificially inflated
        due to the low precision in the model's output and parameters used to
        generate the test cases.
        """

        location = GroundPosition(52.1579, -106.6702, 0.482)
        coor = Coordinate(self.GCRS[210:220], self.GCRS_vel[210:220], "gcrs",
                          self.times[210:220], self.offset)

        true_off_nadir = np.array([66.88, 65.09, 63.90, 63.22, 62.46, 61.67,
                                   58.42, 38.27, 23.73, 56.29])
        test_off_nadir = coor.off_nadir(location)

        self.assertTrue(np.allclose(test_off_nadir, true_off_nadir, atol=0.25))

    def test_WGS84_radius(self):
        """
        Notes
        -----
        Validation against online calculator methodology.[1]_

        References
        ----------
        .. [1] Timur. Earth Radius by Latitude (WGS 84). 2018. url:
           https://planetcalc.com/7721/.
        """

        latitude = np.radians(self.geographical[:, 0])
        numerator = (WGS84_MAJOR_AXIS ** 2 * np.cos(latitude)) ** 2 + \
            (WGS84_MINOR_AXIS ** 2 * np.sin(latitude)) ** 2
        denominator = (WGS84_MAJOR_AXIS * np.cos(latitude)) ** 2 + \
            (WGS84_MINOR_AXIS * np.sin(latitude)) ** 2
        true_radius = np.sqrt(numerator / denominator)

        latitude = Stroke(self.times, self.geographical[:, 0], "cubic")
        test_radius = Coordinate._WGS84_radius(self, latitude)(self.times)

        print(type(test_radius))

        self.assertTrue(np.allclose(true_radius, test_radius, atol=0.01))

    def test_altitude(self):

        test_altitude = self.coor_itrs.altitude()
        self.assertTrue(np.allclose(test_altitude, self.altitude, atol=0.06))

    def test_distance(self):
        """Test `Coordinate.distance`."""

        julian = [30462.5, 30462.50069, 30462.50172, 30462.50279, 30462.50386]
        position = [[6343.81620, -2640.87223, -11.25542802],
                    [6295.64584, -2718.09271, 443.08233],
                    [6173.04658, -2808.91831, 1108.58544],
                    [5980.91111, -2872.96258, 1792.17964],
                    [5724.30203, -2904.09810, 2460.25799]]
        velocity = np.zeros((5, 3))
        location = GroundPosition(52.1579, -106.6702, 0.482)
        coor = Coordinate(position, velocity, "itrs", julian, self.offset)

        true_distance = [9070.49269, 8776.71795, 8330.99543, 7851.70083, 7359.09190]
        test_distance = coor.distance(location)

        self.assertTrue(np.allclose(true_distance, test_distance, atol=0.35))


if __name__ == "__main__":
    unittest.main()
