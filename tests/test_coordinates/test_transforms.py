

from celest.coordinates.gcrs import GCRS
from celest.coordinates.ground_location import GroundLocation
from celest.coordinates.itrs import ITRS
from celest.coordinates.wgs84 import WGS84
from celest.coordinates.transforms import (
    _gcrs_to_itrs,
    _itrs_to_gcrs,
    _itrs_to_wgs84,
    _wgs84_to_itrs,
    _itrs_to_azel
)
from celest import units as u
from unittest import TestCase
import numpy as np


class TestCoordinateTransformations(TestCase):

    def setUp(self):
        file_name = "tests/test_data/coordinate_validation_long.txt"
        cols = (0, 5, 6, 7, 11, 12, 13, 17, 18, 19)
        skiprows = 1
        max_rows = 5000
        data = np.loadtxt(fname=file_name, usecols=cols, skiprows=skiprows,
                          max_rows=max_rows)

        julian = data[:, 0] + 2430000
        wgs84 = data[:, 1:3]
        altitude = data[:, 3]
        gcrs = data[:, 4:7]
        itrs = data[:, 7:]

        self.wgs84 = WGS84(julian, wgs84[:, 0], wgs84[:, 1], altitude, u.deg,
                           u.km)
        self.gcrs = GCRS(julian, gcrs[:, 0], gcrs[:, 1], gcrs[:, 2], u.km)
        self.itrs = ITRS(julian, itrs[:, 0], itrs[:, 1], itrs[:, 2], u.km)

    def test_gcrs_to_itrs_for_validation(self):
        itrs = _gcrs_to_itrs(self.gcrs)

        self.assertIsInstance(itrs, ITRS)
        self.assertTrue(np.allclose(itrs.x.data, self.itrs.x.data, atol=0.35))
        self.assertTrue(np.allclose(itrs.y.data, self.itrs.y.data, atol=0.35))
        self.assertTrue(np.allclose(itrs.z.data, self.itrs.z.data, atol=0.35))

    def test_itrs_to_gcrs_for_validation(self):
        gcrs = _itrs_to_gcrs(self.itrs)

        self.assertTrue(np.allclose(gcrs.x.data, self.gcrs.x.data, atol=0.35))
        self.assertTrue(np.allclose(gcrs.y.data, self.gcrs.y.data, atol=0.35))
        self.assertTrue(np.allclose(gcrs.z.data, self.gcrs.z.data, atol=0.35))

    def test_itrs_to_wgs84_for_validation(self):
        wgs84 = _itrs_to_wgs84(self.itrs)

        self.assertTrue(np.allclose(wgs84.latitude.to(u.deg).data,
                                    self.wgs84.latitude.to(u.deg).data,
                                    atol=0.18))
        self.assertTrue(np.allclose(wgs84.longitude.to(u.deg).data,
                                    self.wgs84.longitude.to(u.deg).data,
                                    atol=0.00001))
        self.assertTrue(np.allclose(wgs84.height.to(u.km).data,
                                    self.wgs84.height.to(u.km).data,
                                    atol=0.06))

    def test_wgs84_to_itrs_for_validation(self):
        itrs = _wgs84_to_itrs(self.wgs84)

        self.assertTrue(np.allclose(itrs.x.data, self.itrs.x.data, atol=0.001))
        self.assertTrue(np.allclose(itrs.y.data, self.itrs.y.data, atol=0.001))
        self.assertTrue(np.allclose(itrs.z.data, self.itrs.z.data, atol=0.001))

    def test_itrs_to_azel_for_validation(self):
        from astropy.coordinates import SkyCoord, ITRS, EarthLocation, AltAz
        from astropy import units as astropy_u
        from astropy.time import Time

        latitude, longitude, height = 52.1579, -106.6702, 0.482
        location = EarthLocation.from_geodetic(longitude * astropy_u.deg,
                                               latitude * astropy_u.deg,
                                               height * astropy_u.km)

        times = Time(self.itrs.time.data, format="jd")
        itrs = SkyCoord(x=self.itrs.x.data, y=self.itrs.y.data,
                        z=self.itrs.z.data, unit='km', frame=ITRS(obstime=times),
                        representation_type='cartesian')
        expected_azel = itrs.transform_to(AltAz(obstime=times, location=location))

        location = GroundLocation(latitude, longitude, height, u.deg, u.km)
        actual_azel = _itrs_to_azel(self.itrs, location)

        self.assertTrue(np.allclose(expected_azel.alt.degree,
                                    actual_azel.elevation.data, atol=0.32))
        self.assertTrue(np.allclose(expected_azel.az.degree,
                                    actual_azel.azimuth.data, atol=7.3))
