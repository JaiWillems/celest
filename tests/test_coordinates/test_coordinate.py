

from celest.coordinates.coordinate import Coordinate
from celest.coordinates.frames.azel import AzEl
from celest.coordinates.frames.gcrs import GCRS
from celest.coordinates.frames.itrs import ITRS
from celest.coordinates.frames.wgs84 import WGS84
from celest.coordinates.ground_location import GroundLocation
from celest import units as u
from unittest import TestCase
import numpy as np


class TestCoordinate(TestCase):

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

        self.wgs84 = WGS84(julian, wgs84[:, 0], wgs84[:, 1], altitude, u.deg, u.km)
        self.gcrs = GCRS(julian, gcrs[:, 0], gcrs[:, 1], gcrs[:, 2], u.km)
        self.itrs = ITRS(julian, itrs[:, 0], itrs[:, 1], itrs[:, 2], u.km)

        self.wgs84_coordinate = Coordinate(self.wgs84)
        self.gcrs_coordinate = Coordinate(self.gcrs)
        self.itrs_coordinate = Coordinate(self.itrs)

        self.location = GroundLocation(43.6532, -79.3832, 76, u.deg, u.m)

    def test_initialization(self):
        self.assertIsInstance(self.wgs84_coordinate, Coordinate)
        self.assertIsInstance(self.gcrs_coordinate, Coordinate)
        self.assertIsInstance(self.itrs_coordinate, Coordinate)

    def test_value_error_raised_when_end_result_frame_is_input(self):
        azel_frame = self.itrs_coordinate.convert_to(AzEl, self.location)
        self.assertRaises(ValueError, Coordinate, azel_frame)

    def test_convert_to_same_returns_same_with_gcrs(self):
        self.assertIsInstance(self.gcrs_coordinate.convert_to(GCRS), GCRS)

    def test_convert_to_same_returns_same_with_itrs(self):
        self.assertIsInstance(self.itrs_coordinate.convert_to(ITRS), ITRS)

    def test_convert_to_same_returns_same_with_wgs84(self):
        self.assertIsInstance(self.wgs84_coordinate.convert_to(WGS84), WGS84)

    def test_gcrs_to_itrs(self):
        self.assertIsInstance(self.gcrs_coordinate.convert_to(ITRS), ITRS)

    def test_itrs_to_gcrs(self):
        self.assertIsInstance(self.itrs_coordinate.convert_to(GCRS), GCRS)

    def test_itrs_to_wgs84(self):
        self.assertIsInstance(self.itrs_coordinate.convert_to(WGS84), WGS84)

    def test_wgs84_to_itrs(self):
        self.assertIsInstance(self.wgs84_coordinate.convert_to(ITRS), ITRS)

    def test_gcrs_to_wgs84(self):
        self.assertIsInstance(self.gcrs_coordinate.convert_to(WGS84), WGS84)

    def test_wgs84_to_gcrs(self):
        self.assertIsInstance(self.wgs84_coordinate.convert_to(GCRS), GCRS)

    def test_itrs_to_azel(self):
        self.assertIsInstance(self.itrs_coordinate.convert_to(AzEl, self.location), AzEl)

    def test_itrs_to_azel_raises_error_with_no_location(self):
        self.assertRaises(ValueError, self.itrs_coordinate.convert_to, AzEl)

    def test_gcrs_to_azel(self):
        self.assertIsInstance(self.gcrs_coordinate.convert_to(AzEl, self.location), AzEl)

    def test_gcrs_to_azel_raises_error_with_no_location(self):
        self.assertRaises(ValueError, self.gcrs_coordinate.convert_to, AzEl)

    def test_wgs84_to_azel(self):
        self.assertIsInstance(self.wgs84_coordinate.convert_to(AzEl, self.location), AzEl)

    def test_wgs84_to_azel_raises_error_with_no_location(self):
        self.assertRaises(ValueError, self.wgs84_coordinate.convert_to, AzEl)
