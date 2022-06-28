

from celest.coordinates.wgs84 import WGS84
from celest import units as u
from unittest import TestCase
import numpy as np
import os


class TestWGS84(TestCase):

    def setUp(self):
        self.julian = np.random.rand(5)
        self.latitude = np.random.rand(5)
        self.longitude = np.random.rand(5)
        self.height = np.random.rand(5)
        self.angular_unit = u.deg
        self.length_unit = u.km
        self.valid_wgs84 = WGS84(self.julian, self.latitude, self.longitude,
                                 self.height, self.angular_unit,
                                 self.length_unit)

    def test_initialization(self):
        self.assertIsInstance(self.valid_wgs84, WGS84)

    def test_parent_class_is_initialized(self):
        self.assertTrue(np.array_equal(self.julian, self.valid_wgs84.time.data))
        self.assertEqual(u.jd2000, self.valid_wgs84.time.unit)

        self.assertTrue(np.array_equal(self.latitude,
                                       self.valid_wgs84.latitude.data))
        self.assertEqual(self.angular_unit, self.valid_wgs84.latitude.unit)

        self.assertTrue(np.array_equal(self.longitude,
                                       self.valid_wgs84.longitude.data))
        self.assertEqual(self.angular_unit, self.valid_wgs84.longitude.unit)

        self.assertTrue(np.array_equal(self.height,
                                       self.valid_wgs84.height.data))
        self.assertEqual(self.length_unit, self.valid_wgs84.height.unit)

    def test_error_raised_with_different_length_inputs(self):
        julian = np.random.rand(1)
        latitude = np.random.rand(2)
        longitude = np.random.rand(3)
        height = np.random.rand(4)
        self.assertRaises(ValueError, WGS84, julian, latitude, longitude,
                          height, self.angular_unit, self.length_unit)

    def test_error_raised_when_wrong_dimension_input_passed_in(self):
        julian = np.random.rand(1, 1)
        latitude = np.random.rand(2, 2)
        longitude = np.random.rand(3, 3)
        height = np.random.rand(4, 4)
        self.assertRaises(ValueError, WGS84, julian, latitude, longitude,
                          height, self.angular_unit, self.length_unit)

    def test_save_to_text_file_saves_file(self):
        file_name = "wgs84_test_file"
        self.valid_wgs84.save_text_file(file_name)
        self.assertTrue(os.path.exists(file_name + ".txt"))
        if os.path.exists(file_name + ".txt"):
            os.remove(file_name + ".txt")
