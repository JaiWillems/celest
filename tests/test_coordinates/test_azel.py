

from celest.coordinates.azel import AzEl
from celest.coordinates.ground_location import GroundLocation
from celest import units as u
from unittest import TestCase
import numpy as np
import os


class TestAzEl(TestCase):

    def setUp(self):
        self.julian = np.random.rand(5)
        self.azimuth = np.random.rand(5)
        self.elevation = np.random.rand(5)
        self.location = GroundLocation(43.6532, -79.3832, 175, u.deg, u.m)
        self.unit = u.deg
        self.valid_azel = AzEl(self.julian, self.azimuth, self.elevation,
                               self.unit, self.location)

    def test_initialization(self):
        self.assertIsInstance(self.valid_azel, AzEl)

    def test_parent_class_is_initialized(self):
        self.assertTrue(np.array_equal(self.julian, self.valid_azel.time.data))
        self.assertEqual(u.jd2000, self.valid_azel.time.unit)

        self.assertTrue(np.array_equal(self.azimuth, self.valid_azel.x.data))
        self.assertEqual(self.unit, self.valid_azel.x.unit)

        self.assertTrue(np.array_equal(self.elevation, self.valid_azel.y.data))
        self.assertEqual(self.unit, self.valid_azel.y.unit)

        self.assertEqual(self.location, self.valid_azel.location)

    def test_error_raised_with_different_length_inputs(self):
        julian = np.random.rand(1)
        azimuth = np.random.rand(2)
        elevation = np.random.rand(3)
        self.assertRaises(ValueError, AzEl, julian, azimuth, elevation,
                          self.unit, self.location)

    def test_error_raised_when_wrong_dimension_input_passed_in(self):
        julian = np.random.rand(1, 1)
        azimuth = np.random.rand(2, 2)
        elevation = np.random.rand(3, 3)
        self.assertRaises(ValueError, AzEl, julian, azimuth, elevation,
                          self.unit, self.location)

    def test_save_to_text_file_saves_file(self):
        file_name = "azel_test_file"
        self.valid_azel.save_text_file(file_name)
        self.assertTrue(os.path.exists(file_name + ".txt"))
        if os.path.exists(file_name + ".txt"):
            os.remove(file_name + ".txt")
