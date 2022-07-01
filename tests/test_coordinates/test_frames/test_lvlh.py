

from celest.coordinates.frames.lvlh import LVLH
from celest.coordinates.ground_location import GroundLocation
from celest import units as u
from unittest import TestCase
import numpy as np
import os


class TestLVLH(TestCase):

    def setUp(self):
        self.julian = np.random.rand(5)
        self.roll = np.random.rand(5)
        self.pitch = np.random.rand(5)
        self.yaw = np.random.rand(5)
        self.unit = u.km
        self.location = GroundLocation(43.6532, -79.3832, 175, u.deg, u.m)
        self.valid_lvlh = LVLH(self.julian, self.roll, self.pitch, self.yaw,
                               self.unit, self.location)

    def test_initialization(self):
        self.assertIsInstance(self.valid_lvlh, LVLH)

    def test_parent_class_is_initialized(self):
        self.assertTrue(np.array_equal(self.julian, self.valid_lvlh.time.data))
        self.assertEqual(u.jd2000, self.valid_lvlh.time.unit)

        self.assertTrue(np.array_equal(self.roll, self.valid_lvlh.x.data))
        self.assertEqual(self.unit, self.valid_lvlh.x.unit)

        self.assertTrue(np.array_equal(self.pitch, self.valid_lvlh.y.data))
        self.assertEqual(self.unit, self.valid_lvlh.y.unit)

        self.assertTrue(np.array_equal(self.yaw, self.valid_lvlh.z.data))
        self.assertEqual(self.unit, self.valid_lvlh.z.unit)

    def test_error_raised_with_different_length_inputs(self):
        julian = np.random.rand(1)
        x = np.random.rand(2)
        y = np.random.rand(3)
        z = np.random.rand(4)
        self.assertRaises(ValueError, LVLH, julian, x, y, z, self.unit,
                          self.location)

    def test_error_raised_when_wrong_dimension_input_passed_in(self):
        julian = np.random.rand(1, 1)
        x = np.random.rand(2, 2)
        y = np.random.rand(3, 3)
        z = np.random.rand(4, 4)
        self.assertRaises(ValueError, LVLH, julian, x, y, z, self.unit,
                          self.location)

    def test_save_to_text_file_saves_file(self):
        file_name = "lvlh_test_file"
        self.valid_lvlh.save_text_file(file_name)
        self.assertTrue(os.path.exists(file_name + ".txt"))
        if os.path.exists(file_name + ".txt"):
            os.remove(file_name + ".txt")
