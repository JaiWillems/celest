

from celest.coordinates.frames.attitude import Attitude
from celest.coordinates.ground_location import GroundLocation
from celest import units as u
from unittest import TestCase
import numpy as np
import os


class TestAttitude(TestCase):

    def setUp(self):
        self.julian = np.random.rand(5)
        self.roll = np.random.rand(5)
        self.pitch = np.random.rand(5)
        self.yaw = np.random.rand(5)
        self.unit = u.deg
        self.location = GroundLocation(52.1579, -106.6702, 0.482, u.deg, u.km)
        self.valid_attitude = Attitude(
            self.julian,
            self.roll,
            self.pitch,
            self.yaw,
            self.unit,
            self.location
        )

    def test_initialization(self):
        self.assertIsInstance(self.valid_attitude, Attitude)

    def test_parent_class_is_initialized(self):
        self.assertTrue(np.array_equal(self.julian, self.valid_attitude.time.data))
        self.assertEqual(u.jd2000, self.valid_attitude.time.unit)

        self.assertTrue(np.array_equal(self.roll, self.valid_attitude.roll.data))
        self.assertEqual(self.unit, self.valid_attitude.roll.unit)

        self.assertTrue(np.array_equal(self.pitch, self.valid_attitude.pitch.data))
        self.assertEqual(self.unit, self.valid_attitude.pitch.unit)

        self.assertTrue(np.array_equal(self.yaw, self.valid_attitude.yaw.data))
        self.assertEqual(self.unit, self.valid_attitude.yaw.unit)

    def test_error_raised_with_different_length_inputs(self):
        julian = np.random.rand(1)
        x = np.random.rand(2)
        y = np.random.rand(3)
        z = np.random.rand(4)
        self.assertRaises(
            ValueError,
            Attitude,
            julian,
            x,
            y,
            z,
            self.unit,
            self.location
        )

    def test_error_raised_when_wrong_dimension_input_passed_in(self):
        julian = np.random.rand(1, 1)
        x = np.random.rand(2, 2)
        y = np.random.rand(3, 3)
        z = np.random.rand(4, 4)
        self.assertRaises(
            ValueError,
            Attitude,
            julian,
            x,
            y,
            z,
            self.unit,
            self.location
        )

    def test_save_to_text_file_saves_file(self):
        file_name = "attitude_test_file"
        self.valid_attitude.save_text_file(file_name)
        self.assertTrue(os.path.exists(file_name + ".txt"))
        if os.path.exists(file_name + ".txt"):
            os.remove(file_name + ".txt")
