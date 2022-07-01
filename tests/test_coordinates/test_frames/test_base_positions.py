

from celest.coordinates.frames.base_positions import Position2d, Position3d
from celest.units.quantity import Quantity
from celest import units as u
from unittest import TestCase
import numpy as np


class TestPosition2d(TestCase):

    def setUp(self):
        self.x = np.random.rand(5)
        self.x_unit = u.m
        self.y = np.random.rand(5)
        self.y_unit = u.m
        self.time = np.random.rand(5)
        self.time_unit = u.s
        self.position = Position2d(self.x, self.x_unit, self.y, self.y_unit,
                                   self.time, self.time_unit)

    def test_initialization(self):
        self.assertIsInstance(self.position, Position2d)

    def test_error_raised_with_wrong_dimensioned_arrays(self):
        x = np.random.rand(5, 5)
        y = np.random.rand(5, 5)
        time = np.random.rand(5, 5)
        self.assertRaises(ValueError, Position2d, x, self.x_unit, y,
                          self.y_unit, time, self.time_unit)

    def test_initialization_with_non_arrays_of_different_lengths(self):
        x = np.random.rand(1)
        y = np.random.rand(2)
        time = np.random.rand(3)
        self.assertRaises(ValueError, Position2d, x, self.x_unit, y,
                          self.y_unit, time, self.time_unit)

    def test_x(self):
        self.assertTrue(np.array_equal(self.x, self.position.x.data))
        self.assertEqual(self.x_unit, self.position.x.unit)

    def test_y(self):
        self.assertTrue(np.array_equal(self.y, self.position.y.data))
        self.assertEqual(self.y_unit, self.position.y.unit)

    def test_time(self):
        self.assertTrue(np.array_equal(self.time, self.position.time.data))
        self.assertEqual(self.time_unit, self.position.time.unit)


class TestPosition3d(TestCase):

    def setUp(self):
        self.x = np.random.rand(5)
        self.x_unit = u.m
        self.y = np.random.rand(5)
        self.y_unit = u.m
        self.z = np.random.rand(5)
        self.z_unit = u.m
        self.time = np.random.rand(5)
        self.time_unit = u.s
        self.position = Position3d(self.x, self.x_unit, self.y, self.y_unit,
                                   self.z, self.z_unit, self.time,
                                   self.time_unit)

    def test_initialization(self):
        self.assertIsInstance(self.position, Position3d)

    def test_error_raised_with_wrong_dimensioned_arrays(self):
        x = np.random.rand(5, 5)
        y = np.random.rand(5, 5)
        z = np.random.rand(5, 5)
        time = np.random.rand(5, 5)
        self.assertRaises(ValueError, Position3d, x, self.x_unit, y,
                          self.y_unit, z, self.z_unit, time, self.time_unit)

    def test_initialization_with_non_arrays_of_different_lengths(self):
        x = np.random.rand(1)
        y = np.random.rand(2)
        z = np.random.rand(3)
        time = np.random.rand(3)
        self.assertRaises(ValueError, Position3d, x, self.x_unit, y,
                          self.y_unit, z, self.z_unit, time, self.time_unit)

    def test_x(self):
        self.assertTrue(np.array_equal(self.x, self.position.x.data))
        self.assertEqual(self.x_unit, self.position.x.unit)

    def test_y(self):
        self.assertTrue(np.array_equal(self.y, self.position.y.data))
        self.assertEqual(self.y_unit, self.position.y.unit)

    def test_z(self):
        self.assertTrue(np.array_equal(self.z, self.position.z.data))
        self.assertEqual(self.z_unit, self.position.z.unit)

    def test_time(self):
        self.assertTrue(np.array_equal(self.time, self.position.time.data))
        self.assertEqual(self.time_unit, self.position.time.unit)
