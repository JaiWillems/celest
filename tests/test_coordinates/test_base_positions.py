

from celest.coordinates.base_positions import Position2d, Position3d
from celest.units.quantity import Quantity
from celest import units as u
from unittest import TestCase


class TestPosition2d(TestCase):

    def setUp(self):
        self.x = Quantity(1, u.m)
        self.y = Quantity(2, u.m)
        self.time = Quantity(3, u.s)
        self.position = Position2d(self.x, self.y, self.time)

    def test_initialization(self):
        self.assertIsInstance(self.position, Position2d)

    def test_initialization_with_non_quantity_parameters(self):
        self.assertRaises(ValueError, Position2d, 1, "2", 3)

    def test_x(self):
        self.assertEqual(self.x, self.position.x)

    def test_y(self):
        self.assertEqual(self.y, self.position.y)

    def test_time(self):
        self.assertEqual(self.time, self.position.time)


class TestPosition3d(TestCase):

    def setUp(self):
        self.x = Quantity(1, u.m)
        self.y = Quantity(2, u.m)
        self.z = Quantity(3, u.m)
        self.time = Quantity(4, u.s)
        self.position = Position3d(self.x, self.y, self.z, self.time)

    def test_initialization(self):
        self.assertIsInstance(self.position, Position3d)

    def test_initialization_with_non_quantity_parameters(self):
        self.assertRaises(ValueError, Position3d, 1, "2", 3, "4")

    def test_x(self):
        self.assertEqual(self.x, self.position.x)

    def test_y(self):
        self.assertEqual(self.y, self.position.y)

    def test_z(self):
        self.assertEqual(self.z, self.position.z)

    def test_time(self):
        self.assertEqual(self.time, self.position.time)
