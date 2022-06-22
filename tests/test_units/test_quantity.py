

from celest.units.quantity import Quantity
from celest import units as u
from unittest import TestCase


class TestQuantity(TestCase):

    def setUp(self):
        self.simple_unit = u.m
        self.simple_quantity = Quantity(5, self.simple_unit)

        self.compound_unit = u.m / u.s
        self.compound_quantity = Quantity(5, self.compound_unit)

    def test_initialization(self):
        self.assertEqual(5, self.simple_quantity.data)
        self.assertEqual(self.simple_unit, self.simple_quantity.unit)

    def test_str(self):
        self.assertEqual("5 m", str(self.simple_quantity))

    def test_repr(self):
        self.assertEqual("Quantity(5, Unit(\"m\"))", repr(self.simple_quantity))

    def test_get_unit(self):
        self.assertEqual(self.simple_unit, self.simple_quantity.get_unit())

    def test_convert_to_same_dimension(self):
        self.simple_quantity.convert_to(u.km)
        self.assertEqual(0.005, self.simple_quantity.data)
        self.assertEqual(u.km, self.simple_quantity.unit)

    def test_convert_to_with_compound_unit(self):
        new_compound_unit = u.km / u.hr
        self.compound_quantity.convert_to(new_compound_unit)
        self.assertEqual(18, self.compound_quantity.data)
        self.assertEqual(new_compound_unit, self.compound_quantity.unit)

    def test_convert_to_different_dimension_raises_ValueError(self):
        self.assertRaises(ValueError, self.simple_quantity.convert_to, u.s)

    def test_convert_to_with_compound_unit_with_different_dimensions(self):
        self.assertRaises(ValueError, self.compound_quantity.convert_to, u.s)
