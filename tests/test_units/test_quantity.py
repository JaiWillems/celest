

from celest.units.quantity import Quantity
from celest import units as u
from unittest import TestCase


class TestQuantity(TestCase):

    def setUp(self):
        self.data = 5

        self.simple_unit = u.m
        self.simple_quantity = Quantity(self.data, self.simple_unit)
        self.simple_quantity_2 = Quantity(2 * self.data, self.simple_unit)

        self.compound_unit = u.m / u.s
        self.compound_quantity = Quantity(self.data, self.compound_unit)

    def test_initialization(self):
        self.assertEqual(self.data, self.simple_quantity.data)
        self.assertEqual(self.simple_unit, self.simple_quantity.unit)

    def test_str(self):
        self.assertEqual(f"{str(self.data)} {str(self.simple_unit)}",
                         str(self.simple_quantity))

    def test_repr(self):
        self.assertEqual(f"Quantity({self.data}, {repr(self.simple_unit)})",
                         repr(self.simple_quantity))

    def test_negation(self):
        result = -self.simple_quantity
        self.assertEqual(-self.data, result.data)
        self.assertEqual(self.simple_unit, result.unit)

    def test_add_scalar(self):
        result = self.simple_quantity.__add__(self.data)
        self.assertEqual(self.data + self.data, result.data)
        self.assertEqual(self.simple_unit, result.unit)

    def test_add_quantity_of_same_dimensions(self):
        result = self.simple_quantity.__add__(self.simple_quantity)
        self.assertEqual(2 * self.data, result.data)
        self.assertEqual(self.simple_unit, result.unit)

    def test_add_quantity_of_different_dimensions_raises_error(self):
        self.assertRaises(ArithmeticError, self.simple_quantity.__add__,
                          self.compound_quantity)

    def test_radd_scalar(self):
        result = self.simple_quantity.__radd__(self.data)
        self.assertEqual(self.data + self.data, result.data)
        self.assertEqual(self.simple_unit, result.unit)

    def test_radd_quantity_of_same_dimensions(self):
        result = self.simple_quantity.__radd__(self.simple_quantity)
        self.assertEqual(2 * self.data, result.data)
        self.assertEqual(self.simple_unit, result.unit)

    def test_radd_quantity_of_different_dimensions_raises_error(self):
        self.assertRaises(ArithmeticError, self.simple_quantity.__radd__,
                          self.compound_quantity)

    def test_sub_scalar(self):
        result = self.simple_quantity - 7
        self.assertEqual(self.data - 7, result.data)
        self.assertEqual(self.simple_unit, result.unit)

    def test_sub_quantity_of_same_dimensions(self):
        result = self.simple_quantity.__sub__(self.simple_quantity)
        self.assertEqual(0, result.data)
        self.assertEqual(self.simple_unit, result.unit)

    def test_sub_quantity_of_different_dimensions_raises_error(self):
        self.assertRaises(ArithmeticError, self.simple_quantity.__sub__,
                          self.compound_quantity)

    def test_rsub_scalar(self):
        result = 7 - self.simple_quantity
        self.assertEqual(7 - self.data, result.data)
        self.assertEqual(self.simple_unit, result.unit)

    def test_rsub_quantity_of_same_dimensions(self):
        result = self.simple_quantity.__rsub__(self.simple_quantity)
        self.assertEqual(0, result.data)
        self.assertEqual(self.simple_unit, result.unit)

    def test_rsub_quantity_of_different_dimensions_raises_error(self):
        self.assertRaises(ArithmeticError, self.simple_quantity.__rsub__,
                          self.compound_quantity)

    def test_mul_scalar(self):
        result = self.simple_quantity.__mul__(self.data)
        self.assertEqual(self.data * self.data, result.data)
        self.assertEqual(self.simple_unit, result.unit)

    def test_mul_quantity_of_same_dimensions(self):
        result = self.simple_quantity.__mul__(self.simple_quantity)
        self.assertEqual(self.data * self.data, result.data)
        self.assertEqual(repr(self.simple_unit ** 2), repr(result.unit))

    def test_mul_quantity_of_different_dimensions(self):
        result = self.simple_quantity.__mul__(self.compound_quantity)
        self.assertEqual(self.data * self.data, result.data)
        self.assertEqual(repr(self.simple_unit * self.compound_unit),
                         repr(result.unit))

    def test_rmul_scalar(self):
        result = self.simple_quantity.__rmul__(self.data)
        self.assertEqual(self.data * self.data, result.data)
        self.assertEqual(self.simple_unit, result.unit)

    def test_rmul_quantity_of_same_dimensions(self):
        result = self.simple_quantity.__rmul__(self.simple_quantity)
        self.assertEqual(self.data * self.data, result.data)
        self.assertEqual(repr(self.simple_unit ** 2), repr(result.unit))

    def test_rmul_quantity_of_different_dimensions(self):
        result = self.simple_quantity.__rmul__(self.compound_quantity)
        self.assertEqual(self.data * self.data, result.data)
        self.assertEqual(repr(self.simple_unit * self.compound_unit),
                         repr(result.unit))

    def test_truediv_scaler(self):
        result = self.simple_quantity.__truediv__(2 * self.data)
        self.assertEqual(0.5, result.data)
        self.assertEqual(self.simple_unit, result.unit)

    def test_truediv_quantity_of_same_dimensions(self):
        result = self.simple_quantity.__truediv__(self.simple_quantity)
        self.assertEqual(1, result.data)
        self.assertEqual(repr(self.simple_unit / self.simple_unit),
                         repr(result.unit))

    def test_truediv_quantity_of_different_dimensions(self):
        result = self.simple_quantity.__truediv__(self.compound_quantity)
        self.assertEqual(1, result.data)
        self.assertEqual(repr(self.simple_unit / self.compound_unit),
                         repr(result.unit))

    def test_rtruediv_scaler(self):
        result = self.simple_quantity.__rtruediv__(2 * self.data)
        self.assertEqual(2, result.data)
        self.assertEqual(repr(1 / self.simple_unit), repr(result.unit))

    def test_rtruediv_quantity_of_same_dimensions(self):
        result = self.simple_quantity.__rtruediv__(self.simple_quantity)
        self.assertEqual(1, result.data)
        self.assertEqual(repr(self.simple_unit / self.simple_unit),
                         repr(result.unit))

    def test_rtruediv_quantity_of_different_dimensions(self):
        result = self.simple_quantity.__rtruediv__(self.compound_quantity)
        self.assertEqual(1, result.data)
        self.assertEqual(repr(self.compound_unit / self.simple_unit),
                         repr(result.unit))

    def test_eq_when_equal_with_same_dimensions(self):
        self.assertTrue(self.simple_quantity == self.simple_quantity)

    def test_eq_when_not_equal_with_same_dimensions(self):
        self.assertFalse(self.simple_quantity == (self.simple_quantity + 1))

    def test_eq_with_different_dimensions_raises_arithmetic_error(self):
        self.assertRaises(ArithmeticError, self.simple_quantity.__eq__,
                          self.compound_quantity)

    def test_neq_when_not_equal_with_same_dimensions(self):
        self.assertTrue(self.simple_quantity != (self.simple_quantity + 1))

    def test_neq_when_equal_with_same_dimensions(self):
        self.assertFalse(self.simple_quantity != self.simple_quantity)

    def test_neq_with_different_dimensions_raises_arithmetic_error(self):
        self.assertRaises(ArithmeticError, self.simple_quantity.__ne__,
                          self.compound_quantity)

    def test_lt_when_less_than_with_same_dimensions(self):
        self.assertTrue(self.simple_quantity < (self.simple_quantity + 1))

    def test_lt_when_greater_than_with_same_dimensions(self):
        self.assertFalse((self.simple_quantity + 1) < self.simple_quantity)

    def test_lt_with_different_dimensions_raises_arithmetic_error(self):
        self.assertRaises(ArithmeticError, self.simple_quantity.__lt__,
                          self.compound_quantity)

    def test_le_when_less_than_with_same_dimensions(self):
        self.assertTrue(self.simple_quantity <= (self.simple_quantity + 1))

    def test_le_when_equal_with_same_dimensions(self):
        self.assertTrue(self.simple_quantity <= self.simple_quantity)

    def test_le_when_greater_than_with_same_dimensions(self):
        self.assertFalse((self.simple_quantity + 1) <= self.simple_quantity)

    def test_le_with_different_dimensions_raises_arithmetic_error(self):
        self.assertRaises(ArithmeticError, self.simple_quantity.__le__,
                          self.compound_quantity)

    def test_gt_when_greater_with_same_dimensions(self):
        self.assertTrue((self.simple_quantity + 1) > self.simple_quantity)

    def test_gt_when_less_with_same_dimensions(self):
        self.assertFalse(self.simple_quantity > (self.simple_quantity + 1))

    def test_gt_with_different_dimensions_raises_arithmetic_error(self):
        self.assertRaises(ArithmeticError, self.simple_quantity.__gt__,
                          self.compound_quantity)

    def test_ge_when_greater_than_with_same_dimensions(self):
        self.assertTrue((self.simple_quantity + 1) >= self.simple_quantity)

    def test_ge_when_equal_with_same_dimensions(self):
        self.assertTrue(self.simple_quantity >= self.simple_quantity)

    def test_ge_when_less_than_with_same_dimensions(self):
        self.assertFalse(self.simple_quantity >= (self.simple_quantity + 1))

    def test_ge_with_different_dimensions_raises_arithmetic_error(self):
        self.assertRaises(ArithmeticError, self.simple_quantity.__ge__,
                          self.compound_quantity)

    def test_to_with_same_dimension(self):
        self.assertEqual(0.005, self.simple_quantity.to(u.km))

    def test_to_with_different_dimension_raises_value_error(self):
        self.assertRaises(ValueError, self.simple_quantity.to, u.s)

    def test_convert_to_same_dimension(self):
        new_quantity = self.simple_quantity._convert_to(u.km)
        self.assertEqual(0.005, new_quantity.data)
        self.assertEqual(u.km, new_quantity.unit)

    def test_convert_to_with_compound_unit(self):
        new_compound_unit = u.km / u.hr
        new_quantity = self.compound_quantity._convert_to(new_compound_unit)
        self.assertEqual(18, new_quantity.data)
        self.assertEqual(new_compound_unit, new_quantity.unit)

    def test_convert_to_different_dimension_raises_ValueError(self):
        self.assertRaises(ValueError, self.simple_quantity._convert_to, u.s)

    def test_convert_to_with_compound_unit_with_different_dimensions(self):
        self.assertRaises(ValueError, self.compound_quantity._convert_to, u.s)
