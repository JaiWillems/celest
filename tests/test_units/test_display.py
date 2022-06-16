

from celest.units.core import NamedUnit
from unittest import TestCase


class TestToString(TestCase):

    def test_namedunit_to_string(self):
        self.assertEqual(str(NamedUnit("m", 'meter')), "m")

    def test_to_string_with_single_unit_positive_power(self):
        self.assertEqual(str(NamedUnit("m", 'meter') ** 2), "m2")

    def test_to_string_with_single_unit_positive_power(self):
        self.assertEqual(str(NamedUnit("m", 'meter') ** 3), "m3")

    def test_to_string_with_single_unit_negative_power(self):
        self.assertEqual(str(NamedUnit("m", 'meter') ** -2), "1 / m2")

    def test_to_string_with_single_unit_negative_unity_power(self):
        self.assertEqual(str(NamedUnit("m", 'meter') ** -1), "1 / m")

    def test_to_string_with_double_numerator_units(self):
        composite_unit = NamedUnit("m", 'meter') * NamedUnit("s", 'second') ** 2
        self.assertEqual(str(composite_unit), "m s2")

    def test_to_string_with_double_denominator_units(self):
        composite_unit = NamedUnit("m", 'meter') ** -1 * NamedUnit("s", 'second') ** -2
        self.assertEqual(str(composite_unit), "1 / m s2")

    def test_to_string_with_ratio_units(self):
        composite_unit = NamedUnit("m", 'meter') / NamedUnit("s", 'second') ** 2
        self.assertEqual(str(composite_unit), "m / s2")

    def test_to_string_with_non_unity_scale(self):
        composite_unit = NamedUnit("m", 'meter') / NamedUnit("s", 'second') ** 2
        composite_unit = 5 * composite_unit
        self.assertEqual(str(composite_unit), "5 m / s2")
