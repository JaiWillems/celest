

from celest.units.core import NamedUnit
from celest import units as u
from unittest import TestCase


class TestToString(TestCase):

    def test_namedunit_to_string(self):
        self.assertEqual(str(NamedUnit("m", 'meter')), "m")

    def test_to_string_with_single_unit_positive_power(self):
        self.assertEqual(str(NamedUnit("m", 'meter') ** 2), "m2")

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


class TestGetUnitDimensionalityString(TestCase):

    def test_get_dimension_string_from_elementary_length_unit(self):
        self.assertEqual("L", u.m.dimension)

    def test_get_dimension_string_from_elementary_time_unit(self):
        self.assertEqual("T", u.s.dimension)

    def test_get_dimension_string_from_elementary_angle_unit(self):
        self.assertEqual("A", u.deg.dimension)

    def test_get_dimension_string_with_no_denominator_unit(self):
        self.assertEqual("L3", (u.m ** 3).dimension)

    def test_get_dimension_string_with_no_numerator_unit(self):
        self.assertEqual("1 / L3", (u.m ** -3).dimension)

    def test_get_dimension_string_from_non_repeating_compound_unit(self):
        self.assertEqual("L A / T3", (u.m * u.deg / u.s ** 3).dimension)

    def test_get_dimension_string_with_repeating_numerator(self):
        self.assertEqual("L2 / T3", (u.m * u.km / u.s ** 3).dimension)

    def test_get_dimension_string_with_repeating_denominator(self):
        self.assertEqual("L / T4", (u.m / (u.s ** 3 * u.min)).dimension)

    def test_get_dimension_string_with_cancelling_dimensions(self):
        self.assertEqual("L / T2", (u.m * u.min / u.s ** 3).dimension)
