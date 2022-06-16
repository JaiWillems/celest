

from celest.units.core import BaseUnit, NamedUnit, Unit, CompoundUnit
from unittest import TestCase


class TestBaseUnit(TestCase):

    def test_scale(self):
        self.assertEqual(BaseUnit().scale, 1.0)

    def test_bases(self):
        base_unit = BaseUnit()
        self.assertListEqual(base_unit.bases, [base_unit])

    def test_powers(self):
        self.assertListEqual(BaseUnit().powers, [1])

    def test_valid_power_returns_composite_unit(self):
        self.assertIsInstance(BaseUnit() ** 2, CompoundUnit)

    def test_valid_power_updates_power_list(self):
        composite_unit = BaseUnit() ** 2
        self.assertListEqual(composite_unit.powers, [2])

    def test_zero_power_returns_unity(self):
        composite_unit = BaseUnit() ** 0

        self.assertEqual(composite_unit.scale, 1.0)
        self.assertListEqual(composite_unit.bases, [])
        self.assertListEqual(composite_unit.powers, [])

    def test_error_raised_for_non_number_powers(self):
        self.assertRaises(TypeError, lambda: BaseUnit() ** "bad idea")

    def test_dividing_by_BaseUnit(self):
        base_unit_1, base_unit_2 = BaseUnit(), BaseUnit()
        composite_unit = base_unit_1.__truediv__(base_unit_2)
        self.assertIsInstance(composite_unit, CompoundUnit)
        self.assertEqual(composite_unit.scale, 1.0)
        self.assertListEqual(composite_unit.bases, [base_unit_1, base_unit_2])
        self.assertListEqual(composite_unit.powers, [1, -1])

    def test_right_dividing_by_BaseUnit(self):
        base_unit_1, base_unit_2 = BaseUnit(), BaseUnit()
        composite_unit = base_unit_1.__rtruediv__(base_unit_2)
        self.assertIsInstance(composite_unit, CompoundUnit)
        self.assertEqual(composite_unit.scale, 1.0)
        self.assertListEqual(composite_unit.bases, [base_unit_2, base_unit_1])
        self.assertListEqual(composite_unit.powers, [1, -1])

    def test_multiply_returns_composite(self):
        base_unit_1, base_unit_2 = BaseUnit(), BaseUnit()
        composite_unit = base_unit_1.__mul__(base_unit_2)
        self.assertIsInstance(composite_unit, CompoundUnit)
        self.assertEqual(composite_unit.scale, 1.0)
        self.assertListEqual(composite_unit.bases, [base_unit_1, base_unit_2])
        self.assertListEqual(composite_unit.powers, [1, 1])

    def test_multiply_by_number_returns_composite(self):
        named_unit = NamedUnit("m", "meter")
        composite_unit = named_unit.__mul__(2.5)
        self.assertIsInstance(composite_unit, CompoundUnit)
        self.assertEqual(composite_unit.scale, 2.5)
        self.assertListEqual(composite_unit.bases, [named_unit])
        self.assertListEqual(composite_unit.powers, [1])

    def test_right_multiply_returns_composite(self):
        base_unit_1, base_unit_2 = BaseUnit(), BaseUnit()
        composite_unit = base_unit_1.__rmul__(base_unit_2)
        self.assertIsInstance(composite_unit, CompoundUnit)
        self.assertEqual(composite_unit.scale, 1.0)
        self.assertListEqual(composite_unit.bases, [base_unit_2, base_unit_1])
        self.assertListEqual(composite_unit.powers, [1, 1])

    def test_right_multiply_by_number_returns_composite(self):
        named_unit = NamedUnit("m", "meter")
        composite_unit = named_unit.__rmul__(2.5)
        self.assertIsInstance(composite_unit, CompoundUnit)
        self.assertEqual(composite_unit.scale, 2.5)
        self.assertListEqual(composite_unit.bases, [named_unit])
        self.assertListEqual(composite_unit.powers, [1])


class TestNamedUnit(TestCase):

    def test_str_returns_string(self):
        self.assertEqual(str(NamedUnit("m", "meter")), "m")

    def test_repr_returns_string(self):
        self.assertEqual(repr(NamedUnit("m", "meter")), "NamedUnit(\"m\")")

    def test_short_name(self):
        self.assertEqual(NamedUnit("m", "meter").short_name, "m")

    def test_long_name(self):
        self.assertEqual(NamedUnit("m", "meter").long_name, "meter")

    def test_insert_into_namespace(self):
        namespace = globals()
        named_unit = NamedUnit("m", "meter", namespace=namespace)
        self.assertEqual(named_unit, m)

    def test_inserting_multiple_instances_of_same_unit(self):
        namespace = globals()
        NamedUnit("m", "meter", namespace=namespace)
        m2 = m * m
        self.assertIsInstance(m2, CompoundUnit)
        self.assertEqual(m2.scale, 1.0)
        self.assertListEqual(m2.bases, [m])
        self.assertListEqual(m2.powers, [2])


class TestUnit(TestCase):

    def test_repr_returns_Unit_types_for_elementary_unit(self):
        self.assertEqual(repr(Unit("m", "meter")), "Unit(\"m\")")

    def test_repr_returns_Unit_types_for_non_elementary_unit(self):
        self.assertEqual(repr(Unit("m", "meter")), "Unit(\"m\")")

    def test_initialization_with_elementary_units(self):

        m = Unit("m", "meter")
        self.assertIsInstance(m, Unit)
        self.assertEqual(m.scale, 1.0)
        self.assertListEqual(m.bases, [m])
        self.assertListEqual(m.powers, [1])

    def test_initialization_with_non_elementary_units(self):

        m = Unit("m", "meter")
        km = Unit("km", "kilometer", base_units=1000*m)
        self.assertIsInstance(km, Unit)
        self.assertEqual(km.scale, 1000.0)
        self.assertListEqual(km.bases, [m])
        self.assertListEqual(km.powers, [1])

    def test_multiplication_retains_units(self):

        m = Unit("m", "meter")
        s = Unit("s", "second")
        ms = m * s
        self.assertIsInstance(ms, CompoundUnit)
        self.assertEqual(ms.scale, 1.0)
        self.assertListEqual(ms.bases, [m, s])
        self.assertListEqual(ms.powers, [1, 1])

    def test_division_retains_units(self):

        m = Unit("m", "meter")
        s = Unit("s", "second")
        m_per_s = m / s
        self.assertIsInstance(m_per_s, CompoundUnit)
        self.assertEqual(m_per_s.scale, 1.0)
        self.assertListEqual(m_per_s.bases, [m, s])
        self.assertListEqual(m_per_s.powers, [1, -1])

    def test_exponential_retains_units(self):
        m = Unit("m", "meter")
        m2 = m ** 2
        self.assertIsInstance(m2, CompoundUnit)
        self.assertEqual(m2.scale, 1.0)
        self.assertListEqual(m2.bases, [m])
        self.assertListEqual(m2.powers, [2])

    def test_division_does_not_simplify(self):

        m = Unit("m", "meter")
        km = Unit("km", "kilometer", base_units=1000 * m)
        km_per_m = km / m
        self.assertIsInstance(km_per_m, CompoundUnit)
        self.assertEqual(km_per_m.scale, 1.0)
        self.assertListEqual(km_per_m.bases, [km, m])
        self.assertListEqual(km_per_m.powers, [1, -1])


class TestCompoundUnit(TestCase):

    def test_str_returns_string(self):
        composite_unit = NamedUnit("m", "meter") ** 1 / NamedUnit("s", "second") ** 3
        self.assertEqual(str(composite_unit), "m / s3")

    def test_repr_returns_string(self):
        composite_unit = NamedUnit("m", "meter") ** 1 / NamedUnit("s", "second") ** 3
        self.assertEqual(repr(composite_unit), "CompoundUnit(\"m / s3\")")

    def test_scale(self):
        self.assertEqual(CompoundUnit(1.0, [BaseUnit(), BaseUnit()],
                                      [1, 1]).scale, 1.0)

    def test_bases(self):
        bases = [BaseUnit(), BaseUnit()]
        self.assertListEqual(CompoundUnit(1.0, bases, [1, 1]).bases, bases)

    def test_powers(self):
        self.assertListEqual(CompoundUnit(1.0, [BaseUnit(), BaseUnit()],
                                          [1, 1]).powers, [1, 1])

    def test_expand_and_collect_composite_base_no_duplicates(self):
        base_unit_1 = BaseUnit()
        base_unit_2 = BaseUnit()
        base_unit_3 = BaseUnit()
        composite_unit = (base_unit_1 * base_unit_2) / base_unit_3

        self.assertEqual(composite_unit.scale, 1.0)
        self.assertListEqual(composite_unit.bases, [base_unit_1, base_unit_2,
                                                    base_unit_3])
        self.assertListEqual(composite_unit.powers, [1, 1, -1])

    def test_expand_and_collect_composite_base_with_duplicates(self):
        base_unit_1 = BaseUnit()
        base_unit_2 = BaseUnit()
        composite_unit = (base_unit_1 * base_unit_1) / base_unit_2

        self.assertEqual(composite_unit.scale, 1.0)
        self.assertListEqual(composite_unit.bases, [base_unit_1, base_unit_2])
        self.assertListEqual(composite_unit.powers, [2, -1])

    def test_expand_and_collect_base_composite_no_duplicates(self):
        base_unit_1 = BaseUnit()
        base_unit_2 = BaseUnit()
        base_unit_3 = BaseUnit()
        composite_unit = base_unit_1 * (base_unit_2 / base_unit_3)

        self.assertEqual(composite_unit.scale, 1.0)
        self.assertListEqual(composite_unit.bases, [base_unit_1, base_unit_2,
                                                    base_unit_3])
        self.assertListEqual(composite_unit.powers, [1, 1, -1])

    def test_expand_and_collect_base_composite_with_duplicates(self):
        base_unit_1 = BaseUnit()
        base_unit_2 = BaseUnit()
        composite_unit = base_unit_1 * (base_unit_1 / base_unit_2)

        self.assertEqual(composite_unit.scale, 1.0)
        self.assertListEqual(composite_unit.bases, [base_unit_1, base_unit_2])
        self.assertListEqual(composite_unit.powers, [2, -1])

    def test_expand_and_collect_composite_composite_no_duplicates(self):
        base_unit_1 = BaseUnit()
        base_unit_2 = BaseUnit()
        base_unit_3 = BaseUnit()
        composite_unit = (base_unit_1 * base_unit_2)
        composite_unit = composite_unit / (base_unit_3 ** 2)

        self.assertEqual(composite_unit.scale, 1.0)
        self.assertListEqual(composite_unit.bases, [base_unit_1, base_unit_2,
                                                    base_unit_3])
        self.assertListEqual(composite_unit.powers, [1, 1, -2])

    def test_expand_and_collect_composite_composite_with_duplicates(self):
        base_unit_1 = BaseUnit()
        base_unit_2 = BaseUnit()
        base_unit_3 = BaseUnit()
        composite_unit = (base_unit_1 * base_unit_2)
        composite_unit = composite_unit * (base_unit_2 / (base_unit_3 ** 3))

        self.assertEqual(composite_unit.scale, 1.0)
        self.assertListEqual(composite_unit.bases, [base_unit_1, base_unit_2,
                                                    base_unit_3])
        self.assertListEqual(composite_unit.powers, [1, 2, -3])

    def test_expand_and_collect_with_zero_power(self):
        base_unit_1 = BaseUnit()
        base_unit_2 = BaseUnit()
        composite_unit = (base_unit_1 * base_unit_2) / base_unit_1

        self.assertEqual(composite_unit.scale, 1.0)
        self.assertListEqual(composite_unit.bases, [base_unit_2])
        self.assertListEqual(composite_unit.powers, [1])

    def test_expand_and_collect_with_zero_power_no_other_unit(self):
        base_unit_1 = BaseUnit()
        composite_unit = base_unit_1 / base_unit_1

        self.assertEqual(composite_unit.scale, 1.0)
        self.assertListEqual(composite_unit.bases, [])
        self.assertListEqual(composite_unit.powers, [])

    def test_multiply_by_unity_composite_returns_other(self):
        base_unit_1 = BaseUnit()
        base_unit_2 = BaseUnit()
        composite_unit = (base_unit_1 / base_unit_1) * base_unit_2

        self.assertEqual(composite_unit.scale, 1.0)
        self.assertListEqual(composite_unit.bases, [base_unit_2])
        self.assertListEqual(composite_unit.powers, [1])

    def test_expand_and_collect_with_unity_length_with_composite(self):
        named_unit_1 = NamedUnit("m", 'meter')
        named_unit_2 = NamedUnit("s", 'second')
        composite_unit = 5 * (named_unit_1 / named_unit_2 ** 2)

        self.assertEqual(composite_unit.scale, 5.0)
        self.assertListEqual(composite_unit.bases, [named_unit_1, named_unit_2])
        self.assertListEqual(composite_unit.powers, [1, -2])
