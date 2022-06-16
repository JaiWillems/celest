

from celest.units.utils import setup_unit
from celest.units.core import Unit
from unittest import TestCase


class TestUtils(TestCase):

    def setUp(self):
        self.namespace = globals()

    def test_initializing_elementary_unit(self):
        setup_unit("m", "meter", self.namespace)
        self.assertIsInstance(m, Unit)
        self.assertEqual(m.scale, 1.0)
        self.assertListEqual(m.bases, [m])
        self.assertListEqual(m.powers, [1])

    def test_initializing_non_elementary_unit(self):
        setup_unit("m", "meter", self.namespace)
        setup_unit("km", "kilometer", self.namespace, 1000 * m)
        self.assertIsInstance(km, Unit)
        self.assertEqual(km.scale, 1000.0)
        self.assertListEqual(km.bases, [m])
        self.assertListEqual(km.powers, [1])
