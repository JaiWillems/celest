

from unittest import TestCase
import celest.units as u


class TestUnitInitialization(TestCase):

    def test_meter(self):
        self.assertIsInstance(u.m, u.core.Unit)
        self.assertEqual(u.m.scale, 1.0)
        self.assertEqual(u.m.bases, [u.m])
        self.assertEqual(u.m.powers, [1])

    def test_millimeter(self):
        self.assertIsInstance(u.mm, u.core.Unit)
        self.assertEqual(u.mm.scale, 0.001)
        self.assertEqual(u.mm.bases, [u.m])
        self.assertEqual(u.mm.powers, [1])

    def test_centimeter(self):
        self.assertIsInstance(u.cm, u.core.Unit)
        self.assertEqual(u.cm.scale, 0.01)
        self.assertEqual(u.cm.bases, [u.m])
        self.assertEqual(u.cm.powers, [1])

    def test_kilometer(self):
        self.assertIsInstance(u.km, u.core.Unit)
        self.assertEqual(u.km.scale, 1000.0)
        self.assertEqual(u.km.bases, [u.m])
        self.assertEqual(u.km.powers, [1])

    def test_inch(self):
        self.assertIsInstance(u.inch, u.core.Unit)
        self.assertEqual(u.inch.scale, 0.0254)
        self.assertEqual(u.inch.bases, [u.m])
        self.assertEqual(u.inch.powers, [1])

    def test_feet(self):
        self.assertIsInstance(u.ft, u.core.Unit)
        self.assertEqual(u.ft.scale, 0.305)
        self.assertEqual(u.ft.bases, [u.m])
        self.assertEqual(u.ft.powers, [1])

    def test_yard(self):
        self.assertIsInstance(u.yd, u.core.Unit)
        self.assertEqual(u.yd.scale, 0.914)
        self.assertEqual(u.yd.bases, [u.m])
        self.assertEqual(u.yd.powers, [1])

    def test_mile(self):
        self.assertIsInstance(u.mi, u.core.Unit)
        self.assertEqual(u.mi.scale, 1609.344)
        self.assertEqual(u.mi.bases, [u.m])
        self.assertEqual(u.mi.powers, [1])

    def test_second(self):
        self.assertIsInstance(u.s, u.core.Unit)
        self.assertEqual(u.s.scale, 1.0)
        self.assertEqual(u.s.bases, [u.s])
        self.assertEqual(u.s.powers, [1])

    def test_minute(self):
        self.assertIsInstance(u.min, u.core.Unit)
        self.assertEqual(u.min.scale, 60.0)
        self.assertEqual(u.min.bases, [u.s])
        self.assertEqual(u.min.powers, [1])

    def test_hour(self):
        self.assertIsInstance(u.hr, u.core.Unit)
        self.assertEqual(u.hr.scale, 3600.0)
        self.assertEqual(u.hr.bases, [u.s])
        self.assertEqual(u.hr.powers, [1])

    def test_day(self):
        self.assertIsInstance(u.dy, u.core.Unit)
        self.assertEqual(u.dy.scale, 86400.0)
        self.assertEqual(u.dy.bases, [u.s])
        self.assertEqual(u.dy.powers, [1])

    def test_julian_day_2000(self):
        self.assertIsInstance(u.jd2000, u.core.Unit)
        self.assertEqual(u.jd2000.scale, 1.0)
        self.assertEqual(u.jd2000.bases, [u.jd2000])
        self.assertEqual(u.jd2000.powers, [1])

    def test_degree(self):
        self.assertIsInstance(u.deg, u.core.Unit)
        self.assertEqual(u.deg.scale, 1.0)
        self.assertEqual(u.deg.bases, [u.deg])
        self.assertEqual(u.deg.powers, [1])

    def test_radian(self):
        self.assertIsInstance(u.rad, u.core.Unit)
        self.assertEqual(u.rad.scale, 0.01745329)
        self.assertEqual(u.rad.bases, [u.deg])
        self.assertEqual(u.rad.powers, [1])

    def test_seconds_of_arc(self):
        self.assertIsInstance(u.arcsec, u.core.Unit)
        self.assertEqual(u.arcsec.scale, 3600)
        self.assertEqual(u.arcsec.bases, [u.deg])
        self.assertEqual(u.arcsec.powers, [1])

    def test_minutes_of_arc(self):
        self.assertIsInstance(u.arcmin, u.core.Unit)
        self.assertEqual(u.arcmin.scale, 60)
        self.assertEqual(u.arcmin.bases, [u.deg])
        self.assertEqual(u.arcmin.powers, [1])

    def test_hourangle(self):
        self.assertIsInstance(u.hourangle, u.core.Unit)
        self.assertEqual(u.hourangle.scale, 1 / 15)
        self.assertEqual(u.hourangle.bases, [u.deg])
        self.assertEqual(u.hourangle.powers, [1])
