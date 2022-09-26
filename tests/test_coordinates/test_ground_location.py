

from celest.coordinates.ground_location import GroundLocation
from celest import units as u
from unittest import TestCase


class TestGroundLocation(TestCase):

    def setUp(self):
        self.latitude = -0.0944
        self.longitude = -22.6015
        self.height = 493.4238
        self.angular_unit = u.deg
        self.height_unit = u.km
        self.location = GroundLocation(
            self.latitude,
            self.longitude,
            self.height,
            self.angular_unit,
            self.height_unit
        )

    def test_initialization(self):
        self.assertIsInstance(self.location, GroundLocation)

    def test_initialization_with_invalid_latitude_raises_value_error(self):
        self.assertRaises(ValueError, GroundLocation, -91, 0, 0, u.deg, u.km)
        self.assertRaises(ValueError, GroundLocation, 91, 0, 0, u.deg, u.km)

    def test_initialization_with_invalid_longitude_raises_value_error(self):
        self.assertRaises(ValueError, GroundLocation, 0, -181, 0, u.deg, u.km)
        self.assertRaises(ValueError, GroundLocation, 0, 181, 0, u.deg, u.km)

    def test_str(self):
        self.assertEqual(
            str(self.location),
            "00°05′39.84″S 22°36′05.40″W 493.42km"
        )

    def test_repr(self):
        self.assertEqual(
            repr(self.location),
            "GroundLocation(-0.0944, -22.6015, 493.4238, Unit(\"deg\"), Unit(\"km\"))"
        )

    def test_latitude(self):
        self.assertEqual(self.location.latitude.data, self.latitude)

    def test_longitude(self):
        self.assertEqual(self.location.longitude.data, self.longitude)

    def test_height(self):
        self.assertEqual(self.location.height.data, self.height)

    def test_radius(self):
        self.assertAlmostEqual(
            self.location.radius.to(u.km),
            6871.555,
            delta=0.01
        )

    def test_itrs_x(self):
        self.assertAlmostEqual(
            self.location.itrs_x.to(u.km),
            6343.8162,
            delta=0.005
        )

    def test_itrs_y(self):
        self.assertAlmostEqual(
            self.location.itrs_y.to(u.km),
            -2640.8722,
            delta=0.005
        )

    def test_itrs_z(self):
        self.assertAlmostEqual(
            self.location.itrs_z.to(u.km),
            -11.2554,
            delta=0.005
        )
