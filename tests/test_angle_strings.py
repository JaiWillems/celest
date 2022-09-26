

from celest.angle_strings import (
    _get_height_string,
    _ISO6709_representation,
    _get_latitude_string,
    _get_longitude_string,
    _get_sexagesimal_string
)
from celest.units.quantity import Quantity
from celest import units as u
from unittest import TestCase


class TestAngleStrings(TestCase):

    def test_ISO6709_representation(self):
        iso_string = _ISO6709_representation(
            Quantity(37.432, u.deg),
            Quantity(37.432, u.deg),
            Quantity(1, u.km)
        )
        self.assertEqual(iso_string, "37°25′55.20″N 37°25′55.20″E 1.00km")

    def test_get_latitude_string_in_northern_hemisphere(self):
        self.assertEqual(
            _get_latitude_string(Quantity(37.432, u.deg)),
            "37°25′55.20″N"
        )

    def test_get_latitude_string_in_southern_hemisphere(self):
        self.assertEqual(
            _get_latitude_string(Quantity(-37.432, u.deg)),
            "37°25′55.20″S"
        )

    def test_get_sexagesimal_string(self):
        self.assertEqual(_get_sexagesimal_string(37.432), "37°25′55.20″")

    def test_get_longitude_string_in_eastern_hemisphere(self):
        self.assertEqual(
            _get_longitude_string(Quantity(37.432, u.deg)),
            "37°25′55.20″E"
        )

    def test_get_longitude_string_in_western_hemisphere(self):
        self.assertEqual(
            _get_longitude_string(Quantity(-37.432, u.deg)),
            "37°25′55.20″W"
        )

    def test_get_height_string(self):
        self.assertEqual(_get_height_string(Quantity(1, u.km)), "1.00km")
