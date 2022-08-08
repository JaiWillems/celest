

from celest.coordinates.coordinate import Coordinate
from celest.coordinates.frames.azel import AzEl
from celest.coordinates.frames.gcrs import GCRS
from celest.coordinates.ground_location import GroundLocation
from celest.encounter.window_generator import (
    generate_vtws,
    _get_day_constraint_indices,
    _get_night_constraint_indices,
    Lighting,
    _get_sun_gcrs
)
from celest.satellite import Satellite
from celest import units as u
from unittest import TestCase
import numpy as np


class TestGenerateVTW(TestCase):

    def setUp(self):
        data = np.loadtxt(
            fname="tests/test_data/coordinate_validation_long.txt",
            usecols=(0, 11, 12, 13, 14, 15, 16),
            skiprows=1,
            max_rows=5000
        )

        self.julian = data[:, 0] + 2430000
        gcrs_position = data[:, 1:4]
        gcrs_velocity = data[:, 4:]

        self.gcrs_position = GCRS(
            self.julian,
            gcrs_position[:, 0],
            gcrs_position[:, 1],
            gcrs_position[:, 2],
            u.km
        )
        self.gcrs_velocity = GCRS(
            self.julian,
            gcrs_velocity[:, 0],
            gcrs_velocity[:, 1],
            gcrs_velocity[:, 2],
            u.m / u.s
        )
        self.satellite = Satellite(self.gcrs_position, self.gcrs_velocity)
        self.location = GroundLocation(52.1579, -106.6702, 0.482, u.deg, u.km)

    def test_generate_vtw_raises_warning_when_satellite_has_no_velocity(self):
        self.assertRaises(Warning, generate_vtws, Satellite(self.gcrs_position),
                          self.location, 10)

    def test_generate_vtw_raises_value_error_for_negative_threshold(self):
        self.assertRaises(ValueError, generate_vtws, self.satellite,
                          self.location, -10)

    def test_generate_vtw_raises_value_error_for_large_threshold(self):
        self.assertRaises(ValueError, generate_vtws, self.satellite,
                          self.location, 100)

    def test_window_constraints_valid_at_day(self):
        vis_threshold = 10
        windows = generate_vtws(
            self.satellite,
            self.location,
            vis_threshold,
            Lighting.DAYTIME
        )
        satellite_elevation = self._get_satellite_elevation().to(u.deg).data
        sun_elevation = self._get_sun_elevation().to(u.deg).data

        for window in windows:
            rise_time, set_time = window.rise_time.data, window.set_time.data
            window_indices = np.where((rise_time <= self.julian) &
                                      (self.julian <= set_time))[0]

            # The first and last indices are always invalid.
            self.assertTrue(np.alltrue(
                satellite_elevation[window_indices][1:-1] > vis_threshold))
            self.assertTrue(np.alltrue(sun_elevation[window_indices][1:-1] > 0))

    def _get_satellite_elevation(self):
        satellite_coordinates = Coordinate(self.satellite.position)
        return satellite_coordinates.convert_to(AzEl, self.location).elevation

    def _get_sun_elevation(self):
        sun_gcrs = _get_sun_gcrs(self.julian)
        sun_azel = Coordinate(sun_gcrs).convert_to(AzEl, self.location)
        return sun_azel.elevation

    def test_window_constraints_valid_at_night(self):
        vis_threshold = 10
        windows = generate_vtws(
            self.satellite,
            self.location,
            vis_threshold,
            Lighting.NIGHTTIME
        )
        satellite_elevation = self._get_satellite_elevation().to(u.deg).data
        sun_elevation = self._get_sun_elevation().to(u.deg).data

        for window in windows:
            rise_time, set_time = window.rise_time.data, window.set_time.data
            window_indices = np.where((rise_time <= self.julian) &
                                      (self.julian <= set_time))[0]

            # The first and last indices are always invalid.
            self.assertTrue(np.alltrue(
                satellite_elevation[window_indices][1:-1] > vis_threshold))
            self.assertTrue(np.alltrue(sun_elevation[window_indices][1:-1] < 0))

    def test_get_night_constraint_indices(self):
        indices = _get_night_constraint_indices(self.julian, self.location)
        sun_elevation = self._get_sun_elevation()
        self.assertTrue(np.alltrue(sun_elevation.to(u.deg).data[indices] < 0))

    def test_get_sun_gcrs(self):
        from astropy.coordinates import get_body
        from astropy import time, units

        julian = np.array([2455368.75000, 2456293.54167, 2459450.85000])
        astropy_time = time.Time(julian, format="jd")

        sun_position = get_body("sun", astropy_time)
        sun_position.representation_type = "cartesian"

        expected_gcrs = sun_position.gcrs
        actual_gcrs = _get_sun_gcrs(julian)

        self.assertTrue(np.allclose(expected_gcrs.x.to(units.km).value,
                                    actual_gcrs.x.to(u.km).data, rtol=0.05))
        self.assertTrue(np.allclose(expected_gcrs.y.to(units.km).value,
                                    actual_gcrs.y.to(u.km).data, rtol=0.05))
        self.assertTrue(np.allclose(expected_gcrs.z.to(units.km).value,
                                    actual_gcrs.z.to(u.km).data, rtol=0.05))

    def test_get_day_constraint_indices(self):
        sun_elevation = self._get_sun_elevation()
        indices = _get_day_constraint_indices(self.julian, self.location)
        self.assertTrue(np.alltrue(sun_elevation.to(u.deg).data[indices] > 0))
