"""Testing module for the window utility functions."""


import numpy as np
import unittest
from unittest import TestCase
from celest.encounter.groundposition import GroundPosition
from celest.encounter._window_utils import _sun_itrs, _get_ang, _window_encounter_ind
from celest.satellite.coordinate import Coordinate
from celest.satellite.satellite import Satellite
from celest.satellite.time import Time


class TestWindowUtils(TestCase):

    def setUp(self):
        """Test fixure for test method execution."""

        fname = "tests/test_data/coordinate_validation_set.txt"
        data = np.loadtxt(fname=fname, delimiter="\t", skiprows=1)
        times, itrs = data[:, 0], data[:, 10:]

        timeData = Time(times, 2430000)
        coor = Coordinate(itrs, "itrs", timeData)

        self.finch = Satellite(coor)

    def test_sun_itrs(self):
        """Test `Encounter._sun_itrs`.

        Notes
        -----
        Test cases are generated using the `Astropy` Python package.
        """
        
        from astropy.coordinates import get_body, ITRS
        from astropy import time, units

        julData = np.array([2455368.75, 2459450.85, 2456293.5416666665])
        astropy_time = time.Time(julData, format="jd")
        
        sun_pos = get_body("sun", astropy_time)
        sun_pos.representation_type="cartesian"

        itrsCoor = sun_pos.transform_to(ITRS(obstime=astropy_time))

        x = itrsCoor.x.to(units.km).value
        y = itrsCoor.y.to(units.km).value
        z = itrsCoor.z.to(units.km).value

        calc_sun_pos = _sun_itrs(julData)
        self.assertTrue(np.allclose(x, calc_sun_pos[:, 0], rtol=0.05))
        self.assertTrue(np.allclose(y, calc_sun_pos[:, 1], rtol=0.05))
        self.assertTrue(np.allclose(z, calc_sun_pos[:, 2], rtol=0.05))

    def test_get_ang(self):
        """Test `Encounter._get_ang`."""

        vec_one = np.array([[56, 92, 76], [9238, 8479, 9387], [2, 98, 23]])
        vec_two = np.array([[36, 29, 38], [2703, 947, 8739], [9827, 921, 1]])
        ang = np.array([16.28, 37, 83.65])

        calc_ang = _get_ang(vec_one, vec_two)

        for i in range(calc_ang.size):
            with self.subTest(i=i):
                self.assertAlmostEqual(ang[i], calc_ang[i], delta=0.01)

    def test_window_encounter_ind(self):
        """Test `Encounter._window_encounter_ind`."""

        julData = self.finch.time.julian()

        # Check data link index generation.
        location = GroundPosition(43.6532, -79.3832)

        calc_enc_ind = _window_encounter_ind(self.finch, location, 30, 1, 0, 1)

        altitude, _ = self.finch.position.horizontal(location=location)
        off_nadir = self.finch.position.off_nadir(location=location)

        sun_itrs = _sun_itrs(julData)
        sun_pos = Coordinate(sun_itrs, "itrs", self.finch.time)
        sun_alt, _ = sun_pos.horizontal(location=location)

        enc_ind = np.arange(0, julData.shape[0], 1)
        
        alt_ind = np.where(altitude >= 0)[0]
        off_nadir_ind = np.where(off_nadir <= 30)
        sun_ind = np.where(sun_alt >= 0)
        
        enc_ind = np.intersect1d(enc_ind, alt_ind)
        enc_ind = np.intersect1d(enc_ind, off_nadir_ind)
        enc_ind = np.intersect1d(enc_ind, sun_ind)

        self.assertTrue(np.array_equal(calc_enc_ind, enc_ind))

        # Check data link index generation.
        location = GroundPosition(43.6532, -79.3832)

        calc_enc_ind = _window_encounter_ind(self.finch, location, 30, 0, 30, 0)

        altitude, _ = self.finch.position.horizontal(location=location)

        enc_ind = np.arange(0, julData.shape[0], 1)
        
        alt_ind = np.where(30 <= altitude)[0]

        sat_itrs = self.finch.position.itrs()

        ground_GEO = [location.lat, location.lon]
        ground_GEO = np.repeat(np.array([ground_GEO]), julData.size, axis=0)
        gnd_itrs = Coordinate(ground_GEO, "geo", self.finch.time).itrs()

        sca_angs = _get_ang(sat_itrs - gnd_itrs, sun_itrs - gnd_itrs)
        sca_ind = np.where(sca_angs > 30)[0]
        
        enc_ind = np.intersect1d(enc_ind, alt_ind)
        enc_ind = np.intersect1d(enc_ind, sca_ind)

        self.assertTrue(np.array_equal(calc_enc_ind, enc_ind))
