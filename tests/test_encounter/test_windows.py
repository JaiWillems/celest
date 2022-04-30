

from celest.encounter.groundposition import GroundPosition
from celest.encounter.windows import _get_ang, _sun_coor, generate
from celest.satellite.satellite import Satellite
from celest.satellite.time import Time
from unittest import TestCase
import numpy as np
import unittest


class TestWindows(TestCase):

    def setUp(self):

        fname = "tests/test_data/coordinate_validation_set.txt"
        data = np.loadtxt(fname=fname, delimiter="\t", skiprows=1)
        times, itrs = data[:, 0], data[:, 10:]

        self.sat = Satellite(itrs, "itrs", times, 2430000)

    def test_sun_coor(self):

        from astropy.coordinates import get_body, ITRS
        from astropy import time, units

        julData = np.array([2455368.75, 2456293.5416666665, 2459450.85])
        astropy_time = time.Time(julData, format="jd")

        sun_pos = get_body("sun", astropy_time)
        sun_pos.representation_type = "cartesian"

        itrsCoor = sun_pos.transform_to(ITRS(obstime=astropy_time))

        x = itrsCoor.x.to(units.km).value
        y = itrsCoor.y.to(units.km).value
        z = itrsCoor.z.to(units.km).value

        calc_x, calc_y, calc_z = _sun_coor(julData).itrs(stroke=False)
        self.assertTrue(np.allclose(x, calc_x, rtol=0.05))
        self.assertTrue(np.allclose(y, calc_y, rtol=0.05))
        self.assertTrue(np.allclose(z, calc_z, rtol=0.05))

    def test_get_ang(self):
        """Test `Encounter._get_ang`."""

        vec_one = np.array([[56, 92, 76], [9238, 8479, 9387], [2, 98, 23]])
        vec_two = np.array([[36, 29, 38], [2703, 947, 8739], [9827, 921, 1]])
        ang = [16.28, 37, 83.65]

        calc_ang = _get_ang(vec_one, vec_two)

        for i in range(calc_ang.size):
            with self.subTest(i=i):
                self.assertAlmostEqual(ang[i], calc_ang[i], delta=0.01)

    def test_windows(self):

        location = GroundPosition(43.6532, -79.3832)

        elevation = self.sat._elevation(location)
        off_nadir = self.sat.off_nadir(location, stroke=True)

        # Test the window generator for day imaging case.

        enc, ang, lighting, tol = "image", 30, 1, 1e-5
        windows = generate(self.sat, location, enc, ang, lighting, tol)

        for window in windows:
            start, end = window.start, window.end
            window_times = np.linspace(start, end, 10)

            tha = Time(window_times).true_hour_angle(location.lon)

            self.assertTrue(np.all(elevation(window_times) > 0))
            self.assertTrue(np.all(off_nadir(window_times) < ang))
            self.assertTrue(np.all((-90 < tha) & (tha < 90)))

        # Test the window generator for all day data link case.

        enc, ang, lighting, tol = "data_link", 10, 0, 1e-5
        windows = generate(self.sat, location, enc, ang, lighting, tol)

        for window in windows:
            start, end = window.start, window.end
            window_times = np.linspace(start, end, 10)

            tha = Time(window_times).true_hour_angle(location.lon)

            self.assertTrue(np.all(elevation(window_times) > 10))


if __name__ == "__main__":
    unittest.main()
