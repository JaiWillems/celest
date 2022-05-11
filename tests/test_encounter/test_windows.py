

from celest.encounter.groundposition import GroundPosition
from celest.encounter.windows import _sun_coor, generate_vtw
from celest.satellite.satellite import Satellite
from unittest import TestCase
import numpy as np
import unittest


class TestWindows(TestCase):

    def setUp(self):

        fname = "tests/test_data/coordinate_validation_long.txt"
        cols = (0, 11, 12, 13, 14, 15, 16)
        skiprows = 1
        max_rows = 5000
        data = np.loadtxt(fname=fname, usecols=cols, skiprows=skiprows,
                          max_rows=max_rows)

        times = data[:, 0]
        GCRS = data[:, 1:4]
        GCRS_vel = data[:, 4:7]

        offset = 2430000
        self.satellite = Satellite(GCRS, GCRS_vel, "gcrs", times, offset)

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

        calc_x, calc_y, calc_z, _, _, _ = _sun_coor(julData).itrs(stroke=False)
        self.assertTrue(np.allclose(x, calc_x, rtol=0.05))
        self.assertTrue(np.allclose(y, calc_y, rtol=0.05))
        self.assertTrue(np.allclose(z, calc_z, rtol=0.05))

    def test_windows(self):

        location = GroundPosition(43.6532, -79.3832)
        elevation, _ = self.satellite.horizontal(location, stroke=True)

        sun_coor = _sun_coor(self.satellite._julian)
        sun_alt, _ = sun_coor.horizontal(location, stroke=True)

        # Test the window generator for all day case.

        vis_threshold, lighting, tol = 10, 0, 1e-5
        windows = generate_vtw(self.satellite, location, vis_threshold, lighting, tol)

        for window in windows:
            rise_time, set_time = window.rise_time, window.set_time
            window_times = np.linspace(rise_time, set_time, 10)

            self.assertTrue(np.all(elevation(window_times) > vis_threshold))

        # Test the window generator for day only case.

        vis_threshold, lighting, tol = 10, 1, 1e-5
        windows = generate_vtw(self.satellite, location, vis_threshold, lighting, tol)

        for window in windows:
            rise_time, set_time = window.rise_time, window.set_time
            window_times = np.linspace(rise_time, set_time, 10)
            sa = sun_alt(window_times)

            self.assertTrue(np.all(elevation(window_times) > vis_threshold))
            self.assertTrue(np.all((0 < sa)))

        # Test the window generator for night only case.

        vis_threshold, lighting, tol = 10, -1, 1e-5
        windows = generate_vtw(self.satellite, location, vis_threshold, lighting, tol)

        for window in windows:
            rise_time, set_time = window.rise_time, window.set_time
            window_times = np.linspace(rise_time, set_time, 10)
            sa = sun_alt(window_times)

            self.assertTrue(np.all(elevation(window_times) > vis_threshold))
            self.assertTrue(np.all((sa <= 0)))


if __name__ == "__main__":
    unittest.main()
