

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
        gcrs_position = data[:, 1:4]
        gcrs_velocity = data[:, 4:7]

        offset = 2430000
        self.satellite = Satellite(gcrs_position, gcrs_velocity, "gcrs", times,
                                   offset)

    def test_sun_coor(self):

        from astropy.coordinates import get_body, ITRS
        from astropy import time, units

        julian = np.array([2455368.75000, 2456293.54167, 2459450.85000])
        astropy_time = time.Time(julian, format="jd")

        sun_position = get_body("sun", astropy_time)
        sun_position.representation_type = "cartesian"

        itrsCoor = sun_position.transform_to(ITRS(obstime=astropy_time))

        true_x = itrsCoor.x.to(units.km).value
        true_y = itrsCoor.y.to(units.km).value
        true_z = itrsCoor.z.to(units.km).value

        test_x, test_y, test_z, _, _, _ = _sun_coor(julian).itrs(False)
        self.assertTrue(np.allclose(true_x, test_x, rtol=0.05))
        self.assertTrue(np.allclose(true_y, test_y, rtol=0.05))
        self.assertTrue(np.allclose(true_z, test_z, rtol=0.05))

    def test_windows(self):

        location = GroundPosition(43.6532, -79.3832)
        elevation, _ = self.satellite.horizontal(location, True)

        sun_coordinates = _sun_coor(self.satellite._julian)
        sun_elevation, _ = sun_coordinates.horizontal(location, True)

        vis_threshold, lighting, tol = 10, 0, 1e-5
        windows = generate_vtw(self.satellite, location, vis_threshold, lighting, tol)

        for window in windows:
            rise_time, set_time = window.rise_time, window.set_time
            window_times = np.linspace(rise_time, set_time, 10)

            self.assertTrue(np.all(elevation(window_times) > vis_threshold))

        vis_threshold, lighting, tol = 10, 1, 1e-5
        windows = generate_vtw(self.satellite, location, vis_threshold, lighting, tol)

        for window in windows:
            rise_time, set_time = window.rise_time, window.set_time
            window_times = np.linspace(rise_time, set_time, 10)
            sun_elev = sun_elevation(window_times)

            self.assertTrue(np.all(elevation(window_times) > vis_threshold))
            self.assertTrue(np.all((0 < sun_elev)))

        vis_threshold, lighting, tol = 10, -1, 1e-5
        windows = generate_vtw(self.satellite, location, vis_threshold, lighting, tol)

        for window in windows:
            rise_time, set_time = window.rise_time, window.set_time
            window_times = np.linspace(rise_time, set_time, 10)
            sun_elev = sun_elevation(window_times)

            self.assertTrue(np.all(elevation(window_times) > vis_threshold))
            self.assertTrue(np.all((sun_elev <= 0)))


if __name__ == "__main__":
    unittest.main()
