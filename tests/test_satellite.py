

from celest.coordinates.frames.attitude import Attitude
from celest.coordinates.frames.gcrs import GCRS
from celest.coordinates.frames.itrs import ITRS
from celest.coordinates.frames.wgs84 import WGS84
from celest.coordinates.ground_location import GroundLocation
from celest.satellite import Satellite
from celest import units as u
from unittest import TestCase
import numpy as np
import os


class TestSatellite(TestCase):

    def setUp(self):
        data = np.loadtxt(
            fname="tests/test_data/coordinate_validation_long.txt",
            usecols=(0, 5, 6, 7, 11, 12, 13, 14, 15, 16, 17, 18, 19),
            skiprows=1,
            max_rows=5000
        )

        self.julian = data[:, 0] + 2430000
        self.wgs84 = data[:, 1:3]
        self.altitude = data[:, 3]
        self.gcrs = data[:, 4:7]
        self.gcrs_velocity = data[:, 7:10]
        self.itrs = data[:, 10:]

        self.gcrs_position = GCRS(
            self.julian,
            self.gcrs[:, 0],
            self.gcrs[:, 1],
            self.gcrs[:, 2],
            u.km
        )
        self.gcrs_velocity = GCRS(
            self.julian,
            self.gcrs_velocity[:, 0],
            self.gcrs_velocity[:, 1],
            self.gcrs_velocity[:, 2],
            u.m / u.s
        )

        self.satellite = Satellite(self.gcrs_position, self.gcrs_velocity)

        self.location = GroundLocation(52.1579, -106.6702, 0.482, u.deg, u.km)

    def test_initialization(self):
        self.assertIsInstance(self.satellite, Satellite)

    def test_value_error_raised_from_improper_data_frames(self):
        wgs84_position = WGS84(
            self.julian,
            self.wgs84[:, 0],
            self.wgs84[:, 1],
            self.altitude,
            u.deg,
            u.km
        )
        self.assertRaises(ValueError, Satellite, wgs84_position)

    def test_attitude(self):
        self.satellite.attitude(self.location)
        self.assertIsInstance(self.satellite.attitude(self.location), Attitude)

    def test_look_angle(self):
        """
        Notes
        -----
        Test cases generated using a geometric model designed by Mingde Yin.

        The tolerance required to pass all test cases is artificially inflated
        due to the low precision in the model's output and parameters used to
        generate the test cases.
        """

        gcrs = GCRS(
            self.julian[210:220],
            self.gcrs[210:220, 0],
            self.gcrs[210:220, 1],
            self.gcrs[210:220, 2],
            u.km
        )
        satellite = Satellite(gcrs)

        expected_look_angle = np.array([66.88, 65.09, 63.90, 63.22, 62.46,
                                        61.67, 58.42, 38.27, 23.73, 56.29])
        actual_look_angle = satellite.look_angle(self.location).to(u.deg).data

        self.assertTrue(np.allclose(actual_look_angle, expected_look_angle,
                                    atol=0.25))

    def test_altitude(self):
        self.assertTrue(np.allclose(self.satellite.altitude().data,
                                    self.altitude, atol=0.1))

    def test_distance(self):
        julian = np.array(
            [30462.5, 30462.50069, 30462.50172, 30462.50279, 30462.50386])
        itrs_x = np.array(
            [6343.81620, 6295.64584, 6173.04658, 5980.91111, 5724.30203])
        itrs_y = np.array(
            [-2640.87223, -2718.09271, -2808.91831, -2872.96258, -2904.09810])
        itrs_z = np.array(
            [-11.25542802, 443.08233, 1108.58544, 1792.17964, 2460.25799])

        itrs = ITRS(julian, itrs_x, itrs_y, itrs_z, u.km)
        satellite = Satellite(itrs)

        expected_distance = np.array(
            [9070.49269, 8776.71795, 8330.99543, 7851.70083, 7359.09190])
        actual_distance = satellite.distance(self.location).to(u.km).data

        self.assertTrue(np.allclose(expected_distance, actual_distance,
                                    atol=0.5))

    def test_save_text_file(self):
        file_name = "satellite_test_file"
        self.satellite.save_text_file(file_name)
        self.assertTrue(os.path.exists(file_name + ".txt"))
        if os.path.exists(file_name + ".txt"):
            os.remove(file_name + ".txt")
