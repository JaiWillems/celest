

from astropy import coordinates, time
from celest.time import Time
from celest import units as u
from unittest import TestCase
import julian as jd
import unittest
import numpy as np


class TestTime(TestCase):

    def setUp(self):
        self.julian = [2455368.75, 2459450.85, 2456293.5416666665]
        self.astropy_time = time.Time(self.julian, format="jd")

    def test_julian(self):
        actual_julian = Time(julian=self.julian).julian.data
        self.assertTrue(np.array_equal(self.julian, actual_julian))

    def test_true_solar_time_for_validation(self):
        """Test `Time.true_solar_time`.

        Notes
        -----
        Test cases are generated from the Global Monitoring Labratory and the
        Government of Australia.[1]_[2]_

        References
        ----------
        .. [1] NOAA US Department of Commerce. ESRL Global Monitoring
           Laboratory -Global Radiation and Aerosols. url:
           https://gml.noaa.gov/grad/solcalc/.
        .. [2] Time Conventions. url:
           http://www.bom.gov.au/climate/data-services/solar/content/data-time.html.
        """

        julian = [2455368.75, 2459450.85, 2456293.5416666665]
        longitude = [-105, -118.24, 147.46]

        expected_true_solar_time = np.array([23.0715, 0.378059, 10.76556])
        actual_true_solar_time = Time(julian).true_solar_time(longitude).to(u.hourangle).data

        self.assertTrue(np.allclose(expected_true_solar_time,
                                    actual_true_solar_time, atol=0.11))

    def test_mean_solar_time_for_validation(self):
        """Test `Time.mean_solar_time`.

        Notes
        -----
        Test cases are generated from the Government of Australia.[1]_

        References
        ----------
        .. [1] Time Conventions. url:
           http://www.bom.gov.au/climate/data-services/solar/content/data-time.html.
        """

        julian = 2456293.5416666665
        longitude = 147.46

        expected_mean_solar_time = np.array([10.83056])
        actual_mean_solar_time = Time(julian).mean_solar_time(longitude).to(u.hourangle).data

        self.assertTrue(np.allclose(expected_mean_solar_time,
                                    actual_mean_solar_time, atol=0.001))

    def test_true_hour_angle_for_validation(self):
        """Test `Time.true_hour_angle`.

        Notes
        -----
        Test cases are taken from Global Monitoring Labratory.[1]_

        References
        ----------
        .. [1] NOAA US Department of Commerce. ESRL Global Monitoring
           Laboratory -Global Radiation and Aerosols. url:
           https://gml.noaa.gov/grad/solcalc/.
        """

        julian = [2455368.75, 2459450.85, 2456293.5416666665]
        longitude = [-105, -118.24, 147.46]

        expected_true_hour_angle = np.array([166.0734, 185.671, 341.53]) / 15 % 24
        actual_true_hour_angle = Time(julian).true_hour_angle(longitude).to(u.hourangle).data

        self.assertTrue(np.allclose(expected_true_hour_angle,
                                    actual_true_hour_angle, atol=1.51))

    def test_mean_hour_angle_for_validation(self):
        julian = 2456293.5416666665
        longitude = 147.46

        true_mean_hour_angle = np.array([342.4584]) / 15 % 24
        test_mean_hour_angle = Time(julian).mean_hour_angle(longitude).to(u.hourangle).data

        self.assertTrue(np.allclose(true_mean_hour_angle, test_mean_hour_angle,
                        atol=0.01))

    def test_ut1_for_validation(self):
        ut1_utc_diff = self.astropy_time.get_delta_ut1_utc().value / 3600

        expected_ut1 = self.astropy_time.to_value("decimalyear") % 1
        expected_ut1 = (expected_ut1 * 365 * 24) % 24 + ut1_utc_diff
        actual_ut1 = Time(self.julian).ut1().to(u.hourangle).data

        self.assertTrue(np.allclose(expected_ut1, actual_ut1, atol=0.001))

    def test_datetime(self):
        test_datetime = Time(julian=self.julian).datetime()

        for julian, test_datetime_instace in zip(self.julian, test_datetime):
            true_datetime_instance = jd.from_jd(julian)

            self.assertEqual(test_datetime_instace.year,
                             true_datetime_instance.year)
            self.assertEqual(test_datetime_instace.month,
                             true_datetime_instance.month)
            self.assertEqual(test_datetime_instace.day,
                             true_datetime_instance.day)
            self.assertEqual(test_datetime_instace.second,
                             true_datetime_instance.second)
            self.assertAlmostEqual(test_datetime_instace.microsecond,
                                   true_datetime_instance.microsecond, delta=1)

    def test_gmst(self):
        expected_gmst = self.astropy_time.sidereal_time("mean", "greenwich")
        expected_gmst = coordinates.Angle(expected_gmst).hour
        actual_gmst = Time(self.julian).gmst().to(u.hourangle).data

        self.assertTrue(np.allclose(expected_gmst, actual_gmst, atol=0.0001))

    def test_lmst(self):
        expected_lmst = self.astropy_time.sidereal_time("mean", longitude="150")
        expected_lmst = coordinates.Angle(expected_lmst).hour
        actual_lmst = Time(self.julian).lmst(150).to(u.hourangle).data

        self.assertTrue(np.allclose(expected_lmst, actual_lmst, atol=0.1))

    def test_gast(self):
        expected_gast = self.astropy_time.sidereal_time("apparent", "greenwich")
        expected_gast = coordinates.Angle(expected_gast).hour
        actual_gast = Time(self.julian).gast().to(u.hourangle).data

        self.assertTrue(np.allclose(expected_gast, actual_gast, atol=0.0001))

    def test_last(self):
        expected_last = self.astropy_time.sidereal_time("apparent", longitude="150")
        expected_last = coordinates.Angle(expected_last).hour
        actual_last = Time(self.julian).last(150).to(u.hourangle).data
        
        self.assertTrue(np.allclose(expected_last, actual_last, atol=0.1))


if __name__ == "__main__":
    unittest.main()
