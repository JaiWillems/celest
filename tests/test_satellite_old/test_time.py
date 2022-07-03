

from astropy import coordinates, time
from celest.satellite.time import Time
from unittest import TestCase
import julian as jd
import unittest
import numpy as np


class TestTime(TestCase):

    def setUp(self):

        self.julian = [2455368.75, 2459450.85, 2456293.5416666665]
        self.astropy_time = time.Time(self.julian, format="jd")

    def test_true_solar_time(self):
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

        true_true_solar_time = np.array([23.0715, 0.378059, 10.76556])
        test_true_solar_time = Time(julian).true_solar_time(longitude)

        self.assertTrue(np.allclose(true_true_solar_time, test_true_solar_time,
                        atol=0.11))

    def test_mean_solar_time(self):
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

        true_mean_solar_time = np.array([10.83056])
        test_mean_solar_time = Time(julian).mean_solar_time(longitude)

        self.assertTrue(np.allclose(true_mean_solar_time, test_mean_solar_time,
                        atol=0.001))

    def test_true_hour_angle(self):
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

        true_true_hour_angle = np.array([166.0734, 185.671, 341.53]) / 15 % 24
        test_true_hour_angle = Time(julian).true_hour_angle(longitude)

        self.assertTrue(np.allclose(true_true_hour_angle, test_true_hour_angle,
                        atol=1.51))

    def test_mean_hour_angle(self):

        julian = 2456293.5416666665
        longitude = 147.46

        true_mean_hour_angle = np.array([342.4584]) / 15 % 24
        test_mean_hour_angle = Time(julian).mean_hour_angle(longitude)

        self.assertTrue(np.allclose(true_mean_hour_angle, test_mean_hour_angle,
                        atol=0.01))

    def test_ut1(self):

        ut1_utc_diff = self.astropy_time.get_delta_ut1_utc().value / 3600

        true_ut1 = self.astropy_time.to_value("decimalyear") % 1
        true_ut1 = (true_ut1 * 365 * 24) % 24 + ut1_utc_diff
        test_ut1 = Time(self.julian).ut1()

        self.assertTrue(np.allclose(true_ut1, test_ut1, atol=0.001))

    def test_julian(self):

        test_julian = Time(julian=self.julian).julian()
        self.assertTrue(np.array_equal(self.julian, test_julian))

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

        true_gmst = self.astropy_time.sidereal_time("mean", "greenwich")
        true_gmst = coordinates.Angle(true_gmst).hour
        test_gmst = Time(self.julian).gmst()

        self.assertTrue(np.allclose(true_gmst, test_gmst, atol=0.0001))

    def test_lmst(self):

        true_lmst = self.astropy_time.sidereal_time("mean", longitude="150")
        true_lmst = coordinates.Angle(true_lmst).hour
        test_lmst = Time(self.julian).lmst(150)

        self.assertTrue(np.allclose(true_lmst, test_lmst, atol=0.1))

    def test_gast(self):

        true_gast = self.astropy_time.sidereal_time("apparent", "greenwich")
        true_gast = coordinates.Angle(true_gast).hour
        test_gast = Time(self.julian).gast()

        self.assertTrue(np.allclose(true_gast, test_gast, atol=0.0001))

    def test_last(self):

        true_last = self.astropy_time.sidereal_time("apparent", longitude="150")
        true_last = coordinates.Angle(true_last).hour
        test_last = Time(self.julian).last(150)
        
        self.assertTrue(np.allclose(true_last, test_last, atol=0.1))


if __name__ == "__main__":
    unittest.main()
