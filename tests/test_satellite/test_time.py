

import julian as jd
import unittest
import numpy as np
from astropy import coordinates, time
from unittest import TestCase
from celest.satellite.time import Time


class TestTime(TestCase):

    def setUp(self):
        """Test fixure for test method execution."""

        self.julData = np.array([2455368.75, 2459450.85, 2456293.5416666665])
        self.astropy_time = time.Time(self.julData, format="jd")

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
        .. [2] Time Conventions. url:http://www.bom.gov.au/climate/data-services/solar/content/data-time.html.
        """

        julian = np.array([2455368.75, 2459450.85, 2456293.5416666665])
        lon = np.array([-105, -118.24, 147.46])
        tst = np.array([23.0715, 0.378059, 10.76556])

        calc_tst = Time(julian=julian).true_solar_time(longitude=lon)

        for i in range(calc_tst.size):
            with self.subTest(i=i):
                self.assertAlmostEqual(tst[i], calc_tst[i], delta=0.11)

    def test_mean_solar_time(self):
        """Test `Time.mean_solar_time`.

        Notes
        -----
        Test cases are generated from the Government of Australia.[1]_

        References
        ----------
        .. [1] Time Conventions. url: http://www.bom.gov.au/climate/data-services/solar/content/data-time.html.
        """

        julData = np.array([2456293.5416666665])
        lon = np.array([147.46])
        tst = np.array([10.83056])

        calc_tst = Time(julian=julData).mean_solar_time(longitude=lon)

        for i in range(calc_tst.size):
            with self.subTest(i=i):
                self.assertAlmostEqual(tst[i], calc_tst[i], delta=0.001)

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

        julData = np.array([2455368.75, 2459450.85, 2456293.5416666665])
        lon = np.array([-105, -118.24, 147.46])
        hour_angle = np.array([166.0734, 185.671, 341.53]) / 15 % 24

        calc_hour_angle = Time(julian=julData).true_hour_angle(longitude=lon)

        for i in range(calc_hour_angle.size):
            with self.subTest(i=i):
                self.assertAlmostEqual(hour_angle[i], calc_hour_angle[i], delta=1.51)

    def test_mean_hour_angle(self):
        """Test `Time.true_hour_angle`.

        Notes
        -----
        Test cases are self generated using the definition of the mean hour
        angle.
        """

        julData = np.array([2456293.5416666665])
        lon = np.array([147.46])
        hour_angle = np.array([342.4584]) / 15 % 24

        calc_hour_angle = Time(julian=julData).mean_hour_angle(longitude=lon)

        for i in range(calc_hour_angle.size):
            with self.subTest(i=i):
                self.assertAlmostEqual(hour_angle[i], calc_hour_angle[i], delta=0.01)

    def test_ut1(self):
        """Test `Time.UT1`.

        Notes
        -----
        Test cases are generated using the `Astropy` Python package.
        """

        dt = self.astropy_time.get_delta_ut1_utc().value / 3600

        ut1 = self.astropy_time.to_value("decimalyear") % 1
        ut1 = (ut1 * 365 * 24) % 24 + dt

        calc_ut1 = Time(julian=self.julData).ut1()

        for i in range(calc_ut1.size):
            with self.subTest(i=i):
                self.assertAlmostEqual(ut1[i], calc_ut1[i], delta=0.001)

    def test_julian(self):
        """Test `Time.julian`."""

        calc_julian = Time(julian=self.julData).julian()

        for i in range(calc_julian.size):
            with self.subTest(i=i):
                self.assertEqual(self.julData[i], calc_julian[i])

    def test_datetime(self):
        """Test `Time.datetime`."""

        calc_datetime = Time(julian=self.julData).datetime()

        for i in range(calc_datetime.size):
            with self.subTest(i=i):
                dt = jd.from_jd(self.julData[i])

                self.assertEqual(calc_datetime[i].year, dt.year)
                self.assertEqual(calc_datetime[i].month, dt.month)
                self.assertEqual(calc_datetime[i].day, dt.day)
                self.assertEqual(calc_datetime[i].second, dt.second)
                self.assertAlmostEqual(calc_datetime[i].microsecond, dt.microsecond, delta=1)

    def test_gmst(self):
        """Test `Time.gmst`.

        Notes
        -----
        Test cases are generated using the `Astropy` Python package.
        """

        gmst = self.astropy_time.sidereal_time("mean", "greenwich")
        gmst = coordinates.Angle(gmst).hour

        calc_gmst = Time(julian=self.julData).gmst()

        for i in range(calc_gmst.size):
            with self.subTest(i=i):
                self.assertAlmostEqual(gmst[i], calc_gmst[i], delta=0.0001)

    def test_lmst(self):
        """Test `Time.lmst`.

        Notes
        -----
        Test cases are generated using the `Astropy` Python package.
        """

        lmst = self.astropy_time.sidereal_time("mean", longitude="150")
        lmst = coordinates.Angle(lmst).hour

        calc_lmst = Time(julian=self.julData).lmst(longitude=150)

        for i in range(calc_lmst.size):
            with self.subTest(i=i):
                self.assertAlmostEqual(lmst[i], calc_lmst[i], delta=0.1)

    def test_gast(self):
        """Test `Time.gast`.

        Notes
        -----
        Test cases are generated using the `Astropy` Python package.
        """

        gast = self.astropy_time.sidereal_time("apparent", "greenwich")
        gast = coordinates.Angle(gast).hour

        calc_gast = Time(julian=self.julData).gast()

        for i in range(calc_gast.size):
            with self.subTest(i=i):
                self.assertAlmostEqual(gast[i], calc_gast[i], delta=0.0001)

    def test_last(self):
        """Test `Time.last`.

        Notes
        -----
        Test cases are generated using the `Astropy` Python package.
        """

        last = self.astropy_time.sidereal_time("apparent", longitude="150")
        last = coordinates.Angle(last).hour

        calc_last = Time(julian=self.julData).last(longitude=150)

        for i in range(calc_last.size):
            with self.subTest(i=i):
                self.assertAlmostEqual(last[i], calc_last[i], delta=0.1)


if __name__ == "__main__":
    unittest.main()
