

from celest.satellite._astronomical_quantities import *
from unittest import TestCase
import julian as jd
import numpy as np
import unittest


class TestAstronomicalQuantities(TestCase):

    def test_nutation_angles(self):
        """Test `AstronomicalQuantities.nutation_angles`.

        Notes
        -----
        Test case is taken from Astronomical Algorithms.[1]_

        References
        ----------
        .. [1] Jean Meeus. Astronomical algorithms. 2nd ed. Willmann-Bell,
           1998, pp. 148. isbn: 9780943396613.
        """

        julian = np.array([2446895.5])
        D, M, N, F, O = nutation_angles(julian)

        self.assertAlmostEqual(D[0], 136.9623, places=4)
        self.assertAlmostEqual(M[0], 94.9792, places=4)
        self.assertAlmostEqual(N[0], 229.2784, places=4)
        self.assertAlmostEqual(F[0], 143.4079, places=4)
        self.assertAlmostEqual(O[0], 11.2531, places=4)

    def test_nutation_components(self):
        """Test `AstronomicalQuantities.nutation_components`.

        Notes
        -----
        Test cases are taken from the PHP Science Labs.[1]_

        References
        ----------
        .. [1] Jay Tanner. NeoProgrammics - Science Computations. 2021.
           url: http://www.neoprogrammics.com/nutations/.
        """

        julian = np.array([2449634.50000, 2453420.56250, 2477418.21181])

        true_lon_nutation = np.array([11.694, -5.993, 2.937])
        true_obliquity_nutation = np.array([-5.946, 8.431, -8.871])
        test_lon_nutation, test_obliquity_nutation = nutation_components(julian)

        self.assertTrue(np.allclose(test_lon_nutation, true_lon_nutation,
                        atol=0.5))
        self.assertTrue(np.allclose(test_obliquity_nutation,
                        true_obliquity_nutation, atol=0.1))

    def test_mean_obliquity(self):
        """Test `AstronomicalQuantities.mean_obliquity`.

        Notes
        -----
        Test cases are taken from PHP Science Labs.[1]_

        References
        ----------
        .. [1] Jay Tanner. Obliquity of the Ecliptic - PHP Science Labs. 2021.
           url: https://www.neoprogrammics.com/obliquity_of_the_ecliptic/Obliq
           uity_Of_The_Ecliptic_Calculator.php
        """

        julian = np.array([2459437.81600, 2477404.57292, 2422327.21875])

        true_mean_obliquity = np.array([23.43648, 23.43008, 23.44969])
        test_mean_obliquity = mean_obliquity(julian)

        self.assertTrue(np.allclose(test_mean_obliquity, true_mean_obliquity,
                        atol=0.0001))

    def test_apparent_obliquity(self):
        """Test `AstronomicalQuantities.apparent_obliquity`.

        Notes
        -----
        Test cases are taken from PHP Science Labs.[1]_

        References
        ----------
        .. [1] Jay Tanner. Obliquity of the Ecliptic - PHP Science Labs. 2021.
           url: https://www.neoprogrammics.com/obliquity_of_the_ecliptic/Obliq
           uity_Of_The_Ecliptic_Calculator.php
        """

        julian = np.array([2459437.81597, 2477404.57292, 2422327.21875])

        true_apparent_obliquity = np.array([23.43763, 23.42763, 23.44798])
        calc_apparent_obliquity = apparent_obliquity(julian)

        self.assertTrue(np.allclose(calc_apparent_obliquity,
                        true_apparent_obliquity, atol=0.0001))

    def test_from_julian(self):
        """Test `AstronomicalQuantities.from_julian`.

        Notes
        -----
        Test cases are generated from the `Julian` Python library.
        """

        julian = np.array([2436116.31000, 2445246.65000, 2456124.09000])
        test_year, test_month, test_day = from_julian(julian)

        true_year, true_month, true_day = [], [], []

        for current_julian in julian:
            dt = jd.from_jd(current_julian)

            seconds = dt.second + dt.microsecond / 1e6
            minutes = dt.minute + seconds / 60
            hours = dt.hour + minutes / 60
            day = dt.day + hours / 24
            month = dt.month
            year = dt.year

            true_year.append(year)
            true_month.append(month)
            true_day.append(day)

        self.assertTrue(np.array_equal(test_year, true_year))
        self.assertTrue(np.array_equal(test_month, true_month))
        self.assertTrue(np.allclose(test_day, true_day, atol=0.0001))

    def test_day_of_year(self):
        """Test `AstronomicalQuantities.day_of_year`.

        Notes
        -----
        Test cases are taken from "Astronomical Algorithms" by Jean Meeus.[1]_

        References
        ----------
        .. [1] Jean Meeus. Astronomical algorithms. 2nd ed. Willmann-Bell,
           1998, pp. 65. isbn: 9780943396613.
        """

        julian = np.array([2443826.5, 2447273.5, 2447273.8, 2447274.4])

        true_day = np.array([318, 113, 113, 113])
        test_day = day_of_year(julian)

        self.assertTrue(np.array_equal(test_day, true_day))

    def test_equation_of_time(self):
        """Test `AstronomicalQuantities.equation_of_time`.

        Notes
        -----
        Test cases are taken from "Astronomical Algorithms" by Jean Meeus,
        PLANETCALC, and the Global Monitoring Laboratory.[1]_[2]_[3]_

        References
        ----------
        .. [1] Jean Meeus. Astronomical algorithms. 2nd ed. Willmann-Bell,
           1998, pp. 185. isbn: 9780943396613.
        .. [2] Anton. Online calculator: Equation of time. url:
           https://planetcalc.com/9235/.
        .. [3] NOAA US Department of Commerce. ESRL Global Monitoring
           Laboratory -Global Radiation and Aerosols. url:
           https://gml.noaa.gov/grad/solcalc/.
        """

        julian = np.array([2455368.75, 2448908.50, 2459448.50])

        true_equation_of_time = np.array([-0.42658, 3.42735, -0.71054])
        test_equation_of_time = equation_of_time(julian)

        self.assertTrue(np.allclose(test_equation_of_time,
                        true_equation_of_time, atol=0.04))

    def test_equation_of_equinoxes(self):
        """Test `AstronomicalQuantities.equation_of_equinoxes`.

        Notes
        -----
        Test cases taken from "Astronomical Algorithms" by Jean Meeus.[1]_

        References
        ----------
        .. [1] Jean Meeus. Astronomical algorithms. 2nd ed. Willmann-Bell,
           1998, pp. 88. isbn: 9780943396613.
        """

        julian = np.array([2446895.5])

        true_equation_of_equinoxes = np.array([-0.23170])
        test_equation_of_equinoxes = equation_of_equinoxes(julian)

        self.assertTrue(np.allclose(test_equation_of_equinoxes,
                        true_equation_of_equinoxes, atol=0.005))

    def test_sun_right_ascension(self):
        """Test `AstronomicalQuantities.sun_right_ascension`.

        Notes
        -----
        The test case is taken from "Astronomical Algorithms" by Jean
        Meeus.[1]_

        References
        ----------
        .. [1] Jean Meeus. Astronomical algorithms. 2nd ed. Willmann-Bell,
           1998, pp. 185. isbn: 9780943396613.
        """

        julian = np.array([2448908.5])

        true_sun_right_ascension = np.array([198.38083])
        test_sun_right_ascension = sun_right_ascension(julian)

        self.assertTrue(np.allclose(true_sun_right_ascension,
                        test_sun_right_ascension, atol=0.001))


if __name__ == "__main__":
    unittest.main()
