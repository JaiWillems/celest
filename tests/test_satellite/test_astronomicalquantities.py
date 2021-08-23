"""Testing module for the `AstronomicalQuantities` class."""


from celest.satellite import AstronomicalQuantities
from unittest import TestCase
import numpy as np
import unittest


class TestAstronomicalQuantities(TestCase):
    """Testing class for the `AstronomicalQuantities` class.

    Attributes
    ----------
    AO : AstronomicalQuantities
        Set up code for test methods.

    Methods
    -------
    setUp()
        Test fixure for test method execution.
    test_nutation_angles(julData)
        Test `AstronomicalQuantities.nutation_angles`.
    test_nutation_components(julData)
        Test `AstronomicalQuantities.nutation_components`.
    test_mean_obliquity(julData)
        Test `AstronomicalQuantities.mean_obliquity`.
    test_apparent_obliquity(julData)
        Test `AstronomicalQuantities.apparent_obliquity`.
    test_from_julian(julData)
        Test `AstronomicalQuantities.from_julian`.
    test_day_of_year(julData)
        Test `AstronomicalQuantities.day_of_year`.
    test_equation_of_time(julData)
        Test `AstronomicalQuantities.equation_of_time`.
    test_equation_of_equinoxes(julData)
        Test `AstronomicalQuantities.equation_of_equinoxes`.
    """

    def setUp(self):
        """Test fixure for test method execution."""

        self.AO = AstronomicalQuantities()

    def test_nutation_angles(self):
        """Test `AstronomicalQuantities.nutation_angles`.

        Notes
        -----
        Test case taken from Astronomical Algorithms.[1]_

        References
        ----------
        .. [1] Jean Meeus. Astronomical algorithms. 2nd ed. Willmann-Bell,
           1998, pp. 148. isbn: 9780943396613.
        """

        julData = np.array([2446895.5])
        D, M, N, F, O = self.AO.nutation_angles(julData)

        self.assertAlmostEqual(D[0], 136.9623, places=4)
        self.assertAlmostEqual(M[0], 94.9792, places=4)
        self.assertAlmostEqual(N[0], 229.2784, places=4)
        self.assertAlmostEqual(F[0], 143.4079, places=4)
        self.assertAlmostEqual(O[0], 11.2531, places=4)

    def test_nutation_components(self):
        """Test `AstronomicalQuantities.nutation_components`.

        Notes
        -----
        Test cases taken from the PHP Science Labs.[1]_

        References
        ----------
        .. [1] Jay Tanner. NeoProgrammics - Science Computations. 2021.
           url: http://www.neoprogrammics.com/nutations/.
        """

        julData = np.array([2449634.5, 2453420.5625, 2477418.211805555555])
        true_d_psi = np.array([11.694, -5.993, 2.937])
        true_d_epsilon = np.array([-5.946, 8.431, -8.871])

        d_psi, d_epsilon = self.AO.nutation_components(julData)

        for i in range(julData.size):
            with self.subTest(i=i):
                self.assertAlmostEqual(d_psi[i], true_d_psi[i], delta=0.5)
                self.assertAlmostEqual(d_epsilon[i], true_d_epsilon[i], delta=0.1)

    def test_mean_obliquity(self):
        """Test `AstronomicalQuantities.mean_obliquity`.

        Notes
        -----
        Test cases taken from PHP Science Labs.[1]_

        References
        ----------
        .. [1] Jay Tanner. Obliquity of the Ecliptic - PHP Science Labs. 2021.
           url: https://www.neoprogrammics.com/obliquity_of_the_ecliptic/Obliquity_Of_The_Ecliptic_Calculator.php
        """

        julData = np.array([2459437.815972222, 2477404.5729166665, 2422327.21875])
        true_mean_obliquity = np.array([23.4364767133, 23.4300808752, 23.4496874486])

        mean_obliquity = self.AO.mean_obliquity(julData=julData)

        for i in range(julData.size):
            with self.subTest(i=i):
                self.assertAlmostEqual(mean_obliquity[i], true_mean_obliquity[i], delta=0.0001)

    def test_apparent_obliquity(self):
        """Test `AstronomicalQuantities.apparent_obliquity`.

        Notes
        -----
        Test cases taken from PHP Science Labs.[1]_

        References
        ----------
        .. [1] Jay Tanner. Obliquity of the Ecliptic - PHP Science Labs. 2021.
           url: https://www.neoprogrammics.com/obliquity_of_the_ecliptic/Obliquity_Of_The_Ecliptic_Calculator.php
        """

        julData = np.array([2459437.815972222, 2477404.5729166665, 2422327.21875])
        true_apparent_obliquity = np.array([23.4376318857, 23.4276258425, 23.4479812709])

        apparent_obliquity = self.AO.apparent_obliquity(julData=julData)

        for i in range(julData.size):
            with self.subTest(i=i):
                self.assertAlmostEqual(apparent_obliquity[i], true_apparent_obliquity[i], delta=0.0001)

    def test_from_julian(self):
        """Test `AstronomicalQuantities.from_julian`.

        Notes
        -----
        Test cases generated from the `Julian` Python library.
        """

        import julian

        julData = np.array([2436116.31, 2445246.65, 2456124.09])
        year, month, day = self.AO.from_julian(julData)

        for i in range(julData.size):
            with self.subTest(i=i):
                dt = julian.from_jd(julData[i])
                true_day = dt.day + (dt.hour + (dt.minute + (dt.second + dt.microsecond / 100000) / 60) / 60) / 24
                true_month = dt.month
                true_year = dt.year

                self.assertAlmostEqual(day[i], true_day, places=3)
                self.assertEqual(month[i], true_month)
                self.assertEqual(year[i], true_year)

    def test_day_of_year(self):
        """Test `AstronomicalQuantities.day_of_year`.

        Notes
        -----
        Test cases taken from "Astronomical Algorithms" by Jean Meeus.[1]_

        References
        ----------
        .. [1] Jean Meeus. Astronomical algorithms. 2nd ed. Willmann-Bell,
           1998, pp. 65. isbn: 9780943396613.
        """

        julData = np.array([2443826.5, 2447273.5, 2447273.8, 2447274.4])
        day = self.AO.day_of_year(julData)

        true_day = np.array([318, 113, 113, 113])

        for i in range(julData.size):
            with self.subTest(i=i):
                self.assertEqual(day[i], true_day[i])

    def test_equation_of_time(self):
        """Test `AstronomicalQuantities.equation_of_time`.

        Notes
        -----
        Test cases taken from "Astronomical Algorithms" by Jean Meeus,
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

        julData = np.array([2455368.75, 2448908.5, 2459448.5])
        true_EOT = np.array([-0.42657696, 3.427351, -0.710537])

        EOT = self.AO.equation_of_time(julData=julData)

        for i in range(julData.size):
            with self.subTest(i=i):
                self.assertAlmostEqual(EOT[i], true_EOT[i], delta=0.04)

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

        julData = np.array([2446895.5])
        true_EOE = np.array([-0.2317])

        EOE = self.AO.equation_of_equinoxes(julData=julData)

        for i in range(julData.size-1):
            with self.subTest(i=i):
                self.assertAlmostEqual(EOE[i], true_EOE[i], delta=0.0001)


if __name__ == "__main__":
    unittest.main()
