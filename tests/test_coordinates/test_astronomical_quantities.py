

from celest.units.quantity import Quantity
from celest.coordinates.astronomical_quantities import (
    earth_rotation_angle,
    nutation_angles,
    nutation_components,
    mean_obliquity,
    apparent_obliquity,
    from_julian,
    day_of_year,
    equation_of_time,
    equation_of_equinoxes,
    sun_right_ascension,
    radius_of_curvature_in_prime_vertical
)
from celest import units as u
from unittest import TestCase


class TestAstronomicalQuantities(TestCase):

    def test_earth_rotation_for_validation(self):
        julian = Quantity(2454545, u.jd2000)
        self.assertAlmostEqual(6.2360075, earth_rotation_angle(julian).to(u.rad),
                               delta=0.01)

    def test_nutation_angles_for_validation(self):
        """
        Notes
        -----
        Test case is taken from Astronomical Algorithms.[1]_

        References
        ----------
        .. [1] Jean Meeus. Astronomical algorithms. 2nd ed. Willmann-Bell,
           1998, pp. 148. isbn: 9780943396613.
        """

        julian = Quantity(2446895.5, u.jd2000)
        d, m, n, f, o = nutation_angles(julian)

        self.assertAlmostEqual(d.to(u.deg), 136.9623, places=4)
        self.assertAlmostEqual(m.to(u.deg), 94.9792, places=4)
        self.assertAlmostEqual(n.to(u.deg), 229.2784, places=4)
        self.assertAlmostEqual(f.to(u.deg), 143.4079, places=4)
        self.assertAlmostEqual(o.to(u.deg), 11.2531, places=4)

    def test_nutation_components_for_validation(self):
        """
        Notes
        -----
        Test cases are taken from the PHP Science Labs.[1]_

        References
        ----------
        .. [1] Jay Tanner. NeoProgrammics - Science Computations. 2021.
           url: http://www.neoprogrammics.com/nutations/.
        """

        julian = Quantity(2449634.50000, u.jd2000)
        nutation_in_longitude, nutation_in_obliquity = nutation_components(julian)

        self.assertAlmostEqual(11.694, nutation_in_longitude.data, delta=0.5)
        self.assertAlmostEqual(-5.946, nutation_in_obliquity.data, delta=0.1)

    def test_mean_obliquity_for_validation(self):
        """
        Notes
        -----
        Test cases are taken from PHP Science Labs.[1]_

        References
        ----------
        .. [1] Jay Tanner. Obliquity of the Ecliptic - PHP Science Labs. 2021.
           url: https://www.neoprogrammics.com/obliquity_of_the_ecliptic/Obliq
           uity_Of_The_Ecliptic_Calculator.php
        """

        julian = Quantity(2459437.81600, u.jd2000)
        self.assertAlmostEqual(23.43648, mean_obliquity(julian).to(u.deg),
                               delta=0.0001)

    def test_apparent_obliquity_for_validation(self):
        """
        Notes
        -----
        Test cases are taken from PHP Science Labs.[1]_

        References
        ----------
        .. [1] Jay Tanner. Obliquity of the Ecliptic - PHP Science Labs. 2021.
           url: https://www.neoprogrammics.com/obliquity_of_the_ecliptic/Obliq
           uity_Of_The_Ecliptic_Calculator.php
        """

        julian = Quantity(2459437.81597, u.jd2000)
        self.assertAlmostEqual(23.43763, apparent_obliquity(julian).to(u.deg),
                               delta=0.0001)

    def test_from_julian_for_validation(self):
        """
        Notes
        -----
        Test cases are generated from the `Julian` Python library.
        """

        julian = Quantity(2436116.31000, u.jd2000)
        year, month, day = from_julian(julian)

        self.assertTrue(1957, year[0])
        self.assertEqual(10, month[0])
        self.assertAlmostEqual(4.81, day[0], delta=0.0001)

    def test_day_of_year_for_validation(self):
        """
        Notes
        -----
        Test cases are taken from "Astronomical Algorithms" by Jean Meeus.[1]_

        References
        ----------
        .. [1] Jean Meeus. Astronomical algorithms. 2nd ed. Willmann-Bell,
           1998, pp. 65. isbn: 9780943396613.
        """

        julian = Quantity(2443826.5, u.jd2000)
        self.assertEqual(318, day_of_year(julian))

    def test_equation_of_time_for_validation(self):
        """
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

        julian = Quantity(2455368.75, u.jd2000)
        self.assertAlmostEqual(-0.42658, equation_of_time(julian).to(u.deg),
                               delta=0.04)

    def test_equation_of_equinoxes_for_validation(self):
        """
        Notes
        -----
        Test cases taken from "Astronomical Algorithms" by Jean Meeus.[1]_

        References
        ----------
        .. [1] Jean Meeus. Astronomical algorithms. 2nd ed. Willmann-Bell,
           1998, pp. 88. isbn: 9780943396613.
        """

        julian = Quantity(2446895.5, u.jd2000)
        self.assertAlmostEqual(-0.23170, equation_of_equinoxes(julian).to(u.arcsec),
                               delta=0.005)

    def test_sun_right_ascension_for_validation(self):
        """
        Notes
        -----
        The test case is taken from "Astronomical Algorithms" by Jean
        Meeus.[1]_

        References
        ----------
        .. [1] Jean Meeus. Astronomical algorithms. 2nd ed. Willmann-Bell,
           1998, pp. 185. isbn: 9780943396613.
        """

        julian = Quantity(2448908.5, u.jd2000)
        self.assertAlmostEqual(198.38083, sun_right_ascension(julian).to(u.deg),
                               delta=0.001)
