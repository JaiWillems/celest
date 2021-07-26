"""Astronomical quantities for astronomical calculations.

The _astronomicalquantities module contains the `AstronomicalQuantities`
class to allow for more advanced astronomical calculations.
"""


from typing import Tuple
import numpy as np
import julian


class AstronomicalQuantities(object):
    """Astronomical quantities for astronomical calculations.

    The `AstronomicalQuantities` class to allow for more advanced astronomical
    calculations.

    Methods
    -------
    nutation_angles(julData)
        Return the five Earth nutation angles.
    nutation_components(julData)
        Return the nutations in longitude and in obliquity.
    mean_obliquity(julData)
        Return the mean obliquity of the ecliptic.
    apparent_obliquity(julData)
        Return the apparent obliquity of the ecliptic.
    equation_of_time(julData)
        Return the Equation of Time.
    equation_of_equinoxes(julData)
        Return the equation of equinoxes.
    """

    def nutation_angles(self, julData: np.array) -> Tuple:
        """Return the five Earth nutation angles.

        Parameters
        ----------
        julData : np.array
            Times for the nutation angle calculations.

        Returns
        -------
        tuple
            Returns a tuple of Numpy arrays of Earth nutation angles in the
            order of: mean elongation of the Moon from the sun, mean anomoly
            of the Sun, mean anomoly of the Moon, the Moon's argument of
            lattitude, and the longitude of the ascending node of the Moon's
            mean orbit on the ecliptic.
        """

        T = (julData - 2451545) / 36525

        D = (297.85036 + 445267.111480 * T - 0.0019142 * T ** 2 + T ** 3 / 189474) % 360
        M = (357.52772 + 35999.050340 * T - 0.0001603 * T ** 2 - T ** 3 / 300000) % 360
        N = (134.96298 + 477198.867398 * T + 0.0086972 * T ** 2 + T ** 3 / 56250) % 360
        F = (93.27191 + 483202.017538 * T - 0.0036825 * T ** 2 + T ** 3 / 327270) % 360
        O = (125.04452 - 1934.136261 * T + 0.0020708 * T ** 2 + T ** 3 / 450000) % 360

        return D, M, N, F, O

    def nutation_components(self, julData: np.array) -> Tuple:
        """Return the nutations in longitude and in obliquity.

        Parameters
        ----------
        julData : np.array
            Times for the nutation component calculations.

        Returns
        -------
        tuple
            Returns a tuple of arrays containing the nutation of longitude and
            nutation of obliquity.
        """

        T = (julData - 2451545) / 36525

        omega = np.radians(125.04452 - 1934.136261 * T + 0.0020708 * T ** 2 + T ** 3 / 450000) % 360)

        L = np.radians(280.4665 + 36000.7698 * T)
        Lp = np.radians(218.3165 + 481267.8813 * T)

        p1 = -17.20 * np.sin(omega)
        p2 = -1.32 * np.sin(2 * L)
        p3 = -0.23 * np.sin(2 * Lp)
        p4 = 0.21 * np.sin(2 * omega)

        e1 = 9.20 * np.cos(omega)
        e2 = 0.57 * np.cos(2 * L)
        e3 = 0.10 * np.cos(2 * Lp)
        e4 = -0.09 * np.cos(2 * omega)

        delta_psi = p1 + p2 +p3 + p4
        delta_epsilon = e1 + e2 + e3 + e4

        return delta_psi, delta_epsilon

    def mean_obliquity(self, julData: np.array) -> np.array:
        """Return the mean obliquity of the ecliptic.

        Parameters
        ----------
        julData : np.array
            Times for the obliquity angle calculations.

        Returns
        -------
        np.array
            Returns an array of shape (n,) containing the mean obliquity of
            the ecliptic.
        """

        T = (julData - 2451545) / 36525
        U = T / 100

        e1 = 84381.448
        e2 = -4680.93 * U
        e3 = -1.55 * U ** 2
        e4 = 1999.25 * U ** 3
        e5 = -51.38 * U ** 4
        e6 = -249.67 * U **5
        e7 = -39.05 * U ** 6
        e8 = 7.12 * U ** 7
        e9 = 27.87 * U **8
        e10 = 5.79 * U ** 9
        e11 = 2.45 * U ** 10

        epsilon_0 = (e1 + e2 + e3 + e4 + e5 + e6 + e7 + e8 + e9 + e10 + e11) / 3600

        return epsilon_0

    def apparent_obliquity(self, julData: np.array) -> np.array:
        """Return the apparent obliquity of the ecliptic.

        Parameters
        ----------
        julData : np.array
            Times for the obliquity angle calculations.

        Returns
        -------
        np.array
            Returns an array of shape (n,) containing the apparent obliquity of
            the ecliptic.
        """

        epsilon_0 = self.mean_obliquity(julData)
        _, delta_epsilon = self.nutation_components(julData)
        epsilon = epsilon_0 + delta_epsilon / 3600

        return epsilon

    def equation_of_time(self, julData: np.array) -> np.array:
        """Return the Equation of Time.

        Parameters
        ----------
        kwargs : Dict, optional
            Optional keyword arguments passed into the `Interpolation()`
            callable.

        Returns
        -------
        np.array
            Array of shape (n,) containing the Equation of Time in decimal
            minutes.

        Notes
        -----
        This method uses an empirical equation to approximate the Equation of
        Time within an accuracy of 0.5 minutes.

        .. math:: EoT = 9.87\sin(2B) - 7.53\cos(B) - 1.5\sin(B)

        where :math:`B = \\frac{360}{365}(d - 81)` and :math:`b` is the number
        of days since the start of the year.
        """

        day_num = np.zeros((self.length,))
        for i, time in enumerate(self._julian):
            day_num[i] = float(julian.from_jd(time).strftime("%j"))

        B = (day_num - 81) * 360 / 365
        B_r, B2_r = np.radians(B), np.radians(2 * B)
        EoT = 9.87 * np.sin(B2_r) - 7.53 * np.cos(B_r) - 1.5 * np.sin(B_r)

        return EoT

    def equation_of_equinoxes(self, julData: np.array) -> np.array:
        """Return the equation of equinoxes.

        Parameters
        ----------
        julData : np.array
            Times for the equation of equinox calculations.

        Returns
        -------
        np.array
            Returns an array of shape (n,) containing equation of equinoxes.
        """

        delta_psi, _ = self.nutation_components(julData)
        epsilon = self.apparent_obliquity(julData)
        EoE = delta_psi * np.cos(np.radians(epsilon)) / 15

        return EoE
