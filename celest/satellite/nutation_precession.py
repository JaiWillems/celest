"""Nutation and precession calculations for coordinate conversions.

This module contains the functionality to incorporate nutation and precession
effects into GCRS and ITRS coordinate conversions.
"""


from celest.satellite._astronomical_quantities import mean_obliquity
from typing import Tuple
import numpy as np


def nutation_angles(julian: np.ndarray) -> Tuple:
    """Return the five Earth nutation angles.

    This method takes a Julian day array and calculates the five nutation
    angles at each time. The calculated angles include the mean elongation of
    the Moon from the Sun (D), mean anomaly of the Sun (M), mean anomaly of
    the Moon (N), Moon's argument of latitude (F), and the longitude of the
    ascending node of the Moon's mean orbit on the ecliptic measured from the
    mean equinox date (O).

    Parameters
    ----------
    jullian : np.ndarray
        Array of shape (n,) containing Julian times in decimal days.

    Returns
    -------
    tuple
        Returns a tuple, `(D, M, N, F, O)`, where the items are Numpy arrays
        of shape (n,) containing the Earth nutation angles in degrees and
        decimals.

    Notes
    -----
    The time in Julian centeries since J2000, :math:`T`, can be calculated
    from the Julian day, :math:`JD`, from the following:

    .. math:: T = \\frac{JD - 2451545}{36525}

    The nutation angles can then be calculated in degrees and decimals. [Mee98b]_

    .. math:: D = 297.85036 + 445267.111480T - 0.0019142T^2 + T^3 / 189474

    .. math:: M = 357.52772 + 35999.050340T - 0.0001603T^2 - T^3 / 300000

    .. math:: N = 134.96298 + 477198.867398T + 0.0086972T^2 + T^3 / 56250

    .. math:: F = 93.27191 + 483202.017538T - 0.0036825T^2 + T^3 / 327270

    .. math:: O = 125.04452 - 1934.136261T + 0.0020708T^2 + T^3 / 450000

    References
    ----------
    .. [Mee98b] Jean Meeus. Astronomical algorithms. 2nd ed. Willmann-Bell,
       1998, pp. 143-144. isbn: 9780943396613.

    Examples
    --------
    >>> julian = np.array([2446895.5])
    >>> D, M, N, F, O = nutation_angles(julian=julian)
    """

    T = (julian - 2451545) / 36525
    T2, T3 = T ** 2, T ** 3

    D = (297.85036 + 445267.111480 * T - 0.0019142 * T2 + T3 / 189474) % 360
    M = (357.52772 + 35999.050340 * T - 0.0001603 * T2 - T3 / 300000) % 360
    N = (134.96298 + 477198.867398 * T + 0.0086972 * T2 + T3 / 56250) % 360
    F = (93.27191 + 483202.017538 * T - 0.0036825 * T2 + T3 / 327270) % 360
    O = (125.04452 - 1934.136261 * T + 0.0020708 * T2 + T3 / 450000) % 360

    return D, M, N, F, O


def nutation_components(julian: np.ndarray) -> Tuple:
    """Return the nutations in longitude and obliquity.

    Parameters
    ----------
    julian : np.ndarray
        Array of shape (n,) containing Julian times in decimal days.

    Returns
    -------
    tuple
        Returns a tuple, `(longitude, obliquity)`, where the items are NumPy
        arrays of shape (n,) containing the nutation components of longitude
        and obliquity in seconds and decimals.

    Notes
    -----
    The time in Julian centeries since J2000, :math:`T`, can be calculated
    from the Julian day, :math:`JD`, from the following:

    .. math:: T = \\frac{JD - 2451545}{36525}

    The longitude of the ascending node of the Moon's orbit on the ecliptic
    measured from the mean equinox can be calculated in decimal degrees
    using the following:

    .. math:: \Omega = 125.04452 - 1934.136261T + 0.0020708T^2 + T^3 / 450000

    The mean longitudes of the Sun, :math:`L`, and Moon, :math:`L'`, can be
    calculated in decimal degrees using the following:

    .. math:: L = 280.4665 + 36000.7698 * T

    .. math:: L' = 218.3165 + 481267.8813 * T

    The nutation in longitude, :math:`\Delta\psi`, and nutation in obliquity,
    :math:`\Delta\epsilon`, can then be calculated in decimal minutes. [Mee98c]_

    .. math:: \Delta\Psi = -17.20\sin\Omega - 1.32\sin 2L - 0.23\sin 2L' + 0.21\sin 2\Omega

    .. math:: \Delta\epsilon = 9.20\cos\Omega + 0.57\cos 2L + 0.10\cos 2L' - 0.09\cos 2\Omega

    References
    ----------
    .. [Mee98c] Jean Meeus. Astronomical algorithms. 2nd ed. Willmann-Bell,
       1998, pp. 143-144. isbn: 9780943396613.

    Examples
    --------
    >>> julian = np.array([2446895.5])
    >>> delta_psi, delta_epsilon = nutation_componenets(julian=julian)
    """

    T = (julian - 2451545) / 36525

    omega = 125.04452 - 1934.136261 * T + 0.0020708 * T ** 2 + T ** 3 / 450000
    omega = np.radians(omega)
    omega_2 = 2 * omega

    L = 2 * np.radians(280.4665 + 36000.7698 * T)
    LP = 2 * np.radians(218.3165 + 481267.8813 * T)

    P1 = -17.20 * np.sin(omega)
    P2 = -1.32 * np.sin(L)
    P3 = -0.23 * np.sin(LP)
    P4 = 0.21 * np.sin(omega_2)

    E1 = 9.20 * np.cos(omega)
    E2 = 0.57 * np.cos(L)
    E3 = 0.10 * np.cos(LP)
    E4 = -0.09 * np.cos(omega_2)

    delta_psi = P1 + P2 + P3 + P4
    delta_epsilon = E1 + E2 + E3 + E4

    return delta_psi, delta_epsilon


def precession_angles(julian: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return precession angles.

    Parameters
    ----------
    julian : np.ndarray
        Array of shape (n,) containing Julian time in the J2000 epoch.

    Returns
    -------
    np.ndarray
        Tuple of length 3 containing the precession angles in arc seconds and
        decimals.

    Notes
    -----
    The convetional precession angles are those derived using the IAU 2000A
    model.[SL13c]_

    References
    ----------
    .. [SL13c] M. Soffel and R. Langhans. Space-Time Reference Systems.
       Astronomy and Astrophysics Library. Springer-Verlag, 2013, pp. 219.
    """

    # Get time constants.
    t = (julian - 2451545.0) / 36525
    t2 = t * t
    t3 = t2 * t
    t4 = t3 * t
    t5 = t4 * t

    # Get precession angles.
    zeta = 2.59796176 + 2306.0809506 * t + 0.3019015 * t2 + 0.0179663 * t3 - \
        0.0000327 * t4 - 0.0000002 * t5
    theta = 2004.1917476 * t - 0.4269353 * t2 - 0.0418251 * t3 - \
        0.0000601 * t4 - 0.0000001 * t5
    z = - 2.5976176 + 2306.0803226 * t + 1.0947790 * t2 + 0.0182273 * t3 + \
        0.0000470 * t4 - 0.0000003 * t5

    return zeta, theta, z


def bias_matrix() -> np.ndarray:
    r"""Generate bias matrix for GCRS and ITRS conversions.

    Returns
    -------
    np.ndarray
        Array of shape (3,3) representing the bias matrix.

    Notes
    -----
    The bias matrixs is a 3-2-1 set of Euler angle rotations about
    :math:`d\\alpha_0=-0.01460`, :math:`\xi_0=-0.0166170`, and
    :math:`-\eta_0=0.0068192`.[SL13c]_

    References
    ----------
    .. [SL13c] M. Soffel and R. Langhans. Space-Time Reference Systems.
       Astronomy and Astrophysics Library. Springer-Verlag, 2013, pp. 197-233.
    """

    xi = -0.0166170
    eta = -0.0068192
    dalpha = -0.01460

    ang1 = np.radians(dalpha / 3600)
    ang2 = np.radians(xi / 3600)
    ang3 = np.radians(-eta / 3600)

    s1, c1 = np.sin(ang1), np.cos(ang1)
    s2, c2 = np.sin(ang2), np.cos(ang2)
    s3, c3 = np.sin(ang3), np.cos(ang3)

    matrix = np.zeros((3, 3))
    matrix[0, 0] = c1 * c2
    matrix[0, 1] = s1 * c2
    matrix[0, 2] = - s2
    matrix[1, 0] = - s1 * c3 + c1 * s2 * s3
    matrix[1, 1] = c1 * c3 + s1 * s2 * s3
    matrix[1, 2] = c2 * s3
    matrix[2, 0] = s1 * s3 + c1 * s2 * c3
    matrix[2, 1] = - c1 * s3 + s1 * s2 * c3
    matrix[2, 2] = c2 * c3

    return matrix


def precession_matrix(julian: np.ndarray) -> np.ndarray:
    """Return precession tensor.

    Parameters
    ----------
    julian : np.ndarray
        Array of shape (n,) containing Julian times in the J2000 epoch.

    Returns
    -------
    np.ndarray
        Array of shape (3,3,n) containing a series of 3x3 precession matrices;
        one for each input time.

    Notes
    -----
    The precession calculations use the methodology put forward in "Space-Time
    Reference Systems".[SL13c]_

    References
    ----------
    .. [SL13c] M. Soffel and R. Langhans. Space-Time Reference Systems.
       Astronomy and Astrophysics Library. Springer-Verlag, 2013, pp. 197-233.
    """

    zeta, theta, z = precession_angles(julian=julian)

    ang1 = np.radians(zeta / 3600)
    ang2 = np.radians(theta / 3600)
    ang3 = - np.radians(z / 3600)

    s1, c1 = np.sin(ang1), np.cos(ang1)
    s2, c2 = np.sin(ang2), np.cos(ang2)
    s3, c3 = np.sin(ang3), np.cos(ang3)

    # Construct matrix.
    matrix = np.zeros((3, 3, julian.size))
    matrix[0, 0, :] = - s1 * s3 + c1 * c2 * c3
    matrix[0, 1, :] = c1 * s3 + s1 * c2 * c3
    matrix[0, 2, :] = - s2 * c3
    matrix[1, 0, :] = - s1 * c3 - c1 * c2 * s3
    matrix[1, 1, :] = c1 * c3 - s1 * c2 * s3
    matrix[1, 2, :] = s2 * s3
    matrix[2, 0, :] = c1 * s2
    matrix[2, 1, :] = s1 * s2
    matrix[2, 2, :] = c2

    return matrix


def nutation_matrix(julian: np.ndarray) -> np.ndarray:
    """Return nutation tensor.

    Parameters
    ----------
    julian : np.ndarray
        Array of shape (n,) containing Julian times in the J2000 epoch.

    Returns
    -------
    np.ndarray
        Array of shape (3,3,n) containing a series of 3x3 nutation matrices;
        one for each input time.

    Notes
    -----
    The Nutation calculations use the methodology put forward in "Space-Time
    Reference Systems".[SL13c]_

    References
    ----------
    .. [SL13c] M. Soffel and R. Langhans. Space-Time Reference Systems.
       Astronomy and Astrophysics Library. Springer-Verlag, 2013, pp. 197-233.
    """

    # Get time constants.
    t = (julian - 2451545.0) / 36525
    t2 = t * t
    t3 = t2 * t

    # Get angles.
    delta_psi, delta_eps = nutation_components(julian)
    eps_A = 3600 * mean_obliquity(julian) - 46.84024 * t - 0.00059 * t2 + \
        0.001813 * t3

    ang1 = np.radians(eps_A / 3600)
    ang2 = - np.radians(delta_psi / 3600)
    ang3 = - np.radians(eps_A / 3600 + delta_eps / 3600)

    s1, c1 = np.sin(ang1), np.cos(ang1)
    s2, c2 = np.sin(ang2), np.cos(ang2)
    s3, c3 = np.sin(ang3), np.cos(ang3)

    # Construct matrix.
    matrix = np.zeros((3, 3, julian.size))
    matrix[0, 0, :] = c2
    matrix[0, 1, :] = s1 * s2
    matrix[0, 2, :] = - c1 * s2
    matrix[1, 0, :] = s2 * s3
    matrix[1, 1, :] = c1 * c3 - s1 * c2 * s3
    matrix[1, 2, :] = s1 * c3 + c1 * c2 * s3
    matrix[2, 0, :] = s2 * c3
    matrix[2, 1, :] = - c1 * s3 - s1 * c2 * c3
    matrix[2, 2, :] = - s1 * s3 + c1 * c2 * c3

    return matrix
