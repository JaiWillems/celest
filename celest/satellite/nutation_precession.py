"""Nutation and precession calculations for coordinate conversions.

This module contains the functionality to incorporate nutation and precession
effects into GCRS and ITRS coordinate conversions.
"""


from celest.satellite._astronomical_quantities import (
    mean_obliquity, nutation_components, precession_angles
)
from typing import Tuple
import numpy as np


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
