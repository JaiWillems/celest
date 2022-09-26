

from celest.coordinates.astronomical_quantities import (
    _calculate_raw_elapsed_jd_century_powers,
    mean_obliquity,
    nutation_components,
    conventional_precession_angles
)
from celest.units.quantity import Quantity
from celest import units as u
import numpy as np


def bias_matrix() -> np.ndarray:
    r"""Generate bias matrix for GCRS and ITRS conversions.

    Returns
    -------
    np.ndarray
        2-D array representing the bias matrix.

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

    xi, eta, dalpha = -0.0166170, -0.0068192, -0.01460

    a1 = np.radians(dalpha / 3600)
    a2 = np.radians(xi / 3600)
    a3 = np.radians(-eta / 3600)

    s1, c1 = np.sin(a1), np.cos(a1)
    s2, c2 = np.sin(a2), np.cos(a2)
    s3, c3 = np.sin(a3), np.cos(a3)

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


def precession_matrix(julian: Quantity) -> np.ndarray:
    """Return precession tensor.

    Parameters
    ----------
    julian : Quantity
        Quantity object containing date data in the JD2000 epoch.

    Returns
    -------
    np.ndarray
        3-D array of shape (n,3,3) containing 3x3 precession matrices
        associated with each time along the 0 axis.

    Notes
    -----
    The precession matrix is implemented as defined in "Space-Time Reference
    Systems". [SL13c]_

    References
    ----------
    .. [SL13c] M. Soffel and R. Langhans. Space-Time Reference Systems.
       Astronomy and Astrophysics Library. Springer-Verlag, 2013, pp. 197-233.
    """

    zeta, theta, z = conventional_precession_angles(julian)

    a1 = zeta.to(u.rad)
    a2 = theta.to(u.rad)
    a3 = - z.to(u.rad)

    s1, c1 = np.sin(a1), np.cos(a1)
    s2, c2 = np.sin(a2), np.cos(a2)
    s3, c3 = np.sin(a3), np.cos(a3)

    matrix = np.zeros((julian.data.size, 3, 3))
    matrix[:, 0, 0] = - s1 * s3 + c1 * c2 * c3
    matrix[:, 0, 1] = c1 * s3 + s1 * c2 * c3
    matrix[:, 0, 2] = - s2 * c3
    matrix[:, 1, 0] = - s1 * c3 - c1 * c2 * s3
    matrix[:, 1, 1] = c1 * c3 - s1 * c2 * s3
    matrix[:, 1, 2] = s2 * s3
    matrix[:, 2, 0] = c1 * s2
    matrix[:, 2, 1] = s1 * s2
    matrix[:, 2, 2] = c2

    return matrix


def nutation_matrix(julian: Quantity) -> np.ndarray:
    """Return nutation tensor.

    Parameters
    ----------
    julian : Quantity
        Quantity object containing date data in the JD2000 epoch.

    Returns
    -------
    np.ndarray
        3-D array of shape (n,3,3) containing 3x3 nutation matrices associated
        with each time along the 0 axis.

    Notes
    -----
    The nutation matrix is implemented as defined in "Space-Time Reference
    Systems". [SL13c]_

    References
    ----------
    .. [SL13c] M. Soffel and R. Langhans. Space-Time Reference Systems.
       Astronomy and Astrophysics Library. Springer-Verlag, 2013, pp. 197-233.
    """

    nutation_in_longitude, nutation_in_obliquity = nutation_components(julian)
    t1, t2, t3 = _calculate_raw_elapsed_jd_century_powers(julian, 3)
    eps_A = mean_obliquity(julian).to(u.arcsec) - 46.84024 * t1 - 0.00059 * t2 + 0.001813 * t3

    a1 = np.radians(eps_A / 3600)
    a2 = - nutation_in_longitude.to(u.rad)
    a3 = - np.radians(eps_A / 3600 + nutation_in_obliquity.to(u.deg))

    s1, c1 = np.sin(a1), np.cos(a1)
    s2, c2 = np.sin(a2), np.cos(a2)
    s3, c3 = np.sin(a3), np.cos(a3)

    matrix = np.zeros((julian.data.size, 3, 3))
    matrix[:, 0, 0] = c2
    matrix[:, 0, 1] = s1 * s2
    matrix[:, 0, 2] = - c1 * s2
    matrix[:, 1, 0] = s2 * s3
    matrix[:, 1, 1] = c1 * c3 - s1 * c2 * s3
    matrix[:, 1, 2] = s1 * c3 + c1 * c2 * s3
    matrix[:, 2, 0] = s2 * c3
    matrix[:, 2, 1] = - c1 * s3 - s1 * c2 * c3
    matrix[:, 2, 2] = - s1 * s3 + c1 * c2 * c3

    return matrix
