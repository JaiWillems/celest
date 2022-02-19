

from typing import Literal
import numpy as np


def _cone_constraint(theta: np.ndarray, U: np.ndarray, X: np.ndarray) -> np.ndarray:
    """Return indices for `X` points falling within a cone.

    The cone geometry is defined with the apex located at `U`, an axis parallel
    to `U` and an aperture angle `theta`.

    Parameters
    ----------
    theta : np.ndarray
        1-D array containing cone aperture angles in decimal degrees.
    U : np.ndarray
        1-D array containing cartesian coordinates of the apex.
    X : np.ndarray
        2-D array containing cartesian coordinates to check.

    Returns
    -------
    np.ndarray
        1-D array of indices corresponding to points that fall within the cone.
    """

    U_norm = np.linalg.norm(U, axis=1)

    d1, d2, d3 = (U / U_norm).T
    y1, y2, y3 = (X - U).T

    ct2 = np.cos(np.radians(theta)) ** 2

    r1 = (d1 ** 2 - ct2) * y1 + d1 * d2 * y2 + d1 * d3 * y3
    r2 = d2 * d1 * y1 + (d2 ** 2 - ct2) * y2 + d2 * d3 * y3
    r3 = d3 * d1 * y1 + d3 * d2 * y2 + (d3 ** 2 - ct2) * y3

    point = y1 * r1 + y2 * r2 + y3 * r3

    ind = np.where(point >= 0)[0]

    return ind


def _plane_constraint(U: np.ndarray, X: np.ndarray) -> np.ndarray:
    """Return indices for `X` points falling above the plane.

    The plane geometry is defined as the plane centered at and normal to the
    `U` vector. 

    Parameters
    ----------
    U : np.ndarray
        1-D array containing cartesian coordinates of the plane normal.
    X : np.ndarray
        2-D array containing cartesian coordinates to check.

    Returns
    -------
    np.ndarray
        1-D array of indices corresponding to points that fall above the plane.

    Notes
    -----
    A point :math:`\textbf{X}=\langle x, y, z\rangle^T` falls falls above the
    plane defined with :math:`\textbf{U}=\langle x_0, y_0, z_0\rangle^T` if it
    satisfies the following inequality:

    .. math:: 0\leq\left(\textbf{X}-\textbf{U}\right)^T\cdot\textbf{U}
    """

    u1, u2, u3 = U.T
    y1, y2, y3 = (X - U).T

    point = y1 * u1 + y2 * u2 + y3 * u3

    ind = np.where(point >= 0)[0]

    return ind


def _aperature_theta(ang: float, form: Literal[0, 1], U: np.ndarray=None,
                     X: np.ndarray=None) -> Literal[np.ndarray, float]:
    """Calculate cone aperature angles from constraint angles.

    The cone aperature angles are linked to the constraint angles to define a
    quasi-valid encounter region.

    Parameters
    ----------
    ang : float
        The constraint angle in degrees.
    form : {0, 1}
        Specifies `ang` as either altitude (`angType=0`) or off-nadir
        (`angType=1`).
    U : np.ndarray, optional
        1-D array containing cartesian coordinates of the apex.
    X : np.ndarray
        2-D array containing cartesian coordinates to check.

    Returns
    -------
    np.ndarray or float
        If the altitude angle is selected, the output will be a scalar value.
        For off-nadir constraint angles, the output will be a 1-D array
        containing the aperture angles in degrees.

    Notes
    -----
    If the constraint angle type is an altitude angle, then :math:`\theta` can
    be calculated from :math:`\theta=90\degree-a\degree` where :math:`a` is the
    altitude constraint angle.

    If the constraint angle type is an off-nadir angle, then :math:`theta` can
    be calculated from:

    .. math::

        \theta = \gamma+\sin^{-1}\left(\frac{\left|\textbf{X}-\textbf{U}
        \right|}{\left|\textbf{U}\right|}\sin{\gamma}\right)

    where :math:`\gamma` is the off-nadir constraint angle.
    """

    if form:
        ang = np.radians(ang)

        num = np.linalg.norm(X - U, axis=1)
        denom = np.linalg.norm(U)

        theta = ang + np.arcsin(np.sin(ang) * num / denom)
        theta = np.degrees(theta)
    else:
        theta = 90 - ang

    return theta


def _analytical_encounter_ind(sat_position: np.ndarray, gnd_position: np.ndarray,
                              ang: float, form: Literal[0, 1]) -> np.ndarray:
    """Return encounter indices.

    The encounter indices are the position/time indices of the satellite where
    it satisfies both the cone and plane geometry constraints.

    Parameters
    ----------
    sat_position : np.ndarray
        2-D array containing x, y, z itrs coordinates of the satellite.
    gnd_position : np.ndarray
        1-D array containing the x, y, z itrs coordinates of the ground
        location.
    ang : float
        Constraint angle in degrees.
    angType : {0, 1}
        Specifies `ang` as either altitude (`angType=0`) or off-nadir
        (`angType=1`).

    Returns
    -------
    np.ndarray
        Array of indices defining valid encounter positions.
    """

    U, X = gnd_position, sat_position
    theta = _aperature_theta(ang, form, U, X)

    ind_1 = _cone_constraint(theta, U, X)
    ind_2 = _plane_constraint(U, X)

    ind_final = np.intersect1d(ind_1, ind_2)

    return ind_final
