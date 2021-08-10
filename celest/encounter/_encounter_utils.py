"""Encounter indices utilities.

This module contains the geometric utilities neccessary for the `Encounter`
modules analytical pass analysis method.
"""


import numpy as np


def cone_constraint(theta, U, X):
    """Return an array of indices of `X` points that fall within a cone.

    This method defines a double sided cone with an apex located at `U` and an
    aperature angle `theta` and returns the indices of `X` points that fall
    within this geometry.

    Parameters
    ----------
    theta: np.array
        Array of size (n,) containing cone aperature angles in degrees.
    U: np.array
        Array of shape (n, 3) containing the cartesian coordinates of the apex.
    X: np.array
        Array of shape (n, 3) containing the cartesian coordinate points to
        check.

    Returns
    -------
    np.array
        Array of points that fall within the cone.

    Notes
    -----
    A point :math:`\textbf{X}=\langle x, y, z\rangle^T` falls within the cone
    where :math:`\textbf{U}=\langle x_0, y_0, z_0\rangle^T` is the apex of the
    cone and :math:`\theta` is the aperature angle if it satisfies the
    following inequality:

    .. math::

        0\leq\left(\textbf{X}-\textbf{U}\right)^T\left(\hat{\textbf{D}}
        \hat{\textbf{D}}^T-\cos^2\left(\theta\right)\textbf{I}_{3\times3}
        \right)\left(\textbf{X}-\textbf{U}\right)
    """

    U_norm = np.repeat(np.linalg.norm(U, axis=1).reshape((-1, 1)), 3, axis=1)
    D = np.divide(U, U_norm)
    Y = X - U

    d1, d2, d3 = D[:, 0], D[:, 1], D[:, 2]
    y1, y2, y3 = Y[:, 0], Y[:, 1], Y[:, 2]

    ct2 = np.cos(np.radians(theta)) ** 2

    r1 = (d1**2 - ct2) * y1 + d1 * d2 * y2 + d1 * d3 * y3
    r2 = d2 * d1 * y1 + (d2**2 - ct2) * y2 + d2 * d3 * y3
    r3 = d3 * d1 * y1 + d3 * d2 * y2 + (d3**2 - ct2) * y3

    point = y1 * r1 + y2 * r2 + y3 * r3

    ind = np.where(point >= 0)[0]

    return ind


def plane_constraint(U, X):
    """Return an array of indices containing points that fall above the plane.

    Parameters
    ----------
    U, X: np.array
        Arrays of shape (n, 3) containing cartesian vectors defining a plane.

    Returns
    -------
    np.array
        Array of points that fall above the plane.

    Notes
    -----
    A point :math:`\textbf{X}=\langle x, y, z\rangle^T` falls falls above the
    plane defined with :math:`\textbf{U}=\langle x_0, y_0, z_0\rangle^T` if it
    satisfies the following inequality:

    .. math:: 0\leq\left(\textbf{X}-\textbf{U}\right)^T\cdot\textbf{U}
    """

    Y = X - U

    u1, u2, u3 = U[:, 0], U[:, 1], U[:, 2]
    y1, y2, y3 = Y[:, 0], Y[:, 1], Y[:, 2]

    point = y1 * u1 + y2 * u2 + y3 * u3

    ind = np.where(point >= 0)[0]

    return ind


def aperature_theta(ang, ang_type, n=None, U=None, X=None):
    """Calculate aperature angle from constraint angles.

    Parameters
    ----------
    ang: float
        The constraint angle in degrees.
    angType: {0, 1}
        Defines the angle as altitude if `angType=0` and off-nadir if
        `angType=1`.
    n: int, optional
        Length of returned theta array.
    U: np.array, optional
        Array of shape (n, 3) containing the cartesian coordinates of the apex.
    X: np.array, optional
        Array of shape (n, 3) containing the cartesian coordinate points to
        check.

    Returns
    -------
    np.array
        Array of shape (n,) containing theta angles in degrees.

    Notes
    -----
    The `n` parameters is neccessary when `angType=0`. The `U` and `X`
    parameters are neccessary when `angtype=1`.

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

    if ang_type:
        ang = np.radians(ang)

        num = np.linalg.norm(X - U, axis=1)
        denom = np.linalg.norm(U, axis=1)

        theta = ang + np.arcsin(np.sin(ang) * num / denom)
        theta = np.degrees(theta)

    else:
        theta = np.full((n,), 90 - ang)

    return theta


def encounter_ind(satECEF, gndECEF, ang, angType):
    """Return encounter indices.

    This method returns the encounter indices corresponding only to a valid
    satellite position. It does not consider the sun constraint angle or
    lightinh conditions.

    Parameters
    ----------
    satECEF: np.array
        Array of shape (n, 3) containing satellite ECEF positions.
    gndECEF: np.array
        Array of shape (n, 3) containing ground location ECEF positions.
    ang: float
        Constraint angle in degrees.
    angType: {0, 1}
        Defines the angle as altitude if `angType=0` and off-nadir if
        `angType=1`.

    Returns
    -------
    np.array
        Array of indices defining valid encounter positions.

    Notes
    -----
    This method finds all indices that fall within a single sided cone
    extending above the Earth's surface with the apex at a ground location of
    interest and its aperature angle controlled by the constraint angle.
    """

    n = satECEF.shape[0]
    U = gndECEF
    X = satECEF
    theta = aperature_theta(ang, angType, n, U, X)

    ind_1 = cone_constraint(theta, U, X)
    ind_2 = plane_constraint(U, X)

    ind_final = np.intersect1d(ind_1, ind_2)

    return ind_final
