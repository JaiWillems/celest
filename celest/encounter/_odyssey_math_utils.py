"""Math utils for the `Satellite.generate_pointing_profiles` functionality."""


from scipy.spatial.transform import Rotation
import numpy as np


def skew(matrix: np.array) -> np.array:
    """Return the skew of a 3x3 input matrix.

    Parameters
    ----------
    matrix : np.array
        3x3 matrix to skew.
    
    Returns
    -------
    np.array
        The skewed version of `matrix`.
    
    Notes
    -----
    The skew of a column vector given by :math:`\vec{V}=\langle a_x, a_y, a_z\rangle^T`
    can be found by:

    .. math::
    
        \vec{v}_{skew} = \begin{bmatrix}
                            0 & -a_z & a_y \\
                            a_z & 0 & -a_x \\
                            -a_y & a_x & 0 \\
                         \end{bmatrix}
    """

    n = matrix.shape[0]
    zeros = np.zeros(n)
    a_x, a_y, a_z = matrix[:, 0], matrix[:, 1], matrix[:, 2]

    matrix_skewed = np.array([[zeros, -a_z, a_y],
                              [a_z, zeros, -a_x],
                              [-a_y, a_x, zeros]])

    return np.transpose(matrix_skewed, axes=(2, 0, 1))


def sat_rotation(sat_z: np.array) -> Rotation:
    """Return satellite rotations for encounter alignment.

    This function defines the satellite orientation as a rotation to be applied
    to the +Z ECI axis to align with the Satellites +Z axis.

    Parameters
    ----------
    z_sat: np.array
        Array of shape (n, 3) containing the Satellite's +Z vectors in the ECI
        frame.
    
    Returns
    -------
    Rotation
        The rotations defined to transform the ECI z-axis to the Satellite's
        z-axis.
    
    Notes
    -----
    Let the ECI z-axis be :math:`\vec{a}` and the satellite's z-axis be
    :math:`\vec{b}`, then the rotation matrix, :math:`R`, is defined as follows

    .. math:: R = I + [\vec{v}]_{\times} + [\vec{v}]_{\times}^2\frac{1 - c}{s^2}

    where :math:`I` is the :math:`3x3` identity matrix, :math:`[\vec{v}]_{\times}`
    is the skew of :math:`\vec{v}=\vec{a}\times\vec{b},
    :math:`s = \left|\left|\vec{v}\right|\right|`, and
    :math:`c=\vec{a}\cdot\vec{b}`.[1]

    References
    ----------
    .. [1] “Calculate rotation matrix to align vector A to vector B in 3d?,”
       Stackexchange.com. [Online]. Available: https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d. [Accessed: 29-Jul-2021].
    """

    # Reference upward direction.
    z_ECI = np.array([0, 0, 1])

    # Normalize the vectors.
    normalized_sat_z = sat_z / np.linalg.norm(sat_z, axis=1)[:, np.newaxis]

    # Cross of difference between the unit vectors.
    v = np.cross(z_ECI[np.newaxis, :], normalized_sat_z, axis=1)

    sin = np.linalg.norm(v, axis=1)  # Sin of angle.
    cos = np.sum(z_ECI[np.newaxis, :] * normalized_sat_z, axis=1)  # Cos of angle.

    identity = np.eye(3, 3)[np.newaxis, :, :]
    v_skew = skew(v)
    v_skew_squared = np.matmul(v_skew, v_skew)
    ratio = ((1 - cos) / sin**2)[:, np.newaxis, np.newaxis]

    matrix = identity + v_skew + v_skew_squared * ratio

    return Rotation.from_matrix(matrix)
