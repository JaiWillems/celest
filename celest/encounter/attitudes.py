"""Generate pointing profiles for a satellite.

This module contains the functionality to determine satellite orientations for
ground tracking.
"""


from celest.satellite.coordinate import Coordinate
from scipy.spatial.transform import Rotation, Slerp
from typing import Any
import numpy as np



def skew(matrix: np.ndarray) -> np.ndarray:
    """Return the skew of a 3x3 input matrix.

    Parameters
    ----------
    matrix : np.ndarray
        The 3x3 matrix to skew.
    
    Returns
    -------
    np.ndarray
        The skewed version of `matrix`.
    
    Notes
    -----
    The skew of a column vector given by
    :math:`\vec{V}=\langle a_x, a_y, a_z\rangle^T` can be found as follows:

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


def sat_rotation(sat_z: np.ndarray) -> Rotation:
    """Return satellite rotations for encounter alignment.

    This function defines the satellite orientation as a rotation to be applied
    to the +Z gcrs axis to align with the Satellites +Z axis.

    Parameters
    ----------
    z_sat: np.ndarray
        Array of shape (n, 3) containing the Satellite's +Z vectors in the gcrs
        frame.
    
    Returns
    -------
    Rotation
        The rotations defined to transform the gcrs z-axis to the Satellite's
        z-axis.
    
    Notes
    -----
    Let the gcrs z-axis be :math:`\vec{a}` and the satellite's z-axis be
    :math:`\vec{b}`, then the rotation matrix, :math:`R`, is defined as follows

    .. math:: R = I + [\vec{v}]_{\times} + [\vec{v}]_{\times}^2\frac{1 - c}{s^2}

    where :math:`I` is the :math:`3x3` identity matrix,
    :math:`[\vec{v}]_{\times}` is the skew of
    :math:`\vec{v}=\vec{a}\times\vec{b},
    :math:`s = \left|\left|\vec{v}\right|\right|`, and
    :math:`c=\vec{a}\cdot\vec{b}`.[1]

    References
    ----------
    .. [1] “Calculate rotation matrix to align vector A to vector B in 3d?,”
       Stackexchange.com. [Online]. Available:
       https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d.
       [Accessed: 29-Jul-2021].
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


def attitudes(satellite: Any, location: Any, encInd: np.ndarray, maneuverTime: int) -> Rotation:
    """Generate satellite rotations for ground tracking.

    This function takes a ground location along and indices associated with
    delested windows at which the spacecraft makes imaging passes over the
    location. It then generates the satellite rotations to track the ground
    location.

    Parameters
    ----------
    location : GroundPosition
        Ground location and encounter information.
    encInd: np.ndarray
        Array of indices that correspond to times and positions where the
        spacecraft is in an imaging encounter window with the given ground
        location.
    maneuverTime: float
        Number of array indices to pad on either side of
        an encounter window to use for maneuvering time.

    Notes
    -----
    The strategy for pointing profile generation is as follows: (1) The default
    orientation is to have the spacecraft camera pointing towards its zenith.
    (2) When the spacecraft is imaging, orient the satellite such that the
    camera is facing the target. (3) Generate an initial set of pointing
    profiles assuming the above. (4) Interpolate the rotations between the
    normal and target-acquired states to smooth out transitions.
    """

    ground_GEO = [location.lat, location.lon, location.radius]
    ground_GEO = np.repeat(np.array([ground_GEO]), satellite.position.length, axis=0)
    target_site = Coordinate(ground_GEO, "geo", satellite.time)

    # Stage 1: preliminary rotations.
    # Get difference vector between spacecraft and target site.
    SC_to_site = target_site.gcrs() - satellite.position.gcrs()

    # Point toward the zenith except when over the imaging site.
    pointing_directions = satellite.position.gcrs()
    pointing_directions[encInd, :] = SC_to_site[encInd, :]

    # Preliminary rotation set.
    # Temporarily represent as quaternion for interpolation.
    rotations = sat_rotation(pointing_directions).as_quat()

    # Set flight modes. By default, point normal.
    flight_ind = 0 * np.ones(SC_to_site.shape[0])

    # Point to target during encounters.
    flight_ind[encInd] = 2

    # Stage 2: interpolation.
    # Generate sets of encounters which are clustered together.
    # This takes individual encInd into clusters which we can use later to
    # figure out when to start interpolation.
    split_ind = np.where(np.diff(encInd) > 1)[0]+1
    print(split_ind)
    encounter_segments = np.split(encInd, split_ind)
    print(encounter_segments)

    for seg in encounter_segments:
        print(seg)

        start_step = seg[0] - maneuverTime
        end_step = seg[-1] + maneuverTime

        # Get starting and ending quaternions.
        start_rotation = rotations[start_step]
        end_rotation = rotations[end_step]

        slerp_1 = Slerp([start_step, seg[0]], Rotation.from_quat([start_rotation, rotations[seg[0]]]))
        interp_rotations_1 = slerp_1(np.arange(start_step, seg[0]))

        rotations[start_step:seg[0]] = interp_rotations_1.as_quat()
        flight_ind[start_step:seg[0]] = 1

        slerp_2 = Slerp([seg[-1], end_step], Rotation.from_quat([rotations[seg[-1]], end_rotation]))
        interp_rotations_2 = slerp_2(np.arange(seg[-1], end_step))

        rotations[seg[-1]:end_step] = interp_rotations_2.as_quat()
        flight_ind[seg[-1]:end_step] = 3

    return Rotation.from_quat(rotations)
