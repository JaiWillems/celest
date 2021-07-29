
import numpy as np
from scipy.spatial.transform import Rotation


def skew(x: np.ndarray) -> np.ndarray:
    '''
    Returns the skew of a 3x3 input matrix
    '''
    return np.transpose(np.array(
        [[np.zeros(x.shape[0]), -x[:, 2], x[:, 1]],
         [x[:, 2], np.zeros(x.shape[0]), -x[:, 0]],
         [-x[:, 1], x[:, 0], np.zeros(x.shape[0])]]), axes=(2, 0, 1))


def sat_rotation(sat_z: np.ndarray) -> Rotation:
    '''
    Calculate Rotation object based on a series of target orientation vectors wrt +z in ECI
    Based on https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
    '''
    ref_vector = np.array((0, 0, 1))
    # reference upward direction

    normalized_sat_z = sat_z / np.linalg.norm(sat_z, axis=1)[:, np.newaxis]
    # normalize the vectors

    v = np.cross(ref_vector[np.newaxis, :], normalized_sat_z, axis=1)
    # cross of difference between the unit vectors

    sine = np.linalg.norm(v, axis=1)  # sin of angle
    cosine = np.sum(ref_vector[np.newaxis, :] * normalized_sat_z, axis=1)
    # cos of angle

    matrix = np.eye(3, 3)[np.newaxis, :, :] + skew(v) + \
        (np.matmul(skew(v), skew(v))) * \
        ((1-cosine)/sine**2)[:, np.newaxis, np.newaxis]

    return Rotation.from_matrix(matrix)
