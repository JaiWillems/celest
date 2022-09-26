

from celest.coordinates.astronomical_quantities import earth_rotation_angle
from celest.coordinates.frames.azel import AzEl
from celest.coordinates.frames.gcrs import GCRS
from celest.coordinates.frames.itrs import ITRS
from celest.coordinates.frames.lvlh import LVLH
from celest.coordinates.frames.wgs84 import WGS84
from celest.coordinates.nutation_precession_matrices import (
    bias_matrix,
    precession_matrix,
    nutation_matrix
)
from celest.coordinates.ground_location import GroundLocation
from celest.constants import (
    WGS84_MINOR_AXIS_KM,
    WGS84_MAJOR_AXIS_KM,
    WGS83_FIRST_ECCENTRICITY
)
from celest import units as u
from copy import deepcopy
from typing import Tuple
import numpy as np


def _gcrs_to_itrs(gcrs: GCRS) -> ITRS:
    """Gcrs to itrs transformation.

    Parameters
    ----------
    gcrs : GCRS
        Gcrs coordinates.

    Returns
    -------
    ITRS
        Itrs coordinates.

    See Also
    --------
    _itrs_to_gcrs : Itrs to gcrs transformation.

    Notes
    -----
    The conversion between the gcrs and itrs frames accounts for Earth
    rotation, precession, and nutation. Polar motion effects are ignored
    due to their poor predictability.
    """

    if not isinstance(gcrs, GCRS):
        raise ValueError(f"Input data is in the {gcrs.__class__} frame and not"
                         " the GCRS frame.")

    unit = gcrs.x.unit

    gcrs = deepcopy(gcrs)
    julian = gcrs.time

    gcrs_x = gcrs.x.data.reshape((-1, 1))
    gcrs_y = gcrs.y.data.reshape((-1, 1))
    gcrs_z = gcrs.z.data.reshape((-1, 1))
    gcrs = np.concatenate((gcrs_x, gcrs_y, gcrs_z), axis=1)

    era = earth_rotation_angle(julian)
    bias_mtrx = bias_matrix()
    precession_mtrx = precession_matrix(julian)
    nutation_mtrx = nutation_matrix(julian)

    itrs = np.einsum('ij, kj -> ki', bias_mtrx, gcrs)
    itrs = np.einsum('ijk, ik -> ij', precession_mtrx, itrs)
    itrs = np.einsum('ijk, ik -> ij', nutation_mtrx, itrs)

    itrs = _rotate_by_earth_rotation_angle(itrs, -era)

    return ITRS(julian.data, itrs[:, 0], itrs[:, 1], itrs[:, 2], unit)


def _rotate_by_earth_rotation_angle(data, degree_era):
    """Rotate vectors about the z-axis by the Earth rotation angle.

    Parameters
    ----------
    data : np.ndarray
        2-D array containing columns of x, y, z coordinate data.
    degree_era : Quantity
        The Earth rotation angles.

    Returns
    -------
    np.ndarray
        Rotated data.
    """

    radian_era = degree_era.to(u.rad).data
    c, s = np.cos(radian_era), np.sin(radian_era)
    c1, c2, c3 = np.copy(data).T

    data[:, 0] = c * c1 - s * c2
    data[:, 1] = s * c1 + c * c2
    data[:, 2] = c3

    return data


def _itrs_to_gcrs(itrs: ITRS) -> GCRS:
    """Itrs to gcrs transformation.

    Parameters
    ----------
    itrs : ITRS
        Itrs coordinates.

    Returns
    -------
    GCRS
        Gcrs coordinates.

    See Also
    --------
    _gcrs_to_itrs : Gcrs to itrs transformation.

    Notes
    -----
    The conversion between the gcrs and itrs frames accounts for Earth
    rotation, precession, and nutation. Polar motion effects are ignored
    due to their poor predictability.
    """

    if not isinstance(itrs, ITRS):
        raise ValueError(f"Input data is in the {itrs.__class__} frame and not"
                         " the ITRS frame.")

    unit = itrs.x.unit

    itrs = deepcopy(itrs)
    julian = itrs.time

    itrs_x = itrs.x.data.reshape((-1, 1))
    itrs_y = itrs.y.data.reshape((-1, 1))
    itrs_z = itrs.z.data.reshape((-1, 1))
    itrs = np.concatenate((itrs_x, itrs_y, itrs_z), axis=1)

    era = earth_rotation_angle(julian)
    bias_mtrx = bias_matrix()
    precession_mtrx = precession_matrix(julian)
    nutation_mtrx = nutation_matrix(julian)

    gcrs = _rotate_by_earth_rotation_angle(itrs, era)

    gcrs = np.einsum('ijk, ik -> ij', np.linalg.inv(nutation_mtrx), gcrs)
    gcrs = np.einsum('ijk, ik -> ij', np.linalg.inv(precession_mtrx), gcrs)
    gcrs = np.einsum('ij, kj -> ki', np.linalg.inv(bias_mtrx), gcrs)

    return GCRS(julian.data, gcrs[:, 0], gcrs[:, 1], gcrs[:, 2], unit)


def _itrs_to_wgs84(itrs: ITRS) -> WGS84:
    """Itrs to wgs84 transformation.

    Parameters
    ----------
    itrs : ITRS
        Itrs coordinates.

    Returns
    -------
    WGS84
        WGS84 coordinates.

    See Also
    --------
    _wgs84_to_itrs : Wgs84 to itrs transformation.
    """

    if not isinstance(itrs, ITRS):
        raise ValueError(f"Input data is in the {itrs.__class__} frame and not"
                         " the ITRS frame.")

    julian = itrs.time.to(u.jd2000)

    itrs_x = itrs.x.to(u.km)
    itrs_y = itrs.y.to(u.km)
    itrs_z = itrs.z.to(u.km)

    latitude = np.arctan(itrs_z / np.sqrt(itrs_x ** 2 + itrs_y ** 2))
    longitude = np.arctan2(itrs_y, itrs_x)
    height = _altitude(itrs)

    return WGS84(julian, latitude, longitude, height, u.rad, u.km)


def _altitude(itrs: ITRS) -> np.ndarray:
    """Return the geodetic altitude in kilometers.

    This method uses the WGS84 reference ellipsoid to calculate the
    geodetic altitude above the Earth's surface.

    Parameters
    ----------
    itrs : ITRS
        Itrs coordinates.

    Returns
    -------
    np.ndarray
        1-D array containing the geodetic altitude in kilometres.

    Notes
    -----
    This method uses an ellipsoid based model of the Earth to calculate
    the ellipsoid height in an iterative manner described in "Coordinate
    Systems in Geodesy" by E. J. Krakiwsky and D.E. Wells. [KW98b]_

    References
    ----------
    .. [KW98b] E. J. Krakiwsky and D. E. Wells. Coordinate Systems in
       Geodesy. Jan. 1998, pp. 31–33.
    """

    if not isinstance(itrs, ITRS):
        raise ValueError(f"Input data is in the {itrs.__class__} frame and not"
                         " the ITRS frame.")

    a, b = WGS84_MAJOR_AXIS_KM, WGS84_MINOR_AXIS_KM

    tol = 10e-10

    itrs_x = itrs.x.to(u.km)
    itrs_y = itrs.y.to(u.km)
    itrs_z = itrs.z.to(u.km)

    e = np.sqrt(1 - b ** 2 / a ** 2)
    p = np.sqrt(itrs_x ** 2 + itrs_y ** 2)

    h = np.sqrt(itrs_x ** 2 + itrs_y ** 2 + itrs_z ** 2) - np.sqrt(a * b)
    phi = np.arctan((itrs_z / p) * (1 - (e ** 2 * a) / (a + h)) ** -1)

    def new_vals(a, b, e, p, phi, z):
        N = a / np.sqrt(np.cos(phi) ** 2 + b ** 2 / a ** 2 * np.sin(phi) ** 2)
        h = p / np.cos(phi) - N
        phi = np.arctan((z / p) * (1 - (e ** 2 * N) / (N + h)) ** -1)

        return h, phi

    hp, phip = new_vals(a, b, e, p, phi, itrs_z)

    while (np.mean(hp - h) > a * tol) or (np.mean(phip - phi) > tol):
        h, phi = hp, phip
        hp, phip = new_vals(a, b, e, p, phi, itrs_z)

    return h


def _wgs84_to_itrs(wgs84: WGS84) -> ITRS:
    """Wgs84 to itrs transformation.

    Parameters
    ----------
    wgs84 : WGS84
        Wgs84 coordinates.

    Returns
    -------
    ITRS
        Itrs coordinates.

    See Also
    --------
    _itrs_to_wgs84 : Itrs to wgs84 transformation.

    Notes
    -----
    An Earth ellipsoid model is used for the wgs84 to itrs conversion using the
    methods described in "Coordinate Systems in Geodesy" by E. J. Krakiwsky and
    D.E. Wells as presented by Christopher Lum. [KW98a]_ [Lum20]_

    References
    ----------
    .. [KW98a] E. J. Krakiwsky and D. E. Wells. Coordinate Systems in
       Geodesy. Jan. 1998.
    .. [Lum20] Christopher Lum. Geodetic Coordinates: Computing Latitude and
       Longitude. June 2020.url:https://www.youtube.com/watch?v=4BJ-GpYbZlU.
    """

    if not isinstance(wgs84, WGS84):
        raise ValueError(f"Input data is in the {wgs84.__class__} frame and not"
                         " the WGS84 frame.")

    julian = wgs84.time.to(u.jd2000)

    latitude = wgs84.latitude.to(u.rad)
    longitude = wgs84.longitude.to(u.rad)
    height = wgs84.height.to(u.km)

    e = WGS83_FIRST_ECCENTRICITY
    n = WGS84_MAJOR_AXIS_KM / np.sqrt(1 - e ** 2 * np.sin(latitude) ** 2)

    itrs_x = (n + height) * np.cos(latitude) * np.cos(longitude)
    itrs_y = (n + height) * np.cos(latitude) * np.sin(longitude)
    itrs_z = (n * (1 - e ** 2) + height) * np.sin(latitude)

    return ITRS(julian, itrs_x, itrs_y, itrs_z, u.km)


def _itrs_to_azel(itrs: ITRS, location: GroundLocation) -> AzEl:
    """Itrs to azel transformation.

    The azimuth angle ranges from 0 to 360 degrees, measured clockwise
    from North. The elevation angle ranges from 0 to 90 degrees measured
    above the local horizon.

    Parameters
    ----------
    itrs : ITRS
        Itrs coordinates.
    location : GroundLocation
        Eath bound origin of the horizontal system.

    Returns
    -------
    AzEl
        Azel coordinates.
    """

    if not isinstance(itrs, ITRS):
        raise ValueError(f"Input data is in the {itrs.__class__} frame and not"
                         " the ITRS frame.")

    azimuth = _azimuth(itrs, location)
    elevation = _elevation(itrs, location)

    return AzEl(itrs.time.to(u.jd2000), azimuth, elevation, u.deg, location)


def _azimuth(itrs: ITRS, location: GroundLocation) -> np.ndarray:
    """Return the azimuth angle of the satellite_old.

    Parameters
    ----------
    itrs : ITRS
        Itrs coordinates.
    location : GroundLocation
        Earth bound origin of the horizontal system.

    Returns
    -------
    np.ndarray
        1-D array containing azimuth angles in degrees.
    """

    sat_itrs = np.array([itrs.x.to(u.km), itrs.y.to(u.km), itrs.z.to(u.km)]).T
    ground_itrs = np.array([
        location.itrs_x.to(u.km),
        location.itrs_y.to(u.km),
        location.itrs_z.to(u.km)
    ])

    latitude = location.latitude.to(u.rad)
    radius = location.radius.to(u.km)

    ground_to_satellite = sat_itrs - ground_itrs

    z_axis_vector = [0, 0, radius / np.sin(latitude)]
    surface_tangent = z_axis_vector - ground_itrs

    sat_to_ground_on_normal = np.dot(ground_to_satellite, ground_itrs) / radius ** 2
    sat_to_ground_on_normal = np.array([sat_to_ground_on_normal *
                                        ground_itrs[i] for i in range(3)]).T
    ground_to_sat_on_plane = ground_to_satellite - sat_to_ground_on_normal

    direction = np.cross(surface_tangent, ground_itrs)
    negative_indices = np.sum(ground_to_sat_on_plane * direction, axis=1) < 0
    azimuth = _get_ang(surface_tangent, ground_to_sat_on_plane)
    azimuth = 360 * negative_indices - 1 * (2 * negative_indices - 1) * azimuth

    return azimuth


def _get_ang(u, v) -> np.ndarray:
    """Return the degree angle between two vectors.

    Parameters
    ----------
    u, v : np.ndarray
        1-D or 2-D arrays containing row vectors.

        If the dimensions of `u` and `v` do not match, the 1-D array will
        be broadcast to have the same number of rows as the 2-D array.

    Returns
    -------
    np.ndarray
        1-D array containing the degree angle between to vector arrays.
    """

    ua = None if u.ndim == 1 else 1
    va = None if v.ndim == 1 else 1

    numerator = np.sum(u * v, axis=(ua or va))
    denominator = np.linalg.norm(u, axis=ua) * np.linalg.norm(v, axis=va)
    ang = np.degrees(np.arccos(numerator / denominator))

    return ang


def _elevation(itrs: ITRS, location: GroundLocation) -> np.ndarray:
    """Return the elevation angle.

    Parameters
    ----------
    itrs : ITRS
        Itrs coordinates.
    location : GroundLocation
        Earth bound origin of the horizontal system.

    Returns
    -------
    np.ndarray
        1-D array containing elevation angles in degrees.
    """

    sat_itrs = np.array([itrs.x.to(u.km), itrs.y.to(u.km), itrs.z.to(u.km)]).T
    ground_itrs = np.array([
        location.itrs_x.to(u.km),
        location.itrs_y.to(u.km),
        location.itrs_z.to(u.km)
    ])

    return 90 - _get_ang(sat_itrs - ground_itrs, ground_itrs)


def _gcrs_to_lvlh(gcrs_position: GCRS, gcrs_velocity: GCRS) -> Tuple[LVLH, LVLH]:
    """Return the LVLH (Hill frame) coordinates.

    Parameters
    ----------
    gcrs_position, gcrs_velocity : GCRS

    Returns
    -------
    Tuple
        Tuple containing the position and velocity `LVLH` objects.

    Notes
    -----
    The LVLH frame definition was taken from NASA's technical memorandum
    on coordinate frames for the space shuttle program. [NASA1974]_

    References
    ----------
    .. [NASA1974] Coordinate Systems for the Space Shuttle Program, Lyndon
       B. Johnson Space Center, Houston, Texas 77058, Oct 1974, no. NASA
       TM X-58153.
    """

    if not isinstance(gcrs_position, GCRS):
        raise ValueError(f"Input position is in the {gcrs_position.__class__} "
                         "frame and not the GCRS frame.")
    if not isinstance(gcrs_velocity, GCRS):
        raise ValueError(f"Input velocity is in the {gcrs_velocity.__class__} "
                         "frame and not the GCRS frame.")

    position_unit = gcrs_position.x.unit
    velocity_unit = gcrs_velocity.x.unit

    gcrs_position = deepcopy(gcrs_position)
    julian = gcrs_position.time.to(u.jd2000)

    gcrs_x = gcrs_position.x.data.reshape((-1, 1))
    gcrs_y = gcrs_position.y.data.reshape((-1, 1))
    gcrs_z = gcrs_position.z.data.reshape((-1, 1))
    position = np.concatenate((gcrs_x, gcrs_y, gcrs_z), axis=1)

    gcrs_vx = gcrs_velocity.x.data.reshape((-1, 1))
    gcrs_vy = gcrs_velocity.y.data.reshape((-1, 1))
    gcrs_vz = gcrs_velocity.z.data.reshape((-1, 1))
    velocity = np.concatenate((gcrs_vx, gcrs_vy, gcrs_vz), axis=1)

    transformation_matrix = _gcrs_to_lvlh_matrix(position, velocity)
    transformed_position = np.einsum(
        'ijk, ik -> ij',
        transformation_matrix,
        position
    )
    transformed_velocity = np.einsum(
        'ijk, ik -> ij',
        transformation_matrix,
        velocity
    )

    lvlh_position = LVLH(
        julian,
        transformed_position[:, 0],
        transformed_position[:, 1],
        transformed_position[:, 2],
        position_unit
    )
    lvlh_velocity = LVLH(
        julian,
        transformed_velocity[:, 0],
        transformed_velocity[:, 1],
        transformed_velocity[:, 2],
        velocity_unit
    )

    return lvlh_position, lvlh_velocity


def _gcrs_to_lvlh_matrix(gcrs_position: np.ndarray, gcrs_velocity:
                         np.ndarray) -> np.ndarray:
    """Return the gcrs-to-lvlh rotation matrix.

    Parameters
    ----------
    gcrs_position, gcrs_velocity : GCRS

    Returns
    -------
    np.ndarray
        3-D rotation matrix.
    """

    norm_position = np.linalg.norm(gcrs_position, axis=1)
    position_cross_velocity = np.cross(gcrs_position, gcrs_velocity)
    norm_position_cross_velocity = np.linalg.norm(position_cross_velocity, axis=1)

    lvlh_z = - gcrs_position / norm_position.reshape((-1, 1)).repeat(3, axis=1)
    lvlh_y = - position_cross_velocity / norm_position_cross_velocity.reshape((-1, 1)).repeat(3, axis=1)
    lvlh_x = np.cross(lvlh_y, lvlh_z)

    transformation_matrix = np.zeros((gcrs_position.shape[0], 3, 3))
    transformation_matrix[:, 0, :] = lvlh_x
    transformation_matrix[:, 1, :] = lvlh_y
    transformation_matrix[:, 2, :] = lvlh_z

    return transformation_matrix
