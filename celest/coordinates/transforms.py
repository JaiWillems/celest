

from celest.coordinates.astronomical_quantities import earth_rotation_angle
from celest.coordinates.azel import AzEl
from celest.coordinates.gcrs import GCRS
from celest.coordinates.itrs import ITRS
from celest.coordinates.wgs84 import WGS84
from celest.coordinates.nutation_precession_matrices import (
    bias_matrix,
    precession_matrix,
    nutation_matrix
)
from celest.coordinates.ground_location import GroundLocation
from celest.constants import WGS84_MINOR_AXIS_KM, WGS84_MAJOR_AXIS_KM
from celest import units as u
from copy import deepcopy
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

    gcrs = deepcopy(gcrs)
    julian = gcrs.time

    gcrs_x = gcrs.x.to(u.km).data.reshape((-1, 1))
    gcrs_y = gcrs.y.to(u.km).data.reshape((-1, 1))
    gcrs_z = gcrs.z.to(u.km).data.reshape((-1, 1))
    gcrs = np.concatenate((gcrs_x, gcrs_y, gcrs_z), axis=1)

    era = earth_rotation_angle(julian)
    bias_mtrx = bias_matrix()
    precession_mtrx = precession_matrix(julian)
    nutation_mtrx = nutation_matrix(julian)

    itrs = np.einsum('ij, kj -> ki', bias_mtrx, gcrs)
    itrs = np.einsum('ijk, ik -> ij', precession_mtrx, itrs)
    itrs = np.einsum('ijk, ik -> ij', nutation_mtrx, itrs)

    # TODO: Refactor negation when Quantity operations are implemented.
    negative_era = deepcopy(era)
    negative_era._data = -negative_era.data

    itrs = _rotate_by_earth_rotation_angle(itrs, negative_era)

    return ITRS(julian.data, itrs[:, 0], itrs[:, 1], itrs[:, 2], u.km)


def _rotate_by_earth_rotation_angle(data, degree_era):
    """Rotate vectors about the z-azis by the Earth rotation angle.

    Parameters
    ----------
    data : np.ndarray
        2-D array containing rows of 3d coordinate data.
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

    itrs = deepcopy(itrs)
    julian = itrs.time

    itrs_x = itrs.x.to(u.km).data.reshape((-1, 1))
    itrs_y = itrs.y.to(u.km).data.reshape((-1, 1))
    itrs_z = itrs.z.to(u.km).data.reshape((-1, 1))
    itrs = np.concatenate((itrs_x, itrs_y, itrs_z), axis=1)

    era = earth_rotation_angle(julian)
    bias_mtrx = bias_matrix()
    precession_mtrx = precession_matrix(julian)
    nutation_mtrx = nutation_matrix(julian)

    gcrs = _rotate_by_earth_rotation_angle(itrs, era)

    gcrs = np.einsum('ijk, ik -> ij', np.linalg.inv(nutation_mtrx), gcrs)
    gcrs = np.einsum('ijk, ik -> ij', np.linalg.inv(precession_mtrx), gcrs)
    gcrs = np.einsum('ij, kj -> ki', np.linalg.inv(bias_mtrx), gcrs)

    # TODO: Return the same unit as passed in.
    return GCRS(julian.data, gcrs[:, 0], gcrs[:, 1], gcrs[:, 2], u.km)


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

    julian = itrs.time.data

    itrs_x = itrs.x.to(u.km).data
    itrs_y = itrs.y.to(u.km).data
    itrs_z = itrs.z.to(u.km).data

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

    a, b = WGS84_MAJOR_AXIS_KM, WGS84_MINOR_AXIS_KM

    tol = 10e-10

    itrs_x = itrs.x.to(u.km).data
    itrs_y = itrs.y.to(u.km).data
    itrs_z = itrs.z.to(u.km).data

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

    julian = wgs84.time.data

    latitude = wgs84.latitude.to(u.rad).data
    longitude = wgs84.longitude.to(u.rad).data
    height = wgs84.height.to(u.km).data

    e = np.sqrt(1 - WGS84_MINOR_AXIS_KM ** 2 / WGS84_MAJOR_AXIS_KM ** 2)
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

    azimuth = _azimuth(itrs, location)
    elevation = _elevation(itrs, location)

    return AzEl(itrs.time.data, azimuth, elevation, u.deg, location)


def _azimuth(itrs: ITRS, location: GroundLocation) -> np.ndarray:
    """Return the azimuth angle of the satellite.

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

    sat_itrs = np.array([itrs.x.to(u.km).data, itrs.y.to(u.km).data,
                         itrs.z.to(u.km).data]).T
    ground_itrs = np.array([location.itrs_x.to(u.km).data,
                            location.itrs_y.to(u.km).data,
                            location.itrs_z.to(u.km).data])

    latitude = location.latitude.to(u.rad).data
    radius = location.radius.to(u.km).data

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
    """Return the degree angle bewteen two vectors.

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

    sat_itrs = np.array([itrs.x.to(u.km).data, itrs.y.to(u.km).data,
                         itrs.z.to(u.km).data]).T
    ground_itrs = np.array([location.itrs_x.to(u.km).data,
                            location.itrs_y.to(u.km).data,
                            location.itrs_z.to(u.km).data])

    return 90 - _get_ang(sat_itrs - ground_itrs, ground_itrs)
