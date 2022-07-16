

from celest.coordinates.ground_location import GroundLocation
from celest.coordinates.frames.attitude import Attitude
from celest.encounter.window_handling import VisibleTimeWindow
from celest.units.quantity import Quantity
from celest import units as u
import numpy as np


toronto = GroundLocation(43.65, -79.38, 0.076, u.deg, u.km)
north_bay = GroundLocation(46.31, -79.46, 0.193, u.deg, u.km)
sudbury = GroundLocation(46.49, -80.99, 0.348, u.deg, u.km)
ottawa = GroundLocation(45.42, -75.70, 0.070, u.deg, u.km)
kingston = GroundLocation(44.23, -76.49, 0.093, u.deg, u.km)
niagara_falls = GroundLocation(43.09, -79.08, 0.099, u.deg, u.km)
london = GroundLocation(42.98, -81.24, 0.251, u.deg, u.km)
mississauga = GroundLocation(43.59, -79.64, 0.156, u.deg, u.km)
timmins = GroundLocation(48.48, -81.33, 0.295, u.deg, u.km)
tobermory = GroundLocation(45.25, -81.66, 0.271, u.deg, u.km)


VTW_FILE_PATH = "tests/test_data/test_vtw_data/"
VTW_FILE_FINALE = "_vtw.csv"
ATTITUDE_FILE_PATH = "tests/test_data/test_attitude_data/"
ATTITUDE_FILE_FINALE = "_attitude.csv"
VTW_DATA = [
    ["toronto", toronto, 2460467, 30, 1, 1, None],
    ["north_bay", north_bay, 2460467, 30, 1, 1, None],
    ["sudbury", sudbury, 2460467, 30, 4, 1, None],
    ["ottawa", ottawa, 2460467, 30, 2, 1, None],
    ["kingston", kingston, 2460467, 30, 7, 1, None],
    ["niagara_falls", niagara_falls, 2460467, 30, 3, 1, None],
    ["london", london, 2460467, 30, 4, 1, None],
    ["mississauga", mississauga, 2460467, 30, 5, 1, None],
    ["timmins", timmins, 2460467, 30, 1, 1, None],
    ["tobermory", tobermory, 2460467, 30, 7, 1, None]
]


def initialize_request_list():
    return [construct_request_from_files(vtw_data_list) for vtw_data_list in VTW_DATA]


def construct_request_from_files(vtw_data_list):
    name, location, deadline, duration, priority, quality, look_ang = vtw_data_list
    vtws = get_vtw_list_from_file(name, location)

    return [
        location,
        Quantity(deadline, u.jd2000),
        Quantity(duration, u.s),
        priority,
        quality,
        Quantity(look_ang, u.deg) if look_ang is not None else None,
        vtws
    ]


def get_vtw_list_from_file(location_string, location):
    attitude = get_roll_pitch_yaw_from_file(location_string, location)
    vtws = load_file(VTW_FILE_PATH + location_string + VTW_FILE_FINALE)
    return [VisibleTimeWindow(vtw[0], vtw[1], attitude) for vtw in vtws]


def get_roll_pitch_yaw_from_file(location_string, location):
    attitude = load_file(ATTITUDE_FILE_PATH + location_string + ATTITUDE_FILE_FINALE)
    julian = np.linspace(0, len(attitude) - 1, len(attitude))
    roll = attitude[:, 0].reshape((-1,))
    pitch = attitude[:, 1].reshape((-1,))
    yaw = attitude[:, 2].reshape((-1,))
    return Attitude(julian, roll, pitch, yaw, u.deg, location)


def load_file(file_name):
    return np.loadtxt(file_name, delimiter=",")
