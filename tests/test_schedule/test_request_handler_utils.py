

from celest.encounter.groundposition import GroundPosition
from celest.encounter._window_handling import VTW
import numpy as np


toronto = GroundPosition(latitude=43.65, longitude=-79.38, height=0.076)
north_bay = GroundPosition(latitude=46.31, longitude=-79.46, height=0.193)
sudbury = GroundPosition(latitude=46.49, longitude=-80.99, height=0.348)
ottawa = GroundPosition(latitude=45.42, longitude=-75.70, height=0.070)
kingston = GroundPosition(latitude=44.23, longitude=-76.49, height=0.093)
niagara_falls = GroundPosition(latitude=43.09, longitude=-79.08, height=0.099)
london = GroundPosition(latitude=42.98, longitude=-81.24, height=0.251)
mississauga = GroundPosition(latitude=43.59, longitude=-79.64, height=0.156)
timmins = GroundPosition(latitude=48.48, longitude=-81.33, height=0.295)
tobermory = GroundPosition(latitude=45.25, longitude=-81.66, height=0.271)


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
    vtws = get_vtw_list_from_file(name)

    return [
        location,
        deadline,
        duration,
        priority,
        quality,
        look_ang,
        vtws
    ]


def get_vtw_list_from_file(location_string):

    roll, pitch, yaw = get_roll_pitch_yaw_from_file(location_string)
    vtws = load_file(VTW_FILE_PATH + location_string + VTW_FILE_FINALE)
    return [VTW(vtw[0], vtw[1], roll, pitch, yaw) for vtw in vtws]


def get_roll_pitch_yaw_from_file(location_string):

    attitude = load_file(ATTITUDE_FILE_PATH + location_string + ATTITUDE_FILE_FINALE)
    return attitude[:, 0], attitude[:, 1], attitude[:, 2]


def load_file(file_name):

    return np.loadtxt(file_name, delimiter=",")
