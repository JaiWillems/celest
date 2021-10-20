

from celest.core.decorators import set_module
from queue import Queue
from time import perf_counter
from typing import Literal
import numpy as np


ORIENTATION_TRANSITION_TIME = 5


def new_enc_start(pass_id, time_dl, time_img):
    if pass_id == 'data link':
        time_dl = perf_counter()
    else:
        time_img = perf_counter()

    return time_dl, time_img


@set_module('celest.encounter.schedule')
def generate(path: str, delimiter: Literal[',', '\t']):

    scheduleData = np.empty((0, 5), dtype="<U25")
    windowQueue = Queue()

    windows = np.loadtxt(path, dtype=str, delimiter=delimiter)
    for row in windows:
        windowQueue.put(row)

    timeDL, timeIMG = perf_counter(), perf_counter()

    item_one = windowQueue.get()
    while not windowQueue.empty():
        time_one = float(item_one[2])
        item_two = windowQueue.get()
        time_two = float(item_two[2])
        dt = time_two - time_one

        startTime_one = round(timeDL if item_one[0] == 'data link' else timeIMG, 4)
        startTime_two = round(timeDL if item_two[0] == 'data link' else timeIMG, 4)

        if time_one + dt + ORIENTATION_TRANSITION_TIME / 86400 < time_two:
            window = np.array(['TYPE: %s - COOR: (%s, %s) - START: %s - END: %s - ELAPSED SECONDS: %s' % (item_one[0], item_one[1], item_one[2], item_one[3], item_one[4], item_one[5])])
            scheduleData = np.append(scheduleData, window, axis=0)
            timeDL, timeIMG = new_enc_start(item_one[0], timeDL, timeIMG)
            item_one = item_two
        elif startTime_one > startTime_two:
            window = np.array(['TYPE: %s - COOR: (%s, %s) - START: %s - END: %s - ELAPSED SECONDS: %s' % (item_one[0], item_one[1], item_one[2], item_one[3], item_one[4], item_one[5])])
            scheduleData = np.append(scheduleData, window, axis=0)
            timeDL, timeIMG = new_enc_start(item_one[0], timeDL, timeIMG)
            item_one = windowQueue.get()
        else:
            item_one = item_two

    window = window = np.array(['TYPE: %s - COOR: (%s, %s) - START: %s - END: %s - ELAPSED SECONDS: %s' % (item_one[0], item_one[1], item_one[2], item_one[3], item_one[4], item_one[5])])
    scheduleData = np.append(scheduleData, window, axis=0)

    np.savetxt('Satellite Encounter Schedule.txt', scheduleData, fmt='%s')
