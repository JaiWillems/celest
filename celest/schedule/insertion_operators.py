

from celest.units.quantity import Quantity
from celest import units as u
import math
import numpy as np


def greedy_insertion(request_handler, number_to_insert):
    request_handler.sort_by_decreasing_priority()
    insert_first_n_requests(request_handler, number_to_insert)


def insert_first_n_requests(request_handler, number_to_insert):
    for request_index in range(request_handler.number_of_requests):

        if not number_to_insert:
            break

        if request_handler.is_request_scheduled(request_index):
            continue

        deadline = request_handler.deadline(request_index)
        duration = request_handler.duration(request_index)
        minimum_image_quality = request_handler.image_quality(request_index)
        look_angle = request_handler.look_angle(request_index)
        request_vtw_list = request_handler.vtw_list(request_index)

        duration_in_days = Quantity(duration.to(u.s) / 86400, u.jd2000)

        for vtw_index, vtw in enumerate(request_vtw_list):

            if look_angle is not None:
                start = _look_angle_time(
                    vtw.attitude.pitch,
                    look_angle,
                    vtw.attitude.time,
                    vtw.rise_time,
                    vtw.set_time
                )
                if start is None:
                    continue
            else:
                start = (vtw.rise_time + vtw.set_time + duration_in_days) / 2

            if (start < vtw.rise_time) or (start + duration_in_days > vtw.set_time):
                continue
            if start + duration_in_days > deadline:
                continue
            if not image_quality_is_met(vtw, start, minimum_image_quality):
                continue
            if does_conflict_exists(request_handler, start, duration_in_days):
                continue

            request_handler.schedule_request(request_index, vtw_index, start, duration)
            number_to_insert -= 1
            break


def _look_angle_time(look_angle, desired_look_angle, julian, start, end):
    julian = julian.to(u.jd2000)
    look_angle = look_angle.to(u.deg) - desired_look_angle.to(u.deg)

    left_value_index = np.where(julian == start.to(u.jd2000))[0][0]
    right_value_index = np.where(julian == end.to(u.jd2000))[0][0]

    left_look_angle = look_angle[left_value_index]
    right_look_angle = look_angle[right_value_index]

    if left_look_angle * right_look_angle > 0:
        return None

    while right_value_index - left_value_index > 1:
        mid_value_index = (left_value_index + right_value_index) // 2
        mid_look_angle = look_angle[mid_value_index]
        if left_look_angle * mid_look_angle > 0:
            left_value_index = mid_value_index
            left_look_angle = mid_look_angle
        else:
            right_value_index = mid_value_index

    return Quantity(julian[left_value_index], u.jd2000)


def image_quality_is_met(vtw, start, minimum_image_quality):
    return image_quality(vtw, start) >= minimum_image_quality


def image_quality(vtw, start):
    nadir_time = ((vtw.rise_time + vtw.set_time) / 2).to(u.jd2000)
    start_time = start.to(u.jd2000)
    rise_time = vtw.rise_time.to(u.jd2000)
    return math.floor(10 - 9 * abs(start_time - nadir_time) / (nadir_time - rise_time))


def does_conflict_exists(request_handler, start, duration_in_days):
    for request in request_handler:
        if not request[0]:
            continue
        scheduled_start = request[2]
        scheduled_duration = Quantity(request[3].to(u.s) / 86400, u.jd2000)
        scheduled_end = scheduled_start + scheduled_duration

        if start <= scheduled_start < start + duration_in_days:
            return True
        if start < scheduled_end <= start + duration_in_days:
            return True
        if scheduled_start <= start and start + duration_in_days <= scheduled_end:
            return True

    return False


def minimum_opportunity_insertion(request_handler, number_to_insert):
    request_handler.sort_by_decreasing_opportunity()
    insert_first_n_requests(request_handler, number_to_insert)


def minimum_conflict_insertion(request_handler, number_to_insert):
    request_handler.sort_vtws_by_increasing_conflict_degree()
    insert_first_n_requests(request_handler, number_to_insert)


_INSERTION_FUNCTIONS = [
    greedy_insertion,
    minimum_opportunity_insertion,
    minimum_conflict_insertion
]
