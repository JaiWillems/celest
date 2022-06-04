

import math


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

        duration_in_days = duration / 86400

        for vtw_index, vtw in enumerate(request_vtw_list):

            if look_angle is not None:
                start = _look_angle_time(vtw.pitch, look_angle, vtw.rise_time, vtw.set_time)
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
            if does_conflict_exists(request_handler, start, duration):
                continue

            request_handler.schedule_request(request_index, vtw_index, start, duration)
            number_to_insert -= 1
            break


def _look_angle_time(look_angle, desired_look_angle, start, end, tol=1e-6):

    look_angle = look_angle - desired_look_angle

    left_value = start
    right_value = end

    left_look_angle = look_angle(left_value)
    right_look_angle = look_angle(right_value)

    if left_look_angle * right_look_angle > 0:
        return None

    while (right_value - left_value > tol):
        mid_value = (left_value + right_value) / 2
        mid_look_angle = look_angle(mid_value)
        if left_look_angle * mid_look_angle > 0:
            left_value = mid_value
            left_look_angle = mid_look_angle
        else:
            right_value = mid_value
            right_look_angle = mid_look_angle

    return (left_value + right_value) / 2


def image_quality_is_met(vtw, start, minimum_image_quality):

    return image_quality(vtw, start) >= minimum_image_quality


def image_quality(vtw, start):

    nadir_time = (vtw.rise_time + vtw.set_time) / 2
    return math.floor(10 - 9 * abs(start - nadir_time) / (nadir_time - vtw.rise_time))


def does_conflict_exists(request_handler, start, duration_in_days):

    for request in request_handler:
        if not request[0]:
            continue
        scheduled_start = request[2]
        scheduled_end = scheduled_start + request[3] / 86400

        if start < scheduled_start and scheduled_start < start + duration_in_days:
            return True
        if start < scheduled_end and scheduled_end < start + duration_in_days:
            return True
        if scheduled_start < start and start + duration_in_days < scheduled_end:
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
