

from enum import IntEnum


class RequestIndices(IntEnum):

    is_scheduled = 0
    vtw_index = 1
    scheduled_start_time = 2
    scheduled_duration = 3
    location = 4
    deadline = 5
    duration = 6
    priority = 7
    quality = 8
    look_angle = 9
    vtw_list = 10


class RequestHandler:
    """RequestHandler()

    Container for handling satellite request scheduling.

    This class is the central data structure that stores and manipulates task
    fulfillment opportunities and scheduling. It is used in the scheduling
    algorithm to determine the most optimal satellite operations schedule.
    """

    def __init__(self):
        self.requests = []
        self.number_of_requests = 0
        self.number_of_scheduled_requests = 0

    def __len__(self):
        return self.number_of_requests

    def __getitem__(self, index):
        return self.requests[index]

    def __iter__(self):
        self.iteration_index = 0
        return self

    def __next__(self):
        current_index = self.iteration_index
        self.iteration_index += 1
        if self.iteration_index < self.number_of_requests:
            return self.requests[current_index]
        raise StopIteration

    def add_request(self, location, deadline, duration, priority, quality, look_ang, vtws):
        self.number_of_requests += 1
        self.requests.append([
            False,
            None,
            None,
            None,
            location,
            deadline,
            duration,
            priority,
            quality,
            look_ang,
            vtws
        ])

    def schedule_request(self, request_index, vtw_index, start, duration):
        self.requests[request_index][RequestIndices.is_scheduled] = True
        self.requests[request_index][RequestIndices.vtw_index] = vtw_index
        self.requests[request_index][RequestIndices.scheduled_start_time] = start
        self.requests[request_index][RequestIndices.scheduled_duration] = duration
        self.number_of_scheduled_requests += 1

    def unschedule_all_requests(self):
        for index in range(self.number_of_requests):
            self.unschedule_request(index)
        self.number_of_scheduled_requests = 0

    def unschedule_request(self, request_index):
        self.requests[request_index][RequestIndices.is_scheduled] = False
        self.requests[request_index][RequestIndices.vtw_index] = None
        self.requests[request_index][RequestIndices.scheduled_start_time] = None
        self.requests[request_index][RequestIndices.scheduled_duration] = None
        self.number_of_scheduled_requests -= 1

    def is_request_scheduled(self, request_index):
        return self.requests[request_index][RequestIndices.is_scheduled]

    def vtw_index(self, request_index):
        return self.requests[request_index][RequestIndices.vtw_index]

    def start_time(self, request_index):
        return self.requests[request_index][RequestIndices.scheduled_start_time]

    def scheduled_duration(self, request_index):
        return self.requests[request_index][RequestIndices.scheduled_duration]

    def deadline(self, request_index):
        return self.requests[request_index][RequestIndices.deadline]

    def duration(self, request_index):
        return self.requests[request_index][RequestIndices.duration]

    def priority(self, request_index):
        return self.requests[request_index][RequestIndices.priority]

    def image_quality(self, request_index):
        return self.requests[request_index][RequestIndices.quality]

    def look_angle(self, request_index):
        return self.requests[request_index][RequestIndices.look_angle]

    def vtw_list(self, request_index):
        return self.requests[request_index][RequestIndices.vtw_list]

    def sort_by_decreasing_priority(self):
        key = lambda x: x[RequestIndices.priority]
        self.requests = sorted(self.requests, key=key, reverse=True)

    def sort_by_decreasing_opportunity(self):
        key = lambda x: len(x[RequestIndices.vtw_list])
        self.requests = sorted(self.requests, key=key, reverse=True)

    def sort_by_increasing_opportunity(self):
        key = lambda x: len(x[RequestIndices.vtw_list])
        self.requests = sorted(self.requests, key=key, reverse=False)

    def sort_by_decreasing_conflict_degree(self):
        conflict_degree = []
        conflict_degree_index = []
        for request_index, request in enumerate(self.requests):
            if request[RequestIndices.is_scheduled]:
                vtw_index = request[RequestIndices.vtw_index]
                vtw = request[RequestIndices.vtw_list][vtw_index]
                conflict_degree.append(self._conflict_degree(vtw))
            else:
                conflict_degree.append(0)

            conflict_degree_index.append(request_index)

        indices = [index for _, index in
                   sorted(zip(conflict_degree, conflict_degree_index))]
        indices.reverse()

        self.requests = [self.requests[index] for index in indices]

    def _conflict_degree(self, vtw) -> float:
        overlapping_vtws = self._set_of_overlapping_vtws(vtw)
        total_overlap_time = 0
        for other_vtw in overlapping_vtws:
            total_overlap_time +=self._julian_day_overlap_between_vtws(vtw, other_vtw)

        return total_overlap_time / len(overlapping_vtws)

    def _set_of_overlapping_vtws(self, current_vtw) -> list:
        overlap_vtws = []
        for request in self.requests:
            for vtw in request[RequestIndices.vtw_list]:
                if self._are_vtws_overlaping(current_vtw, vtw):
                    overlap_vtws.append(vtw)

        return overlap_vtws

    def _are_vtws_overlaping(self, vtw1, vtw2) -> bool:
        return not (vtw1.set_time < vtw2.rise_time) | (vtw2.set_time < vtw1.rise_time)

    def _julian_day_overlap_between_vtws(self, vtw1, vtw2) -> float:
        s1, s2 = vtw1.rise_time, vtw2.rise_time
        e1, e2 = vtw1.set_time, vtw2.set_time

        if (e1 < s2) | (e2 < s1):
            return 0
        elif (s1 < s2) & (e2 < e1):
            return e2 - s2
        elif (s2 < s1) & (e1 < e2):
            return e1 - s1
        elif s2 < e1:
            return e1 - s2
        elif s1 < e2:
            return e2 - s1

    def sort_vtws_by_increasing_conflict_degree(self):
        for request_index in range(self.number_of_requests):
            conflict_degree = []
            vtw_indices = []
            vtw_list = self.requests[request_index][RequestIndices.vtw_list]
            for vtw_index, vtw in enumerate(vtw_list):
                conflict_degree.append(self._conflict_degree(vtw))
                vtw_indices.append(vtw_index)

            self.requests[request_index][RequestIndices.vtw_list] = \
                [vtw for _, vtw in sorted(zip(conflict_degree, vtw_list))]

    def sort_vtws_by_increasing_rise_time(self):
        for request in self.requests:
            request[RequestIndices.vtw_list] = sorted(
                request[RequestIndices.vtw_list],
                key=lambda x: x.rise_time.data
            )
