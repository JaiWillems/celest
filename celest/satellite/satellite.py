"""Satellite orbital representations and coordinate conversions.

This module contains the `Satellite` class to localize satellite-related
information and functionality.
"""
#NAT was here hehe :)

from celest.core.decorators import set_module
from celest.satellite import Coordinate, sat_rotation
from scipy.spatial.transform import Rotation, Slerp
from typing import Any, List, Literal
import pandas as pd
import numpy as np


@set_module('celest.satellite')
class Satellite(object):
    """Localize satellite information and functionality.

    The `Satellite` class represents a satellite, be it artificial or natural,
    and allows for the position to be represented with time through multiple
    representations.

    Parameters
    ----------
    position : Coordinate
        Satellite time-evolving positions.

    Attributes
    ----------
    time : Time
        Times associated with the satellite positions.
    position : Coordinate
        Position of the satellite.

    Methods
    -------
    generate_pointing_profiles(groundPos, encInd, maneuverTime)
        Generate satellite rotations for ground tracking.
    save_data(fileName, delimiter, posTypes)
        Save satellite time and position data.
    """

    def __init__(self, position: Any) -> None:
        """Initialize attributes."""

        self.times = position.times
        self.position = position

    def generate_pointing_profiles(self, groundPos: Any, encInd: np.array,
                                   maneuverTime: int) -> Rotation:
        """Generate satellite rotations for ground tracking.

        This function is intended to take in a single ground location along
        with the windows at which the spacecraft makes imaging passes over the
        location. This method uses the Odyssey Pointing Profile
        determination system created by Mingde Yin.

        NOTE: This should be generalized in the future to many sites.

        Parameters
        ----------
        groundPos : GroundPosition
            Ground location and encounter information.
        encInd: np.array
            Array of arrays of indices that correspond to times and positions
            where the spacecraft is in an imaging encounter window with the
            given ground location.
        maneuverTime: float
            Number of array indices to pad on either side of
            an encounter window to use for maneuvering time.
            TODO: Turn into a standard time unit since setting a fixed number
            of array indices is limited.

        Notes
        -----
        The strategy for pointing profile generation is as follows:
        1. The default orientation is to have the spacecraft camera pointing
        towards its zenith.
        2. When the spacecraft is imaging, orient the satellite such that the
        camera is facing the target.
        3. Generate an initial set of pointing profiles assuming the above.
        4. Interpolate the rotations between the normal and target-acquired
           states to smooth out transitions.
        """

        ground_GEO = [groundPos.coor[0], groundPos.coor[1], groundPos.radius]
        ground_GEO = np.repeat(np.array([ground_GEO]), self.position.length, axis=0)
        target_site = Coordinate(ground_GEO, "GEO", self.times)

        # Stage 1: preliminary rotations.
        # Get difference vector between spacecraft and target site.
        SC_to_site = target_site.ECI() - self.position.ECI()

        # Point toward the zenith except when over the imaging site.
        pointing_directions = self.position.ECI()
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
        encounter_segments = np.split(encInd, split_ind)

        for encInd in encounter_segments:

            start_step = encInd[0] - maneuverTime
            end_step = encInd[-1] + maneuverTime

            # Get starting and ending quaternions.
            start_rotation = rotations[start_step]
            end_rotation = rotations[end_step]

            slerp_1 = Slerp([start_step, encInd[0]], Rotation.from_quat([start_rotation, rotations[encInd[0]]]))
            interp_rotations_1 = slerp_1(np.arange(start_step, encInd[0]))

            rotations[start_step:encInd[0]] = interp_rotations_1.as_quat()
            flight_ind[start_step:encInd[0]] = 1

            slerp_2 = Slerp([encInd[-1], end_step], Rotation.from_quat([rotations[encInd[-1]], end_rotation]))
            interp_rotations_2 = slerp_2(np.arange(encInd[-1], end_step))

            rotations[encInd[-1]:end_step] = interp_rotations_2.as_quat()
            flight_ind[encInd[-1]:end_step] = 3

        return Rotation.from_quat(rotations)

    def save_data(self, fileName: str, delimiter: Literal[",", "\\t"], posTypes: List) -> None:
        """Save satellite time and position data.

        Parameters
        ----------
        fileName : str
            File name of the output file as either a .txt or .csv file.
        delimiter : str
            String of length 1 representing the field delimiter for the output
            file.
        posTypes : List
            List containing the types of position data to store. Possible
            list values include "GEO", "ECI", "ECEF", and "Altitude".

        Notes
        -----
        It is recommended to use a tab delimiter for .txt files and comma
        delimiters for .csv files. The method will return an error if the
        fileName already exists in the current working directory.

        Examples
        --------
        >>> posTypes = ["ECI", "ECEF"]
        >>> finch.save_data(fileName="data.csv", delimiter=",", posTypes=posTypes)
        """

        data = {}
        data["Time (julian)"] = pd.Series(self.times.julian())
        if "GEO" in posTypes:
            GEO_pos = self.position.GEO()
            data["GEO.lat"] = pd.Series(GEO_pos[:, 0])
            data["GEO.lon"] = pd.Series(GEO_pos[:, 1])
            data["GEO.height"] = pd.Series(GEO_pos[:, 2])
        if "ECI" in posTypes:
            ECI_pos = self.position.ECI()
            data["ECI.X"] = pd.Series(ECI_pos[:, 0])
            data["ECI.y"] = pd.Series(ECI_pos[:, 1])
            data["ECI.z"] = pd.Series(ECI_pos[:, 2])
        if "ECEF" in posTypes:
            ECEF_pos = self.position.ECEF()
            data["ECEF.X"] = pd.Series(ECEF_pos[:, 0])
            data["ECEF.y"] = pd.Series(ECEF_pos[:, 1])
            data["ECEF.z"] = pd.Series(ECEF_pos[:, 2])
        if "Altitude" in posTypes:
            data["Altitude"] = pd.Series(self.position.altitude())

        df = pd.DataFrame(data)
        df.to_csv(fileName, sep=delimiter)
