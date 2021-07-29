"""Satellite orbital representations and coordinate conversions.

The `Satellite` object localizes satellite related information and
functionality, and is passed into the `Encounter` class for encounter
planning.
"""
#NAT was here hehe :)

from celest.core.decorators import set_module
from celest.satellite import Coordinate
from typing import List, Literal
import pandas as pd
from scipy.spatial.transform import Rotation, Slerp
import numpy as np
from ._odyssey_math_utils import sat_rotation
from ._odyssey_enums import FlightStates


@set_module('celest.satellite')
class Satellite(object):
    """Localize satellite information and functionality.

    The `Satellite` class represents a satellite, be it artificial or natural,
    and allows for the position to be represented with time through multiple
    representations.

    Parameters
    ----------
    coordinates : Coordinate
        `Coordinate` object containing the position and time evolution of the
        satellite.

    Attributes
    ----------
    time : Time
        Times associated with the satellite positions.
    position : Coordinate
        Position of the satellite.

    Methods
    -------
    solar_power(pointProfiles)
        Calculate the solar power generated from the satellite solar cells.
    solar_radiation_pressure(pointProfiles)
        Calculate the solar radiation pressure experienced by the satellite.
    save_data(fileName, delimiter, posTypes)
        Save the time and position data of the satellite.
    """
    
    def __init__(self, coordinates: Coordinate) -> None:
        """Initialize attributes."""

        self.time = coordinates.timeData
        self.position = coordinates
    
    def save_data(self, fileName: str, delimiter: Literal[",", "\\t"], posTypes: List) -> None:
        """Save satellite data to local directory.

        Parameters
        ----------
        fileName : str
            File name of the output file as wither a .txt or .csv file.
        delimiter : str
            String of length 1 representing the feild delimiter for the output
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
        data["Time (julian)"] = pd.Series(self.time)
        if "GEO" in posTypes:
            GEO_pos = self.position.GEO()
            data["GEO.lat"] = pd.Series(GEO_pos[:, 0])
            data["GEO.lon"] = pd.Series(GEO_pos[:, 1])
            data["GEO.radius"] = pd.Series(GEO_pos[:, 2])
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

    def generate_pointing_profiles(self, target_site: Coordinate,
                                   encounter_indices: np.ndarray,
                                   maneuver_gap: int) -> Rotation:
        """Generates rotations using the Odyssey Pointing profile determination system

        This function is intended to take in a single ground site, along with
        the windows at which the spacecraft makes **IMAGING** passes
        over the site.

        NOTE: This should be generalized in the future to many sites.

        Parameters
        ----------
        target_site: Coordinate
            coordinate representation of ground site location.

        encounter_indices: np.ndarray
            array indices for which spacecraft is in
            an imaging encounter window with the given target site

        maneuver_gap: float
            Number of array indices to pad on either side of
            an encounter window to use for maneuvering time.
            TODO: turn into a standard time unit, like seconds,
            since setting a fixed number of array indices is bad
        """

        '''
        Strategy
        1. Normally, the spacecraft points towards Zenith (up)
        2. When the spacecraft can image, make it point camera (down) facing
           the target
        3. Generate an initial set of pointing profiles assuming this
        4. Interpolate between normal and target pointing to smooth
           things out after the matter
        '''

        # Stage 1: Preliminary rotations

        SC_to_site: np.ndarray = target_site.ECI() - self.position.ECI()
        # Get difference vector between spacecraft and target site

        pointing_directions = self.position.ECI()
        pointing_directions[encounter_indices, :] =\
            SC_to_site[encounter_indices, :]
        # Set the spacecraft to be zenith pointing, EXCEPT when over the site.

        rotations: np.ndarray = sat_rotation(pointing_directions).as_quat()
        # Preliminary rotation set
        # Temporarily represent as quaternion for interpolation

        # Set flight modes
        flight_indices = FlightStates.NORMAL_POINTING *\
            np.ones(SC_to_site.shape[0])
        # By default, point normal

        flight_indices[encounter_indices] = FlightStates.ENCOUNTER
        # Point to target during encounters

        # Stage 2: Interpolation
        encounter_segments = np.split(encounter_indices,
                                      np.where(
                                          np.diff(encounter_indices) > 1)[0]+1)
        # Generate sets of encounters which are clustered together
        # This takes individual indices into clusters which we can use
        # later to figure out when to start interpolation.

        for encounter_indices in encounter_segments:
            starting_step = encounter_indices[0] - maneuver_gap
            ending_step = encounter_indices[-1] + maneuver_gap

            starting_rotation = rotations[starting_step]
            # Starting quat
            ending_rotation = rotations[ending_step]
            # Ending quat

            slerp_1 = Slerp([starting_step, encounter_indices[0]], Rotation.from_quat([starting_rotation, rotations[encounter_indices[0]]]))
            interpolated_rotations_1: Rotation = slerp_1(np.arange(starting_step, encounter_indices[0]))

            rotations[starting_step:encounter_indices[0]] = interpolated_rotations_1.as_quat()
            flight_indices[starting_step:encounter_indices[0]] = FlightStates.PRE_ENCOUNTER

            slerp_2 = Slerp([encounter_indices[-1], ending_step], Rotation.from_quat([rotations[encounter_indices[-1]], ending_rotation]))
            interpolated_rotations_2: Rotation = slerp_2(np.arange(encounter_indices[-1], ending_step))

            rotations[encounter_indices[-1]:ending_step] = interpolated_rotations_2.as_quat()
            flight_indices[encounter_indices[-1]:ending_step] = FlightStates.POST_ENCOUNTER

        return Rotation.from_quat(rotations)
