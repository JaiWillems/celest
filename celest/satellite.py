

from celest.coordinates.frames.gcrs import GCRS
from celest.coordinates.frames.itrs import ITRS
from celest.coordinates.ground_location import GroundLocation
from celest.coordinates.transforms import (
    _gcrs_to_itrs,
    _altitude,
    _get_ang,
    _itrs_to_gcrs,
    _gcrs_to_lvlh,
    _gcrs_to_lvlh_matrix
)
from celest.file_save import _save_data_as_txt
from celest.units.quantity import Quantity
from celest import units as u
from typing import Union, Tuple
import numpy as np


class Satellite:
    """Satellite(position, velocity=None)

    Representation of an Earth orbiting satellite.

    Parameters
    ----------
    position : {GCRS, ITRS}
        Position of the satellite in the GCRS or ITRS frame.
    velocity : {GCRS, ITRS}, optional
        Velocity of the satellite in the GCRS or ITRS frame.

        The velocity is only necessary for the attitude calculations.

    Methods
    -------
    attitude(location)
        Return the satellite's attitude to a ground location.
    look_angle(location)
        Return the satellite's look angle to a ground location.
    altitude()
        Return the satellite's geodetic altitude.
    distance(location)
        Return the satellite's distance to a ground location.
    save_text_file(file_name)
        Save the satellite's position and velocity to a text file.
    """

    # TODO: Make type that represents the different coordinate frames.
    def __init__(self, position: Union[GCRS, ITRS], velocity:
                 Union[GCRS, ITRS]=None) -> None:

        if isinstance(position, ITRS):
            self.position = position
        elif isinstance(position, GCRS):
            self.position = _gcrs_to_itrs(position)
        else:
            raise ValueError(
                "Input position must be in the GCRS or ITRS frame.")

        if velocity is not None:
            if isinstance(velocity, ITRS):
                self.velocity = velocity
            elif isinstance(velocity, GCRS):
                self.velocity = _gcrs_to_itrs(velocity)
            else:
                raise ValueError(
                    "Input velocity must be in the GCRS or ITRS frame.")
        else:
            self.velocity = None

    def attitude(self, location: GroundLocation) -> Tuple[Quantity, Quantity,
                                                          Quantity]:
        """Return satellite roll, pitch, and yaw angles.

        This method returns the attitude angles required to align the
        satellite's nadir with a ground location for ground imaging.

        Parameters
        ----------
        location : GroundLocation
            Location for satellite attitude pointing.

        Returns
        -------
        Tuple[Quantity, Quantity, Quantity]
            Tuple containing roll, pitch, and yaw angles.

        Notes
        -----
        The methods of attitude determination were taken from [adcs1]_.

        References
        ----------
        .. [adcs1] G. H. J. van Vuuren, “The design and simulation analysis of
           an attitude determination and control system for a small earth
           observation satellite_old,” Master of Engineering, Stellenbosch
           University, Stellenbosch, South Africa, Mar 2015.
        """

        if self.velocity is None:
            raise ValueError("Velocity is required for attitude determination.")

        gcrs_position = _itrs_to_gcrs(self.position)
        gcrs_velocity = _itrs_to_gcrs(self.velocity)

        gcrs_position_array = gcrs_position.to_numpy(u.km)
        gcrs_velocity_array = gcrs_velocity.to_numpy(u.m/u.s)

        lvlh_position, lvlh_velocity = _gcrs_to_lvlh(gcrs_position,
                                                     gcrs_velocity)
        lvlh_position_array = lvlh_position.to_numpy(u.km)

        ground_itrs = ITRS(
            self.position.time.data,
            np.full((len(lvlh_position.x.data),),
                    location.itrs_x.to(u.km).data),
            np.full((len(lvlh_position.x.data),),
                    location.itrs_y.to(u.km).data),
            np.full((len(lvlh_position.x.data),),
                    location.itrs_z.to(u.km).data),
            u.km
        )
        ground_gcrs = _itrs_to_gcrs(ground_itrs)
        ground_gcrs = np.array([
            ground_gcrs.x.to(u.km).data,
            ground_gcrs.y.to(u.km).data,
            ground_gcrs.z.to(u.km).data
        ]).T

        transformation_matrix = _gcrs_to_lvlh_matrix(
            gcrs_position_array,
            gcrs_velocity_array
        )
        ground_lvlh = np.einsum('jki, ji -> jk', transformation_matrix,
                                ground_gcrs)

        satellite_lvlh_norm = np.linalg.norm(lvlh_position_array, axis=1)
        s = - lvlh_position_array / satellite_lvlh_norm[:, None]

        ground_norm = np.linalg.norm(ground_lvlh - lvlh_position_array, axis=1)
        g = (ground_lvlh - lvlh_position_array) / ground_norm[:, None]

        v, c = np.cross(s, g, axis=1), np.sum(s * g, axis=1)

        a32 = v[:, 0] + v[:, 1] * v[:, 2] / (1 + c)
        a31 = - v[:, 1] + v[:, 0] * v[:, 2] / (1 + c)
        a33 = 1 - (v[:, 0] ** 2 + v[:, 1] ** 2) / (1 + c)
        a12 = - v[:, 2] + v[:, 0] * v[:, 1] / (1 + c)
        a22 = 1 - (v[:, 0] ** 2 + v[:, 2] ** 2) / (1 + c)

        roll = Quantity(-np.arcsin(a32), u.rad)
        pitch = Quantity(np.arctan2(a31, a33), u.rad)
        yaw = Quantity(np.arctan2(a12, a22), u.rad)

        return roll, pitch, yaw

    def look_angle(self, location: GroundLocation) -> Quantity:
        """Return look-angles to a ground location.

        The look angle is the angular distance of a ground location from
        the satellites nadir.

        Parameters
        ----------
        location : GroundLocation
            Ground location of interest.

        Returns
        -------
        Quantity
            Satellite's look angle to `location`.
        """

        satellite_itrs = np.array([
            self.position.x.to(u.km).data,
            self.position.y.to(u.km).data,
            self.position.z.to(u.km).data
        ]).T
        ground_to_satellite_itrs = np.array([
            self.position.x.to(u.km).data - location.itrs_x.to(u.km).data,
            self.position.y.to(u.km).data - location.itrs_y.to(u.km).data,
            self.position.z.to(u.km).data - location.itrs_z.to(u.km).data
        ]).T
        angles = _get_ang(ground_to_satellite_itrs, satellite_itrs)

        return Quantity(angles, u.deg)

    def altitude(self) -> Quantity:
        """Return geodetic altitude.

        This method uses the WGS84 reference ellipsoid to calculate the
        geodetic altitude above the Earth's surface.

        Returns
        -------
        Quantity
            The satellite_old's geodetic altitude.

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

        return Quantity(_altitude(self.position), u.km)

    def distance(self, location: GroundLocation) -> Quantity:
        """Return distance to a ground location.

        Parameters
        ----------
        location : GroundLocation

        Returns
        -------
        Quantity
            The satellite_old's distance to `location`.
        """

        distance = np.linalg.norm(np.array([
            self.position.x.to(u.km).data - location.itrs_x.to(u.km).data,
            self.position.y.to(u.km).data - location.itrs_y.to(u.km).data,
            self.position.z.to(u.km).data - location.itrs_z.to(u.km).data
        ]), axis=0)

        return Quantity(distance, u.km)

    def save_text_file(self, file_name: str) -> None:
        gcrs_position = _itrs_to_gcrs(self.position)
        gcrs_velocity = _itrs_to_gcrs(self.velocity)

        header = "Satellite position and velocity."
        data = [
            ["Time", self.position.time],
            ["ITRS X", self.position.x],
            ["ITRS Y", self.position.y],
            ["ITRS Z", self.position.z],
            ["ITRS VX", self.velocity.x],
            ["ITRS VY", self.velocity.y],
            ["ITRS VZ", self.velocity.z],
            ["GCRS X", gcrs_position.x],
            ["GCRS Y", gcrs_position.y],
            ["GCRS Z", gcrs_position.z],
            ["GCRS VX", gcrs_velocity.x],
            ["GCRS VY", gcrs_velocity.y],
            ["GCRS VZ", gcrs_velocity.z]
        ]
        _save_data_as_txt(file_name, header, data=data)
