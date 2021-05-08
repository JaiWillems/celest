"""Ground Position Instantiation.

Notes
-----
This module contains the GroundPosition class to instantiate ground position
data for use in the Satellite class.

Class
-----
GroundPosition: Object used to combine ground position data.
    getRadius : Used to instantiate the radius attribute.
"""


from math import sin, cos, sqrt, radians
import numpy as np


class GroundPosition:
    """
    GoundPosition object stores ground location data. Used in the Satellite
    class for altitude, azimuth, and nadir angle calculations.

    Methods
    -------
    getRadius : Used to instantiate the radius attribute.

    Instance Variables
    ------------------
    name : Tag for ground station location.
    coor : Ground position coordinates as (latitude, longitude) in degrees.
    radius : Radius of earths surface at given coordinates.
    ECEFpos : GroundPosition location in ECEF frame.
    alt : Altitude data.
    az : Azimuth data.
    nadirAng : Nadir angle data.
    length : Length of data attributes.
    """

    def __init__(self, name, coor):
        """
        Define instance variables.

        Parameters
        ----------
        name : str
            A tag given to the location for indexing purposes.
        coor : tuple of floats
            Specifies observer position in decimal degrees,
            (latitude, longitude).
        """
        self.name = name
        self.coor = coor
        self.radius = self.getRadius(coor)
        self.ECEFpos = self.getECEF(coor, self.radius)
        self.alt = None
        self.az = None
        self.nadirAng = None
        self.length = None

    def __str__(self):
        """Defines GroundPosition information string."""
        title = 'Celest.GoundPosition Object\n'
        name = f'Name: {self.name}\t'
        coor = f'Coordinates: {self.coor}\t'
        radius = f'Radius: {self.radius}'

        return title + name + coor + radius

    def getRadius(self, obsCoor):
        """
        Instantiates radius attribute.

        Uses the World Geodetic System WGS84 to calculate the radius at the
        given coordinates.

        Parameters
        ----------
        obsCoor : tuple of floats
            Specifies observer position in decimal degrees,
            (latitude, longitude).

        Returns
        -------
        radius : float
            The radius given in km.
        """
        phi = radians(obsCoor[0])

        # Define WGS84 Parameters.
        semiMajor = 6378.137**2
        semiMinor = 6356.752314245**2

        numerator = semiMajor * semiMinor
        denominator = semiMajor * sin(phi)**2 + semiMinor * cos(phi)**2

        return sqrt(numerator / denominator)

    def getECEF(self, obsCoor, radius):
        """
        Instantiates ECEFpos attribute.

        Converts the ground positions geographical coordinates and radius into
        the ECEF cartesian reference frame.

        Parameters
        ----------
        obsCoor : tuple of floats
            Specifies observer position in decimal degrees,
            (latitude, longitude).
        radius : float
            Earth radius at the given coordinates.

        Returns
        -------
        ECEFpos : ndarray
            Array of shape (3,) with X, Y, Z ECEF position data.
        """
        if obsCoor[1] < 0:
            theta = radians(360 + obsCoor[1])
        else:
            theta = radians(obsCoor[1])
        phi = np.radians(90 - obsCoor[0])
        x = radius*cos(theta)*sin(phi)
        y = radius*sin(theta)*sin(phi)
        z = radius*cos(phi)
        ECEFpos = np.array([x, y, z])

        return ECEFpos
