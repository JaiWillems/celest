"""Satellite orbital representations and coordinate conversions.

The satellite module contains the Satellite class to calculate different
orbital position representations including ECI, ECEF, and Horizontal systems.
The Satellite class also forms the bases for the Encounter module to perform
encounter planning.

Notes
-----
Units of the satellite module are represented by the metric system. Specific
units will be detailed in method documentation strings.
"""


import numpy as np
import pandas as pd
from datetime import datetime
import julian
from typing import Union, Tuple, Literal
from celest.groundposition import GroundPosition


class Satellite(object):
    """Store position representations and compute orbital conversions,

    The Satellite class represents a satellite object, be it artificial or
    natural, and allows for the position to be represented with time through
    multiple representations.

    Attributes
    ----------
    times : np.array
        Array of shape (n,) of orbital time dependencies in UTC.
    ERAdata : np.array
        Array of shape (n,) of Earth rotation angles.
    ECIdata : np.array
        Array of shape (n,3) of Earth centered inertial representations where
        the columns are X, Y, Z data.
    ECEFdata : np.array
        Array of shape (n,3) of Earth centered Earth fixed representations
        where the columns are X, Y, Z data.
    gs : Dict
        Dictionary with GroundPosition.name as the key and the GroundPosition
        object as the corresponding value.
    length : int
        Length n of the data attributes.

    Methods
    -------
    time_data(timeData)
        Instantiate time attribute with orbital time dependency.
    position_data(posData, type)
        Instantiate ECIdata or ECEFdata attribute with orbital position data.
    ERA(**kwargs)
        Instantiate ERAdata attribute.
    ECI(**kwargs)
        Instantiate ECIdata attribute.
    ECEF(**kwargs)
        Instantiate ECEFdata attribute.
    horizontal(groundPos, **kwargs)
        Instantiates the GroundPosition object's .alt and .az attributes.
    nadir_ang(groundPos, **kwargs)
        Instantiates the GroundPosition object's .nadirAng attribute.
    distance(groundPos)
        Instantiates the GroundPosition object's .distance attribute.
    save_data(fileName, delimiter)
        Save class data in local directory.
    """

    def __init__(self) -> None:
        """Define instance variables."""

        self.times = None
        self.ERAdata = None
        self.ECIdata = None
        self.ECEFdata = None
        self.elevation = None
        self.gs = {}
        self.length = None

    def time_data(self, timeData: np.array, jul: bool=False, julOffset:
                  float=0) -> None:
        """Instantiate time attribute with orbital time dependency.

        Parameters
        ----------
        timeData : np.array
            Array of shape (n,) containing time data.
        jul : bool, optional
            Indicates the input times are Julian dates if True or UTC datetime
            strings if False.
        julOffset : float, optional
            Offset to be added to the inputed Julian dates.

        Notes
        -----
        For more coincise code, the timeData can be passed directly into all
        conversion methods.

        Examples
        --------
        >>> finch = Satellite()
        >>> UTCTimeData = np.array(['2020-06-01 12:00:00.0340', ...,
        ...                         '2020-06-01 12:01:00.0340'])
        >>> finch.time_data(timeData=UTCtimeData)
        """
        self.length = timeData.shape[0]

        if jul:
            self.times = timeData + julOffset
        else:
            self.times = np.zeros((self.length,))
            for i in range(self.length):
                try:
                    self.times[i] = julian.to_jd(datetime.strptime(
                        timeData[i], "%Y-%m-%d %H:%M:%S.%f"))
                except:
                    self.times[i] = julian.to_jd(datetime.strptime(
                        timeData[i], "%Y-%m-%d %H:%M:%S"))

    def position_data(self, posData: np.array, type: Literal["ECI", "ECEF"]) -> None:
        """Instantiate orbital position data.

        Parameters
        ----------
        posData : np.array
            Array of shape (n,3) with columns of X, Y, Z position data. Units
            in meters.
        type : {"ECI", "ECEF"}
            Specifies posData as either "ECI" or "ECEF".

        Notes
        -----
        For more coincise code, the posData can be passed directly into all
        conversion methods.

        Examples
        --------
        >>> finch = Satellite()
        >>> ECIvec = np.array([[-4.46e+03, -5.22e+03, 1.75e-04], ...,
        ...                    [2.73e+03, 2.08e+03, -6.02e+03]])
        >>> finch.position_data(posData=ECIvec, type="ECI")
        """
        if type == 'ECI':
            self.ECIdata = posData
        elif type == 'ECEF':
            self.ECEFdata = posData

    def ERA(self, timeData: np.array=None, **kwargs:
            Union[bool, float]) -> np.array:
        """Instantiate ERAdata attribute.

        Parameters
        ----------
        timeData : np.array, optional
            Array of shape (n,) containing time data.
        **kwargs : dict, optional
            Extra arguments to `getERA`: refer to getERA documentation for a
            list of all possible arguments.

        Returns
        -------
        np.array
            Array of shape (n,) containing radian earth rotation angles.

        Notes
        -----
        The Satellite class instance must have the times attribute initiated or
        a timeData input passed in under \*\*kwargs.

        The method implements the earth rotation angle formula:

        .. math:: \gamma = 360.99(\Delta T) + 280.46

        Examples
        --------
        >>> finch = Satellite()
        >>> UTCTimeData = np.array(['2020-06-01 12:00:00.0340', ...,
        ...                         '2020-06-01 12:01:00.0340'])
        >>> ERAangles = finch.ERA(timeData=UTCtimeData)
        """
        if type(self.times) == type(None):
            self.time_data(timeData, ** kwargs)

        angArr = np.zeros((self.length,))
        Julian = self.times

        # Multiply time elapsed since J2000 by ERA and add J2000 orientation.
        dJulian = Julian - 2451545
        angArr = (360.9856123035484 * dJulian + 280.46) % 360

        angArr = np.radians(angArr)
        self.ERAdata = angArr

        return self.ERAdata

    def ECI(self, posData: np.array=None, timeData: np.array=None, **kwargs:
            Union[bool, float]) -> np.array:
        """Instantiate ECIdata attribute.

        Parameters
        ----------
        posData : np.array, optional
            Array of shape (n,3) with columns of X, Y, Z ECEF position data.
        timeData : np.array, optional
            Array of shape (n,) containing time data.
        **kwargs : dict, optional
            Extra arguments to `getECI`: refer to getECI documentation for a
            list of all possible arguments.

        Returns
        -------
        np.array
            Array of shape (n,3) with columns of X, Y, Z ECI position data.

        See Also
        --------
        getECEF : Initiate ECEFdata attribute.

        Notes
        -----
        The Satellite class instance must have ECEFdata and ERAdata attributes
        initiated or posData and timeData inputs passed in.

        Examples
        --------
        >>> UTCTimeData = np.array(['2020-06-01 12:00:00.0340', ...,
        ...                         '2020-06-01 12:01:00.0340'])
        >>> ECEFvec = np.array([[-4.46e+03, -5.22e+03, 1.75e-04], ...,
        ...                     [2.73e+03, 2.08e+03, -6.02e+03]])
        >>> finch = Satellite()
        >>> ECIvec = finch.ECI(posData=ECEFvec, timeData=UTCTimeData)
        """
        if type(self.ERAdata) == type(None):
            self.ERA(timeData=timeData, **kwargs)
        if type(self.ECEFdata) == type(None):
            self.position_data(posData=posData, type="ECEF")

        # Rotate ECEFdata around z-axis by ERA.
        ECIvec = np.zeros((self.length, 3))
        A11 = np.cos(self.ERAdata)
        A12 = -np.sin(self.ERAdata)
        A21 = np.sin(self.ERAdata)
        A22 = np.cos(self.ERAdata)

        ECIvec[:, 0] = np.add(np.multiply(
            A11, self.ECEFdata[:, 0]), np.multiply(A12, self.ECEFdata[:, 1]))
        ECIvec[:, 1] = np.add(np.multiply(
            A21, self.ECEFdata[:, 0]), np.multiply(A22, self.ECEFdata[:, 1]))
        ECIvec[:, 2] = self.ECEFdata[:, 2]

        self.ECIdata = ECIvec

        return self.ECIdata

    def ECEF(self, posData: np.array=None, timeData: np.array=None, **kwargs:
             Union[bool, float]) -> np.array:
        """Instantiate ECEFdata attribute.

        Parameters
        ----------
        posData : np.array, optional
            Array of shape (n,3) with columns of X, Y, Z ECI position data.
        timeData : np.array, optional
            Array of shape (n,) containing time data.
        **kwargs : dict, optional
            Extra arguments to `getECEF`: refer to getECEF documentation for a
            list of all possible arguments.

        Returns
        -------
        np.array
            Array of shape (n,3) with columns of X, Y, Z ECEF position data.

        See Also
        --------
        getECI : Initiate ECIdata attribute.

        Notes
        -----
        The Satellite class instance must have ECIdata and ERAdata attributes
        initiated or posData and timeData inputs passed in.

        Examples
        --------
        >>> UTCTimeData = np.array(['2020-06-01 12:00:00.0340', ...,
        ...                         '2020-06-01 12:01:00.0340'])
        >>> ECIvec = np.array([[-4.46e+03, -5.22e+03, 1.75e-04], ...,
        ...                    [2.73e+03, 2.08e+03 -6.02e+03]])
        >>> finch = Satellite()
        >>> ECEFvec = finch.ECEF(posData=ECIvec, timeData=UTCTimeData)
        """
        if type(self.ERAdata) == type(None):
            self.ERA(timeData=timeData, **kwargs)
        if type(self.ECIdata) == type(None):
            self.position_data(posData=posData, type="ECI")

        # Rotate ECIdata around z-axis by -ERA.
        ECEFvec = np.zeros((self.length, 3))
        A11 = np.cos(-self.ERAdata)
        A12 = -np.sin(-self.ERAdata)
        A21 = np.sin(-self.ERAdata)
        A22 = np.cos(-self.ERAdata)

        ECEFvec[:, 0] = np.add(np.multiply(
            A11, self.ECIdata[:, 0]), np.multiply(A12, self.ECIdata[:, 1]))
        ECEFvec[:, 1] = np.add(np.multiply(
            A21, self.ECIdata[:, 0]), np.multiply(A22, self.ECIdata[:, 1]))
        ECEFvec[:, 2] = self.ECIdata[:, 2]

        self.ECEFdata = ECEFvec

        return self.ECEFdata

    def _get_ang(self, vecOne: np.array, vecTwo: np.array) -> float:
        """Calculate degree angle bewteen two vectors.
        
        Parameters
        ----------
        vecOne, vecTwo : np.array
            Arrays of shape (n,3) with rows of ECEF data.
        
        Returns
        -------
        float
            Degree angle between the two arrays.
        """
        # Use simple linalg formula.
        dividend = np.einsum('ij, ij->i', vecOne, vecTwo)
        divisor = np.multiply(np.linalg.norm(
            vecOne, axis=1), np.linalg.norm(vecTwo, axis=1))
        arg = np.divide(dividend, divisor)
        ang = np.degrees(np.arccos(arg))

        return ang

    def _geo_to_ECEF(self, obsCoor: Tuple[float, float], radius:
                     float) -> np.array:
        """Convert geographical coordinates to ECEF.
        
        Parameters
        ----------
        obsCoor : Tuple
            Coordinates of a ground location in decimal degrees
            `(lattitude, longitude)`.
        radius : float
            Radius of the Earth at position `obsCoor`.
        
        Returns
        -------
        np.array
            Array of shape (n,3) with rows of ECEF data.
        """
        if obsCoor[1] < 0:
            theta = np.radians(360 + obsCoor[1])
        else:
            theta = np.radians(obsCoor[1])
        phi = np.radians(90 - obsCoor[0])
        x = radius*np.cos(theta)*np.sin(phi)
        y = radius*np.sin(theta)*np.sin(phi)
        z = radius*np.cos(phi)

        return np.array([x, y, z])

    def horizontal(self, groundPos: GroundPosition, posData: np.array=None,
                   timeData: np.array=None, **kwargs:
                   Union[bool, float]) -> np.array:
        """Instantiate GroundPositions .alt and .az attributes.

        This method takes in a GroundPosition object and instantiates its .alt
        and .az attributes using the Satellite instance's position data. The
        GroundPosition object is then stored as a value in the gs attribute's
        dictionary with the corresponding key being its name attribute.

        Parameters
        ----------
        groundPos : GroundPosition object
            GroundPosition object instantiated as per its documentation.
        posData : np.array, optional
            Array of shape (n,3) with columns of X, Y, Z position data
            assumed ECI.
        timeData : np.array, optional
            Array of shape (n,) containing time data.
        **kwargs : dict, optional
            Extra arguments to `horizontal`: refer to horizontal documentation
            for a list of all possible arguments.

        Returns
        -------
        tuple
            The altitude/azimuth data is returned in a tuple where both the
            altitude and azimuth are arrays of shape (n,).

        Notes
        -----
        The Satellite class instance must have the ECEFdata attribute initiated
        or have the posData and timeData inputs passed in.

        Examples
        --------
        >>> UTCTimeData = np.array(['2020-06-01 12:00:00.0340', ...,
        ...                         '2020-06-01 12:01:00.0340'])
        >>> ECIvec = np.array([[-4.46e+03, -5.22e+03, 1.75e-04], ...,
        ...                    [2.73e+03, 2.08e+03, -6.02e+03]])
        >>> toronto = GroundPosition(name="Toronto",
        ...                          coor=(43.662300, -79.394530))
        >>> finch = Satellite()
        >>> Alt, Az = finch.horizontal(groundPos=toronto, posData=ECIvec,
        ...                          timeData=UTCTimeData)
        """
        if type(self.ECEFdata) == type(None):
            self.ECEF(posData=posData, timeData=timeData, **kwargs)

        if groundPos.name not in self.gs:
            groundPos.length = self.length
            self.gs[groundPos.name] = groundPos

        # Convert observer position into spherical then cartesian.
        obsCoor = groundPos.coor
        radius = groundPos.radius
        xyzObs = np.full((self.length, 3), self._geo_to_ECEF(obsCoor, radius))

        # Determine line of sight vector then altitude.
        ECEFvec = self.ECEFdata
        xyzLOS = np.subtract(ECEFvec, xyzObs)
        self.gs[groundPos.name].alt = 90 - self._get_ang(xyzLOS, xyzObs)

        # Find surface tangent vector passing through z-axis.
        kHat = np.full((self.length, 3), np.array([0, 0, 1]))
        beta = np.pi/2 - np.radians(90 - obsCoor[0])
        tangentVec = np.subtract((kHat.T * radius/np.sin(beta)).T, xyzObs)

        # Find LOS projection on tangent plane.
        coeff = np.einsum('ij, ij->i', xyzLOS, xyzObs)/radius**2
        normProj = (xyzObs.T * coeff).T
        projLOS = np.subtract(xyzLOS, normProj)

        # Determing azimuth.
        vecOne = np.cross(tangentVec, xyzObs)
        normOne = 1/np.linalg.norm(vecOne, axis=1).reshape((self.length, 1))
        vecOneUnit = normOne*vecOne
        vecTwo = (vecOneUnit.T * np.einsum('ij, ij->i', projLOS, vecOneUnit)).T
        normTwo = 1/np.linalg.norm(vecTwo, axis=1).reshape((self.length, 1))
        vecTwoUnit = normTwo*vecTwo

        posEq = np.isclose(vecOneUnit, vecTwoUnit, atol=0.01)
        negEq = np.isclose(vecOneUnit, -vecTwoUnit, atol=0.01)
        posInd = np.where(np.all(posEq, axis=1))[0]
        negInd = np.where(np.all(negEq, axis=1))[0]

        Az = np.zeros((self.length,))
        Az[posInd] = self._get_ang(tangentVec[posInd], projLOS[posInd])
        Az[negInd] = 360 - self._get_ang(tangentVec[negInd], projLOS[negInd])

        self.gs[groundPos.name].az = Az

        return self.gs[groundPos.name].alt, self.gs[groundPos.name].az

    def nadir_ang(self, groundPos: GroundPosition, posData: np.array=None,
                  timeData: np.array=None, **kwargs:
                  Union[bool, float]) -> np.array:
        """Instantiate GroundPositions.nadirAng attribute.

        This method takes in a GroundPosition object and instantiates its
        nadirAng attribute with the angle made between the satellites nadir and
        line-of-sight vectors. The GroundPosition object is then stored as a
        value in the gs attribute dictionary with the corresponding key being
        its name attribute.

        Parameters
        ----------
        groundPos : GroundPosition object
            GroundPosition object instantiated as per its documentation.
        posData : np.array, optional
            Array of shape (n,3) with columns of X, Y, Z position data
            assumed ECI.
        timeData : np.array, optional
            Array of shape (n,) containing time data.
        **kwargs : dict, optional
            Extra arguments to `nadir_ang`: refer to nadir_ang documentation
            for a list of all possible arguments.

        Returns
        -------
        np.array
            Array of shape (n,) with nadir-LOS angle data.

        Notes
        -----
        The Satellite instance must have the ECEFdata attribute initiated or
        the posData and timeData inputs passed in.

        Examples
        --------
        >>> UTCTimeData = np.array(['2020-06-01 12:00:00.0340', ...,
        ...                         '2020-06-01 12:01:00.0340'])
        >>> ECIvec = np.array([[-4.46e+03, -5.22e+03, 1.75e-04], ...,
        ...                    [2.73e+03, 2.08e+03 -6.02e+03]])
        >>> toronto = GroundPosition(name="Toronto",
        ...                          coor=(43.662300, -79.394530))
        >>> finch = Satellite()
        >>> NdrAng = finch.nadir_ang(groundPos=toronto, posData=ECIvec,
        ...                          timeData=UTCTimeData)
        """
        if type(self.ECEFdata) == type(None):
            self.ECEF(posData=posData, timeData=timeData, **kwargs)

        if groundPos.name not in self.gs:
            groundPos.length = self.length
            self.gs[groundPos.name] = groundPos

        obsCoor = groundPos.coor
        radius = groundPos.radius

        xyzObs = np.full((self.length, 3), self._geo_to_ECEF(obsCoor, radius))

        ECEFvec = self.ECEFdata
        xyzLOS = np.subtract(ECEFvec, xyzObs)

        ang = self._get_ang(xyzLOS, ECEFvec)
        self.gs[groundPos.name].nadirAng = ang

        return self.gs[groundPos.name].nadirAng

    def distance(self, groundPos: GroundPosition) -> np.array:
        """Get distance from satellite to ground location.

        Parameters
        ----------
        groundPos : GroundPosition
            GroundPosition object instantiated as per its documentation.

        Returns
        -------
        np.array
            Array of shape (n,) containing time varying distances between the
            satellite and a ground location.
        """
        if groundPos.name not in self.gs:
            groundPos.length = self.length
            self.gs[groundPos.name] = groundPos

        n = self.length
        satECEF = self.ECEFdata
        groundECEF = np.full((n, 3), groundPos.ECEFpos)

        # Find LOS vector norm.
        LOSvec = satECEF - groundECEF
        distances = np.linalg.norm(LOSvec, axis=1)

        self.gs[groundPos.name].distance = distances

        return distances

    def _WGS84_radius(self, lattitude: float) -> float:
        """Calculate radius from WGS84.

        Parameters
        ----------
        latitude : float
            Lattitude of a ground location.
        
        Returns
        -------
        float
            Earth's radius at `lattitude` using WGS84.
        """
        phi = np.radians(lattitude)

        # Define WGS84 Parameters.
        semiMajor = 6378.137**2
        semiMinor = 6356.752314245**2

        numerator = semiMajor * semiMinor
        denominator = semiMajor * np.sin(phi)**2 + semiMinor * np.cos(phi)**2

        return np.sqrt(numerator / denominator)

    def altitude(self) -> np.array:
        """Get the altitude/elevation of the satellite.

        This method implements WGS84 to calculate the radius of the earth and
        then computes the satellites altitude/elevation.

        Returns
        -------
        np.array
            Array of shape (n,) containing time varying altitudes of the
            satellite.
        """
        ECEFdata = self.ECEFdata
        x = ECEFdata[:, 0]
        y = ECEFdata[:, 1]
        z = ECEFdata[:, 2]
        arg = np.sqrt(x**2 + y**2) / z
        lattitude = np.arctan(arg)

        earthRadius = self._WGS84_radius(lattitude)
        satRadius = np.linalg.norm(ECEFdata, axis=1)

        altitude = satRadius - earthRadius
        self.elevation = altitude

        return altitude

    def save_data(self, fileName: str, delimiter: str) -> None:
        """Save satellite data to local directory.

        Parameters
        ----------
        fileName : str
            File name of the output file as wither a .txt or .csv file.
        delimiter : str
            String of length 1 representing the feild delimiter for the output
            file.

        Notes
        -----
        It is recommended to use a tab delimiter for .txt files and comma
        delimiters for .csv files. The method will return an error if the
        fileName already exists in the current working directory.

        Examples
        --------
        >>> finch.save_data(fileName="data.csv", delimiter=",")
        """
        data = {}

        if type(self.times) != type(None):
            data['Time (UTC)'] = pd.Series(self.times)

        if type(self.ERAdata) != type(None):
            data['ERA (rad)'] = pd.Series(self.ERAdata)

        if type(self.ECIdata) != type(None):
            data['ECI.X'] = pd.Series(self.ECIdata[:, 0])
            data['ECI.Y'] = pd.Series(self.ECIdata[:, 1])
            data['ECI.Z'] = pd.Series(self.ECIdata[:, 2])

        if type(self.ECEFdata) != type(None):
            data['ECEF.X'] = pd.Series(self.ECEFdata[:, 0])
            data['ECEF.Y'] = pd.Series(self.ECEFdata[:, 1])
            data['ECEF.Z'] = pd.Series(self.ECEFdata[:, 2])

        for key in self.gs:
            if type(self.gs[key].alt) != type(None):
                data[f'{key}.alt'] = self.gs[key].alt
                data[f'{key}.az'] = self.gs[key].az
            if type(self.gs[key].nadirAng) != type(None):
                data[f'{key}.NdrAng'] = self.gs[key].nadirAng

        df = pd.DataFrame(data)
        df.to_csv(fileName, sep=delimiter)
