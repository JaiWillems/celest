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
from math import pi, cos, sin, radians
from datetime import datetime
import julian
from groundposition import GroundPosition


class Satellite:
    """Store position representations and compute orbital conversions,

    The Satellite class represents a satellite object, be it artificial or
    natural, and allows for the position to be represented with time through
    multiple representations.

    Parameters
    ----------
    None

    Attributes
    ----------
    times : np.ndarray
        Array of shape (n,) of orbital time dependencies in UTC.
    ERAdata : np.ndarray
        Array of shape (n,) of Earth rotation angles.
    ECIdata : np.ndarray
        Array of shape (n,3) of Earth centered inertial representations where
        the columns are X, Y, Z data.
    ECEFdata : np.ndarray
        Array of shape (n,3) of Earth centered Earth fixed representations
        where the columns are X, Y, Z data.
    gs : Dict
        Dictionary with GroundPosition.name as the key and the GroundPosition
        object as the corresponding value.
    length : int
        Length n of the data attributes.

    Methods
    -------
    timeData(timeData)
        Instantiate time attribute with orbital time dependency.
    positionData(posData, type)
        Instantiate ECIdata or ECEFdata attribute with orbital position data.
    getERA(**kwargs)
        Instantiate ERAdata attribute.
    getECI(**kwargs)
        Instantiate ECIdata attribute.
    getECEF(**kwargs)
        Instantiate ECEFdata attribute.
    getAltAz(groundPos, **kwargs)
        Instantiates the GroundPosition object's .alt and .az attributes.
    getNdrAng(groundPos, **kwargs)
        Instantiates the GroundPosition object's .nadirAng attribute.
    saveData(fileName, delimiter)
        Save class data in local directory.
    """

    def __init__(self) -> None:
        """Define instance variables."""

        self.times = None
        self.ERAdata = None
        self.ECIdata = None
        self.ECEFdata = None
        self.gs = {}
        self.length = None

    def timeData(self, timeData: np.ndarray) -> None:
        """Instantiate time attribute with orbital time dependency.

        Parameters
        ----------
        timeData : np.ndarray
            Array of shape (n,) containing datetime objects in UTC.

        Returns
        -------
        None

        Notes
        -----
        For more coincise code, the timeData can be passed directly into all
        conversion methods.

        Examples
        --------
        >>> finch = Satellite()
        >>> UTCTimeData = np.array(['2020-06-01 12:00:00.0340', ...,
        ...                         '2020-06-01 12:01:00.0340'])
        >>> finch.timeData(timeData=UTCtimeData)
        """
        self.times = timeData
        self.length = timeData.shape[0]

    def positionData(self, posData: np.ndarray, type: str) -> None:
        """Instantiate orbital position data.

        Parameters
        ----------
        posData : np.ndarray
            Array of shape (n,3) with columns of X, Y, Z position data. Units
            in meters.
        type : {"ECI", "ECEF"}
            Specifies posData as either "ECI" or "ECEF".

        Returns
        -------
        None

        Notes
        -----
        For more coincise code, the posData can be passed directly into all
        conversion methods.

        Examples
        --------
        >>> finch = Satellite()
        >>> ECIvec = np.array([[-4.46e+03, -5.22e+03, 1.75e-04], ...,
        ...                    [2.73e+03, 2.08e+03, -6.02e+03]])
        >>> finch.positionData(posData=ECIvec, type="ECI")
        """
        if type == 'ECI':
            self.ECIdata = posData
        elif type == 'ECEF':
            self.ECEFdata = posData

    def getERA(self, **kwargs: np.ndarray) -> np.ndarray:
        """Instantiate ERAdata attribute.

        Parameters
        ----------
        **kwargs : dict, optional
            Extra arguments to `getERA`: refer to getERA documentation for a
            list of all possible arguments.

        Returns
        -------
        np.ndarray
            Array of shape (n,) containing radian earth rotation angles.

        Notes
        -----
        The Satellite class instance must have the times attribute initiated or
        a timeData input passed in under **kwargs.

        The method implements the earth rotation angle formula:

        .. math:: \gamma = 360.99(\Delta T) + 280.46

        Examples
        --------
        >>> finch = Satellite()
        >>> UTCTimeData = np.array(['2020-06-01 12:00:00.0340', ...,
        ...                         '2020-06-01 12:01:00.0340'])
        >>> ERAangles = finch.getERA(timeData=UTCtimeData)
        """
        if type(self.times) == type(None):
            timeData = kwargs['timeData']
            self.timeData(timeData)

        angArr = np.zeros((self.length,))
        Julian = np.zeros((self.length,))

        for i in range(self.length):
            Julian[i] = julian.to_jd(datetime.strptime(self.times[i][:19],
                                                       '%Y-%m-%d %H:%M:%S'))

        # Multiply the elaped time since J2000 by Earth rotation rate and add
        # orientation at J2000.
        dJulian = Julian - 2451545
        angArr = (360.9856123035484 * dJulian + 280.46) % 360

        angArr = np.radians(angArr)
        self.ERAdata = angArr

        return self.ERAdata

    def getECI(self, **kwargs: np.ndarray) -> np.ndarray:
        """Instantiate ECIdata attribute.

        Parameters
        ----------
        **kwargs : dict, optional
            Extra arguments to `getECI`: refer to getECI documentation for a
            list of all possible arguments.

        Returns
        -------
        np.ndarray
            Array of shape (n,3) with columns of X, Y, Z ECI position data.

        See Also
        --------
        getECEF : Initiate ECEFdata attribute.

        Notes
        -----
        The Satellite class instance must have ECEFdata and ERAdata attributes
        initiated or posData and timeData inputs passed in under **kwargs.

        Examples
        --------
        >>> UTCTimeData = np.array(['2020-06-01 12:00:00.0340', ...,
        ...                         '2020-06-01 12:01:00.0340'])
        >>> ECEFvec = np.array([[-4.46e+03, -5.22e+03, 1.75e-04], ...,
        ...                     [2.73e+03, 2.08e+03, -6.02e+03]])
        >>> finch = Satellite()
        >>> ECIvec = finch.getECI(posData=ECEFvec, timeData=UTCTimeData)
        """
        if type(self.ERAdata) == type(None):
            self.getERA(**kwargs)
        if type(self.ECEFdata) == type(None):
            ECEFdata = kwargs['posData']
            self.positionData(ECEFdata, 'ECEF')

        # Rotate ECEFdata around z-axis by ERA.
        ECIvec = np.zeros((self.length, 3))
        A11 = np.cos(self.ERAdata)
        A12 = -np.sin(self.ERAdata)
        A21 = np.sin(self.ERAdata)
        A22 = np.cos(self.ERAdata)

        ECIvec[:, 0] = np.add(np.multiply(A11, self.ECEFdata[:, 0]),
                              np.multiply(A12, self.ECEFdata[:, 1]))
        ECIvec[:, 1] = np.add(np.multiply(A21, self.ECEFdata[:, 0]),
                              np.multiply(A22, self.ECEFdata[:, 1]))
        ECIvec[:, 2] = self.ECEFdata[:, 2]

        self.ECIdata = ECIvec

        return self.ECIdata

    def getECEF(self, **kwargs: np.ndarray) -> np.ndarray:
        """Instantiate ECEFdata attribute.

        Parameters
        ----------
        **kwargs : dict, optional
            Extra arguments to `getECEF`: refer to getECEF documentation for a
            list of all possible arguments.

        Returns
        -------
        np.ndarray
            Array of shape (n,3) with columns of X, Y, Z ECEF position data.

        See Also
        --------
        getECI : Initiate ECIdata attribute.

        Notes
        -----
        The Satellite class instance must have ECIdata and ERAdata attributes
        initiated or posData and timeData inputs passed in under **kwargs.

        Examples
        --------
        >>> UTCTimeData = np.array(['2020-06-01 12:00:00.0340', ...,
        ...                         '2020-06-01 12:01:00.0340'])
        >>> ECIvec = np.array([[-4.46e+03, -5.22e+03, 1.75e-04], ...,
        ...                    [2.73e+03, 2.08e+03 -6.02e+03]])
        >>> finch = Satellite()
        >>> ECEFvec = finch.getECEF(posData=ECIvec, timeData=UTCTimeData)
        """
        if type(self.ERAdata) == type(None):
            self.getERA(**kwargs)
        if type(self.ECIdata) == type(None):
            ECIdata = kwargs['posData']
            self.positionData(ECIdata, 'ECI')

        # Rotate ECIdata around z-axis by -ERA.
        ECEFvec = np.zeros((self.length, 3))
        A11 = np.cos(-self.ERAdata)
        A12 = -np.sin(-self.ERAdata)
        A21 = np.sin(-self.ERAdata)
        A22 = np.cos(-self.ERAdata)

        ECEFvec[:, 0] = np.add(np.multiply(A11, self.ECIdata[:, 0]),
                               np.multiply(A12, self.ECIdata[:, 1]))
        ECEFvec[:, 1] = np.add(np.multiply(A21, self.ECIdata[:, 0]),
                               np.multiply(A22, self.ECIdata[:, 1]))
        ECEFvec[:, 2] = self.ECIdata[:, 2]

        self.ECEFdata = ECEFvec

        return self.ECEFdata

    def getAltAz(self, groundPos: GroundPosition, **kwargs: np.ndarray) -> np.ndarray:
        """Instantiate GroundPositions .alt and .az attributes.

        This method takes in a GroundPosition object and instantiates its .alt
        and .az attributes using the Satellite instance's position data. The
        GroundPosition object is then stored as a value in the gs attribute's
        dictionary with the corresponding key being its name attribute.

        Parameters
        ----------
        groundPos : GroundPosition object
            GroundPosition object instantiated as per its documentation.
        **kwargs : dict, optional
            Extra arguments to `getAltAz`: refer to getAltAz documentation for
            a list of all possible arguments.

        Returns
        -------
        tuple
            The altitude/azimuth data is returned in a tuple where both the
            altitude and azimuth are arrays of shape (n,).

        Notes
        -----
        The Satellite class instance must have the ECEFdata attribute initiated
        or have the posData and timeData inputs passed in under **kwargs.

        Examples
        --------
        >>> UTCTimeData = np.array(['2020-06-01 12:00:00.0340', ...,
        ...                         '2020-06-01 12:01:00.0340'])
        >>> ECIvec = np.array([[-4.46e+03, -5.22e+03, 1.75e-04], ...,
        ...                    [2.73e+03, 2.08e+03, -6.02e+03]])
        >>> toronto = GroundPosition(name="Toronto",
        ...                          coor=(43.662300, -79.394530))
        >>> finch = Satellite()
        >>> Alt, Az = finch.getAltAz(groundPos=toronto, posData=ECIvec,
        ...                          timeData=UTCTimeData)
        """
        def getAngle(vecOne: np.ndarray, vecTwo: np.ndarray) -> float:
            """Calculate degree angle bewteen two vectors.

            Takes two multidimensional cartesian coordinate vectors and returns
            degree angle between the two arrays element-wise.

            Parameters
            ----------
            vecOne, vecTwo: np.ndarray
                Array of shape (n,3) representing n vectors of columns X, Y, Z.

            Returns
            -------
            float
                Degree angle between vecOne and vecTwo.
            """
            # Use simple linalg formula.
            dividend = np.einsum('ij, ij->i', vecOne, vecTwo)
            divisor = np.multiply(np.linalg.norm(vecOne, axis=1),
                                  np.linalg.norm(vecTwo, axis=1))
            arg = np.divide(dividend, divisor)
            ang = np.degrees(np.arccos(arg))

            return ang

        if type(self.ECEFdata) == type(None):
            self.getECEF(**kwargs)

        if groundPos.name not in self.gs:
            groundPos.length = self.length
            self.gs[groundPos.name] = groundPos

        # Convert observer position into spherical then cartesian.
        obsCoor = groundPos.coor
        radius = groundPos.radius
        if obsCoor[1] < 0:
            theta = radians(360 + obsCoor[1])
        else:
            theta = radians(obsCoor[1])
        phi = np.radians(90 - obsCoor[0])
        x = radius*cos(theta)*sin(phi)
        y = radius*sin(theta)*sin(phi)
        z = radius*cos(phi)
        xyzObs = np.full((self.length, 3), np.array([x, y, z]))

        # Determine line of sight vector then altitude.
        ECEFvec = self.ECEFdata
        xyzLOS = np.subtract(ECEFvec, xyzObs)
        self.gs[groundPos.name].alt = 90 - getAngle(xyzLOS, xyzObs)

        # Find surface tangent vector passing through z-axis.
        kHat = np.full((self.length, 3), np.array([0, 0, 1]))
        beta = pi/2 - phi
        tangentVec = np.subtract((kHat.T * radius/np.sin(beta)).T, xyzObs)

        # Find LOS projection on tangent plane.
        coeff = np.einsum('ij, ij->i', xyzLOS, xyzObs)/radius**2
        normProj = (xyzObs.T * coeff).T
        projLOS = np.subtract(xyzLOS, normProj)

        # Determing aECIzimuth.
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
        Az[posInd] = getAngle(tangentVec[posInd], projLOS[posInd])
        Az[negInd] = 360 - getAngle(tangentVec[negInd], projLOS[negInd])

        self.gs[groundPos.name].az = Az

        return self.gs[groundPos.name].alt, self.gs[groundPos.name].az

    def getNdrAng(self, groundPos: GroundPosition, **kwargs: np.ndarray) -> np.ndarray:
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
        **kwargs : dict, optional
            Extra arguments to `getNdrAng`: refer to getNdrAng documentation
            for a list of all possible arguments.

        Returns
        -------
        np.ndarray
            Array of shape (n,) with nadir-LOS angle data.

        Notes
        -----
        The Satellite instance must have the ECEFdata attribute initiated or
        the posData and timeData inputs passed in as **kwargs.

        Examples
        --------
        >>> UTCTimeData = np.array(['2020-06-01 12:00:00.0340', ...,
        ...                         '2020-06-01 12:01:00.0340'])
        >>> ECIvec = np.array([[-4.46e+03, -5.22e+03, 1.75e-04], ...,
        ...                    [2.73e+03, 2.08e+03 -6.02e+03]])
        >>> toronto = GroundPosition(name="Toronto",
        ...                          coor=(43.662300, -79.394530))
        >>> finch = Satellite()
        >>> NdrAng = finch.getNdrAng(groundPos=toronto, posData=ECIvec,
        ...                          timeData=UTCTimeData)
        """
        if type(self.ECEFdata) == type(None):
            self.getECEF(**kwargs)

        if groundPos.name not in self.gs:
            groundPos.length = self.length
            self.gs[groundPos.name] = groundPos

        obsCoor = groundPos.coor
        radius = groundPos.radius
        if obsCoor[1] < 0:
            theta = radians(360 + obsCoor[1])
        else:
            theta = radians(obsCoor[1])
        phi = np.radians(90 - obsCoor[0])
        x = radius*cos(theta)*sin(phi)
        y = radius*sin(theta)*sin(phi)
        z = radius*cos(phi)
        xyzObs = np.full((self.length, 3), np.array([x, y, z]))

        ECEFvec = self.ECEFdata
        xyzLOS = np.subtract(ECEFvec, xyzObs)

        # Use simple linalg formula.
        dividend = np.einsum('ij, ij->i', xyzLOS, ECEFvec)
        divisor = np.multiply(np.linalg.norm(xyzLOS, axis=1),
                              np.linalg.norm(ECEFvec, axis=1))
        arg = np.divide(dividend, divisor)
        ang = np.degrees(np.arccos(arg))

        self.gs[groundPos.name].nadirAng = ang

        return self.gs[groundPos.name].nadirAng

    def saveData(self, fileName: str, delimiter: str) -> None:
        """Save satellite data to local directory.

        Parameters
        ----------
        fileName : str
            File name of the output file as wither a .txt or .csv file.
        delimiter : str
            String of length 1 representing the feild delimiter for the output
            file.

        Returns
        -------
        None

        Notes
        -----
        It is recommended to use a tab delimiter for .txt files and comma
        delimiters for .csv files. The method will return an error if the
        fileName already exists in the current working directory.

        Examples
        --------
        >>> finch.saveData(fileName="data.csv", delimiter=",")
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
