"""Satellite orbital coordinate conversions.

Notes
-----
This module contains the Satellite class to perform coordinate
system calculations.

Class's
-------
Satellite: Object used to collect and calculate orbital representations.
    timeData : Instantiate time attribute with orbital time dependency.
    positionData : Instantiate ECIdata or ECEFdata attribute with orbital
        position data.
    getERA : Instantiate ERAdata attribute.
    getECI : Instantiate ECIdata attribute.
    getECEF : Instantiate ECEFdata attribute.
    getAltAz : Instantiate horizontal attribute.
    getNdrAng : Instantiate nadirAng attribute.
    saveData : Save class data in local directory.
"""


import numpy as np
import pandas as pd
from math import pi, cos, sin, radians
from datetime import datetime, timedelta
import julian


class Satellite:
    """
    Satellite object stores and computes orbital coordinate data. Units are in
    metric.

    Methods
    -------
    timeData : Instantiate time attribute with orbital time dependency.
    positionData : Instantiate ECIdata or ECEFdata attribute with orbital
        position data.
    getERA : Instantiate ERAdata attribute.
    getECI : Instantiate ECIdata attribute.
    getECEF : Instantiate ECEFdata attribute.
    getAltAz : Instantiate horizontal attribute.
    getNdrAng : Instantiate nadirAng attribute.
    saveData : Save class data in local directory.

    Instance Variables
    ------------------
    times : Orbital time dependencies.
    ERAdata : Earth rotation angles.
    ECIdata : Earth centered inertial position data.
    ECEFdata : Earth centered earth fixed position data.
    horizontal : Altitude and azimuth data.
    nadirAng : Nadir-LOS angle data.
    length : Length of data attributes.
    """

    def __init__(self):
        """Define instance variables."""

        self.times = None
        self.ERAdata = None
        self.ECIdata = None
        self.ECEFdata = None
        self.horizontal = None
        self.nadirAng = None
        self.length = None

    def timeData(self, timeData):
        """
        Instantiate time attribute with orbital time dependency.

        Takes UTC times.

        Parameters
        ----------
        timeData : ndarray
            Array of shape (n,). Contains datetime objects in UTC.

        Returns
        -------
        None
        """
        self.times = timeData
        self.length = timeData.shape[0]

    def positionData(self, posData, type):
        """
        Instantiate ECIdata or ECEFdata attribute with orbital position data.

        Takes in ECI or ECEF positon data in cartesian coordinates and a data
        type specifier.

        Parameters
        ----------
        posData : ndarray
            Array of shape (n,3) with columns of X, Y, Z position data.
        type : str
            Specifies posData as "ECI" or "ECEF".

        Returns
        -------
        None
        """
        if type == "ECI":
            self.ECIdata = posData
        elif type == "ECEF":
            self.ECEFdata = posData

    def getERA(self, **kwargs):
        """
        Instantiate ERAdata attribute.

        Uses times attribute to calculate the Earth Rotation Angle in radians.

        Parameters
        ----------
        None

        **kwargs
        --------
        timeData : ndarray
            Array of shape (n,). Contains datetime objects in UTC.

        Returns
        -------
        self.ERAdata : ndarray
            Array of shape (n,). Contains radian earth rotation angles.

        Usage
        -----
        Must have times attribute initiated or timeData input.
        """
        if type(self.times) == type(None):
            timeData = kwargs["timeData"]
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

    def getECI(self, **kwargs):
        """
        Instantiate ECIdata attribute.

        Computes ECIdata attribute from ECEFdata and ERAdata attribites.

        Parameters
        ----------
        None

        **kwargs
        --------
        posData : ndarray
            Array of shape (n,3) with columns of X, Y, Z position data assumed
            ECEF.
        timeData : ndarray
            Array of shape (n,). Contains datetime objects in UTC.

        Returns
        -------
        self.ECIdata : ndarray
            Array of shape (n,3) with columns of X, Y, Z ECI position data.

        Usage
        -----
        Must have ECEFdata and ERAdata attributes initiated or posData and/or
        timeData inputs.
        """
        if type(self.ERAdata) == type(None):
            self.getERA(**kwargs)
        if type(self.ECEFdata) == type(None):
            ECEFdata = kwargs["posData"]
            self.positionData(ECEFdata, "ECEF")

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

    def getECEF(self, **kwargs):
        """
        Instantiate ECEFdata attribute.

        Computes ECEFdata attribute from ECIdata and ERAdata attribites.

        Parameters
        ----------
        None

        **kwargs
        --------
        posData : ndarray
            Array of shape (n,3) with columns of X, Y, Z position data assumed
            ECI.
        timeData : ndarray
            Array of shape (n,). Contains datetime objects in UTC.

        Returns
        -------
        self.ECEFdata : ndarray
            Array of shape (n,3) with columns of X, Y, Z ECEF position data.

        Usage
        -----
        Must have ECIdata and ERAdata attributes initiated or posData and/or
        timeData inputs.
        """
        if type(self.ERAdata) == type(None):
            self.getERA(**kwargs)
        if type(self.ECIdata) == type(None):
            ECIdata = kwargs["posData"]
            self.positionData(ECIdata, "ECI")

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

    def getAltAz(self, obsCoor, radius, **kwargs):
        """
        Instantiate horizontal attribute.

        This method takes in an observers latitude and longitude in degrees as
        well as the surface radius in meters to instantiate the horisontal
        attribute.

        Parameters
        ----------
        obsCoor : tuple of floats
            Specifies observer position in degrees, (latitude, longitude).
        radius : int or float
            Fixed radius of Earth.

        **kwargs
        --------
        posData : ndarray
            Array of shape (n,3) with columns of X, Y, Z position data assumed
            ECEF.
        timeData : ndarray
            Array of shape (n,). Contains datetime objects in UTC.

        Returns
        -------
        self.horizontal : ndarray
            Array of shape (n,2) with columns of alt, az position data.

        Usage
        -----
        Must have ECEFdata attribute initiated or the posData and timeData
        inputs.
        """
        def getAngle(vecOne, vecTwo):
            """
            Calculate degree angle bewteen two vectors.

            Takes two multidimensional cartesian coordinate vectors and returns
            degree angle between the two arrays element-wise.

            Parameters
            ----------
            vecOne : ndarray
                Array of shape (n,3) representing n XYZ vectors.
            vecTwo : ndarray
                Array of shape (n,3) representing n XYZ vectors.

            Returns
            -------
            ang : float
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

        # Convert observer position into spherical then cartesian.
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
        altAz = np.zeros((self.length, 2))
        altAz[:, 0] = 90 - getAngle(xyzLOS, xyzObs)

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

        altAz[posInd, 1] = getAngle(tangentVec[posInd], projLOS[posInd])
        altAz[negInd, 1] = 360 - getAngle(tangentVec[negInd], projLOS[negInd])

        self.horizontal = altAz

        return self.horizontal

    def getNdrAng(self, obsCoor, radius, **kwargs):
        """
        Return Nadir-LOS angle array.

        This method takes in an observers latitude and longitude in degrees as
        well as the surface radius in meters to return the angle made between
        the satellites nadir and the line of light from the observer to the
        ground station.

        Parameters
        ----------
        obsCoor : tuple of floats
            Specifies observer position in degrees, (latitude, longitude).
        radius : int or float
            Fixed radius of Earth.

        **kwargs
        --------
        posData : ndarray
            Array of shape (n,3) with columns of X, Y, Z position data assumed
            ECI.
        timeData : ndarray
            Array of shape (n,). Contains datetime objects in UTC.

        Returns
        -------
        self.nadirAng : ndarray
            Array of shape (n,) with nadir-LOS angle data.

        Usage
        -----
        Must have ECEFdata attribute initiated or the posData and timeData
        inputs.
        """
        if type(self.ECEFdata) == type(None):
            self.getECEF(**kwargs)

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

        self.nadirAng = ang

        return self.nadirAng

    def saveData(self, fileName, delimiter):
        """
        Save class data in local directory.

        Parameters
        ----------
        fileName : str
            File name of output file. Can be a .txt or .csv file.
        delimiter : str
            Feild delimiter for the output file. String of length 1.

        Returns
        -------
        None

        Usage
        -----
        It is recommended to use a tab delimiter for .txt files and comma
        delimiters for .csv files.
        """
        data = {}

        if type(self.times) != type(None):
            data['Time (UTC)'] = pd.Series(self.times)

        if type(self.ERAdata) != type(None):
            data['ERA (Deg)'] = pd.Series(self.ERAdata)

        if type(self.ECIdata) != type(None):
            data['ECI.X'] = pd.Series(self.ECIdata[:, 0])
            data['ECI.Y'] = pd.Series(self.ECIdata[:, 1])
            data['ECI.Z'] = pd.Series(self.ECIdata[:, 2])

        if type(self.ECEFdata) != type(None):
            data['ECEF.X'] = pd.Series(self.ECEFdata[:, 0])
            data['ECEF.Y'] = pd.Series(self.ECEFdata[:, 1])
            data['ECEF.Z'] = pd.Series(self.ECEFdata[:, 2])

        if type(self.horizontal) != type(None):
            data['Alt (Deg)'] = pd.Series(self.horizontal[:, 0])
            data['Az (Deg)'] = pd.Series(self.horizontal[:, 1])

        if type(self.nadirAng) != type(None):
            data['NadirAng (Deg)'] = pd.Series(self.nadirAng)

        df = pd.DataFrame(data)
        df.to_csv(fileName, sep=delimiter)
