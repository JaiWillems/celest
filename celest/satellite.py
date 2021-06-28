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
from scipy.interpolate import interp1d
from datetime import datetime
import julian
from typing import Union, Tuple
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
    GEOdata : np.array
        Array of shape (n, 3) containing Geographical coordinates where the
        columns are lattitude, longitude, radius data.
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
        self.GEOdata = None
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
        >>> UTCTimeData = np.array(["2020-06-01 12:00:00.0340", ...,
        ...                         "2020-06-01 12:01:00.0340"])
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

    def position_data(self, posData: np.array, type: str, factor: int=0, timeData: np.array=None, jul: bool=False, julOffset: float=0) -> None:
        """Instantiate orbital position data.

        Parameters
        ----------
        posData : np.array
            If `type="ECI"` or `type="ECEF", input will be an array of shape
            (n,3) with columns of X, Y, Z position data. If `type="GEO"`, input
            can be of shape (n, 2) with columns of latitude, longitude data or
            of shape (n, 3) with the last column being a radius component. Note
            that angular data must be passed in as decimal degrees.
        type : {"GEO", "ECI", "ECEF"}
            Specifies posData as either "GEO", "ECI" or "ECEF".
        factor : int, optional
            Interpolate the inputted data by the factor `factor`.
        timeData : np.array
            Array of shape (n,) containing time data.
        jul : bool, optional
            Indicates the input times are Julian dates if True or UTC datetime
            strings if False.
        julOffset : float, optional
            Offset to be added to the inputed Julian dates.

        Notes
        -----
        For more coincise code, the posData can be passed directly into all
        conversion methods.

        If `type="GEO" and posData is of shape (n, 2), the radius componenet
        will be filled with the Earths radius at the corresponding lattitude,
        longitude position.

        To interpolate the position data upon input, the `timeData` attribute
        must be instantiated or passed in as a parameter.

        Examples
        --------
        >>> finch = Satellite()
        >>> ECIvec = np.array([[-4.46e+03, -5.22e+03, 1.75e-04], ...,
        ...                    [2.73e+03, 2.08e+03, -6.02e+03]])
        >>> finch.position_data(posData=ECIvec, type="ECI")
        """
        if type(self.timeData) == type(None) and type(timeData) != type(None):
            self.time_data(timeData=timeData, jul=jul, julOffset=julOffset)
    
        if type == "GEO" and posData.shape[1] == 2:
            radius = self._WGS84_radius(posData[:, 0])
            posData = np.concatenate((posData, radius), axis=1)

        if factor > 0:
            times = self.timeData
            timesP = np.linspace(min(times), max(times), factor * times.shape[0])
            self.time_data(timeData=timesP)
            newData = np.zeros((factor * self.length, 3))
            for i in range(3):
                func = interp1d(times, posData[:, i])
                newData[:, i] = func(self.timeData)
        else:
            newData = timeData

        if type == "GEO":
            self.GEOdata = newData
        elif type == "ECI":
            self.ECIdata = newData
        elif type == "ECEF":
            self.ECEFdata = newData

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
        >>> UTCTimeData = np.array(["2020-06-01 12:00:00.0340", ...,
        ...                         "2020-06-01 12:01:00.0340"])
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
    
    def _ECEF_to_GEO(self, ECEFpos: np.array) -> np.array:
        """Convert from ECEF to geographical coordinates.

        Parameters
        ----------
        ECEFpos : np.array
            Array of shape (n, 3) with rows containing XYZ ECEF position data.
        
        Returns
        -------
        np.array
            Array of shape (n, 3) with columns of lattitude, longitude, radius
            data in decimal degrees and km.
        """
        # Cartesian coordinates.
        x = ECEFpos[:, 0]
        y = ECEFpos[:, 1]
        z = ECEFpos[:, 2]

        # Convert to spherical coordinates.
        radius = np.sqrt(x**2 + y**2 + z**2)
        theta = np.degrees(np.arccos(np.divide(z, radius)))
        phi = np.degrees(np.arctan(np.divide(y, x)))
        phi[np.where(x < 0)[0]] = phi[np.where(x < 0)[0]] - 180

        # Convert to geographical coordinates.
        lat = 90 - theta
        lon = phi
        lon[np.where(lon < 180)[0]] = phi[np.where(lon < 180)[0]] + 360
        lon[np.where(lon >= 180)[0]] = phi[np.where(lon >= 180)[0]] - 360

        # Fomulate output array.
        lat = lat.reshape(-1, 1)
        lon = lon.reshape(-1, 1)
        radius = radius.reshape(-1, 1)
        geo = np.concatenate((lat, lon, radius), axis=1)

        return geo

    def _GEO_to_ECEF(self, geoPos: np.array) -> np.array:
        """Convert from geographical to ECEF coordinates.

        Parameters
        ----------
        geoPos : np.array
            Array of shape (n, 3) containing Geodetic coordinates for a location
            with columns of lattitude, longitude, radius given in decimal
            degrees and km.
        radius : np.array, optional
            Radius of the position from the Earth's surface.
        
        Returns
        -------
        np.array
            Array of shape (n, 3) containing the XYZ ECEF position data.
        """
        radius = geoPos[:, 2]
        
        # Convert geographical to spherical.
        theta = geoPos[:, 1]
        negLon = np.where(theta < 0)
        posLon = np.where(theta >= 0)
        theta[negLon] = np.radians(theta[negLon] + 360)
        theta[posLon] = np.radians(theta[posLon])
        phi = np.radians(90 - geoPos[:, 0])

        # Convert spherical to cartesian.
        x = radius * np.cos(theta) * np.sin(phi).reshape(-1, 1)
        y = radius * np.sin(theta) * np.sin(phi).reshape(-1, 1)
        z = radius * np.cos(phi).reshape(-1, 1)
        ECEFpos = np.concatenate((x, y, z), axis=1)

        return ECEFpos

    def GEO(self, posData: np.array=None) -> np.array:
        """Instantiate the GEOdata attribute.

        Parameters
        ----------
        posData : np.array
            Array of shape (n, 3) with coumns of X, Y, Z ECEF position data.
        
        Returns
        -------
        np.array
            Array of shape (n, 3) with rows of lattitude, longitude, radius
            position data.
        """
        if type(self.ECEFdata) == type(None):
            self.position_data(posData=posData, type="ECEF")

        self.GEOdata = self._ECEF_to_GEO(ECEFpos=self.ECEFdata)

        return self.GEOdata

    def _ECI_and_ECEF(self, posData: np.array, type: str) -> np.array:
        """Convert between ECI and ECEF positions.

        Parameters
        ----------
        posData : np.array
            Array of shape (n, 3) representing the input data as XYZ cartesian
            data.
        type : {"ECI", "ECEF"}
            The type of inputted data.
        
        Returns
        -------
        np.array
            Array of shape (n, 3) representing the output data as XYZ cartesian
            data.
        """
        if type == "ECI":
            theta = -self.ERAdata
        else:
            theta = self.ERAdata
        
        # Construct rotational matrix.
        outVec = np.zeros((self.length, 3))
        A11 = np.cos(theta)
        A12 = -np.sin(theta)
        A21 = np.sin(theta)
        A22 = np.cos(theta)

        # Rotate position data around z-axis by ERA.
        outVec[:, 0] = np.add(np.multiply(A11, self.posData[:, 0]),
                              np.multiply(A12, self.posData[:, 1]))
        outVec[:, 1] = np.add(np.multiply(A21, self.posData[:, 0]),
                              np.multiply(A22, self.posData[:, 1]))
        outVec[:, 2] = self.posData[:, 2]

        return outVec

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
        ECEF : Initiate ECEFdata attribute.

        Notes
        -----
        The Satellite class instance must have ECEFdata and ERAdata attributes
        initiated or posData and timeData inputs passed in.

        Examples
        --------
        >>> UTCTimeData = np.array(["2020-06-01 12:00:00.0340", ...,
        ...                         "2020-06-01 12:01:00.0340"])
        >>> ECEFvec = np.array([[-4.46e+03, -5.22e+03, 1.75e-04], ...,
        ...                     [2.73e+03, 2.08e+03, -6.02e+03]])
        >>> finch = Satellite()
        >>> ECIvec = finch.ECI(posData=ECEFvec, timeData=UTCTimeData)
        """
        if type(self.ERAdata) == type(None):
            self.ERA(timeData=timeData, **kwargs)
        if type(self.ECEFdata) == type(None):
            self.position_data(posData=posData, type="ECEF")

        self.ECIdata = self._ECI_and_ECEF(posData=self.ECEFdata, type="ECEF")

        return self.ECIdata

    def ECEF(self, posData: np.array=None, timeData: np.array=None, type="ECI", **kwargs:
             Union[bool, float]) -> np.array:
        """Instantiate ECEFdata attribute.

        Parameters
        ----------
        posData : np.array, optional
            Array of shape (n,3) with columns of X, Y, Z ECI position data.
        timeData : np.array, optional
            Array of shape (n,) containing time data.
        type : {"ECI", "GEO"}, optional
            Specifies posData as either "ECI" or "GEO".
        **kwargs : dict, optional
            Extra arguments to `getECEF`: refer to getECEF documentation for a
            list of all possible arguments.

        Returns
        -------
        np.array
            Array of shape (n,3) with columns of X, Y, Z ECEF position data.

        See Also
        --------
        ECI : Initiate ECIdata attribute.

        Notes
        -----
        The Satellite class instance must have ECIdata and ERAdata attributes
        initiated or posData and timeData inputs passed in.

        Examples
        --------
        >>> UTCTimeData = np.array(["2020-06-01 12:00:00.0340", ...,
        ...                         "2020-06-01 12:01:00.0340"])
        >>> ECIvec = np.array([[-4.46e+03, -5.22e+03, 1.75e-04], ...,
        ...                    [2.73e+03, 2.08e+03 -6.02e+03]])
        >>> finch = Satellite()
        >>> ECEFvec = finch.ECEF(posData=ECIvec, timeData=UTCTimeData)
        """
        if type(self.ECIdata) == type(None) and type(self.GEOdata) == type(None):
            if type == "ECI":
                self.position_data(posData=posData, type="ECI")

            elif type == "GEO":
                self.position_data(posData=posData, type="GEO")

        if type(self.GEOdata) != type(None):
            self.ECEFdata = self._GEO_to_ECEF(posData=self.GEOdata)

        elif type(self.ECIdata) != type(None):
            if type(self.ERAdata) == type(None):
                self.ERA(timeData=timeData, **kwargs)
            self.ECEFdata = self._ECI_and_ECEF(posData=self.ECIdata, type="ECI")

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
        dividend = np.einsum("ij, ij->i", vecOne, vecTwo)
        divisor = np.multiply(np.linalg.norm(
            vecOne, axis=1), np.linalg.norm(vecTwo, axis=1))
        arg = np.divide(dividend, divisor)
        ang = np.degrees(np.arccos(arg))

        return ang

    def horizontal(self, groundPos: GroundPosition, posData: np.array=None,
                   timeData: np.array=None, **kwargs:
                   Union[bool, float]) -> Tuple[np.array, np.array]:
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
        Tuple
            The altitude/azimuth data is returned in a tuple where both the
            altitude and azimuth are arrays of shape (n,).

        Notes
        -----
        The Satellite class instance must have the ECEFdata attribute initiated
        or have the posData and timeData inputs passed in.

        Examples
        --------
        >>> UTCTimeData = np.array(["2020-06-01 12:00:00.0340", ...,
        ...                         "2020-06-01 12:01:00.0340"])
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
        coeff = np.einsum("ij, ij->i", xyzLOS, xyzObs)/radius**2
        normProj = (xyzObs.T * coeff).T
        projLOS = np.subtract(xyzLOS, normProj)

        # Determing azimuth.
        vecOne = np.cross(tangentVec, xyzObs)
        normOne = 1/np.linalg.norm(vecOne, axis=1).reshape((self.length, 1))
        vecOneUnit = normOne*vecOne
        vecTwo = (vecOneUnit.T * np.einsum("ij, ij->i", projLOS, vecOneUnit)).T
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
        >>> UTCTimeData = np.array(["2020-06-01 12:00:00.0340", ...,
        ...                         "2020-06-01 12:01:00.0340"])
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

    def _WGS84_radius(self, lattitude: np.array) -> np.array:
        """Calculate the Earth's radius using WGS84.

        Parameters
        ----------
        latitude : np.array
            Array of shape (n,) representing the lattitude of a ground
            location.

        Returns
        -------
        np.array
            Earth's radius at each row in `lattitude` using WGS84.

        Notes
        -----
        The Earth can be modeled as an ellipsoid given by the following
        equation:

        .. math:: r = \sqrt{\\frac{(6378.14)^2(6356.75)^2}{(6378.14)^2\sin{\phi}^2+(6356.75)^2\cos{\phi}^2}}

        where :math:`\phi` is the observers lattitude.
        """
        # Get lattidue parameter.
        phi = np.radians(lattitude)

        # Define WGS84 Parameters.
        semiMajor = 6378.137**2
        semiMinor = 6356.752314245**2

        numerator = semiMajor * semiMinor
        denominator = semiMajor * np.sin(phi)**2 + semiMinor * np.cos(phi)**2
        radius = np.sqrt(numerator / denominator)

        return radius

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
    
    def _times_broadening(self, times, indices, factor, dt):
        """Broaden the times data.

        Parameters
        ----------
        times : np.array
            Array of times to broaden.
        indices : np.array
            Array of indices defining the `times` region to broaden.
        factor : int
            The factor increase in the number of steps in the broadening region.
        dt : int, optional
            The number of data points adjacent to the broadening region to include
            in the times broadening process.
        
        Returns
        -------
        np.array
            New `times` array with the region defined by `indices` broadened.
        np.array
            An array of the broadened times.
        np.array
            Updated indices of the broadened region in the `times` array.
        """
        # Reform the viable region.
        indices = np.append(indices, np.arange(indices[-1] + 1, indices[-1] + 1 + dt, 1))
        indices = np.insert(indices, 0, np.arange(indices[0] - dt, indices[0], 1))

        minI = indices[0]
        maxI = indices[-1]

        # Determine new times and insert into current time data.
        timeInsert = np.linspace(times[minI], times[maxI], int(factor * (maxI - minI + 1)))
        times = np.delete(times, indices)
        times = np.insert(times, minI, timeInsert)

        return times, timeInsert, indices
    
    def interp_ECI(self, factor, dt=0, indices=None):
        """Interpolate ECI data.
        
        If `indices=None`, this method will interpolate all ECI data, otherwise
        it will interpolate the regions defined by `indices`.
        
        Parameters
        ----------
        factor : int
            The factor increase in the number of steps in interpolated regions.
        dt : int, optional
            The number of data points adjacent to a valid encounter region to
            interpolate within. Only used if `indices!=None`.
        indices : np.array, optional
            An array of arrays of indices defining encounter regions.
            
        returns
        -------
        np.array
            An array of shape (n, 3) containing X, Y, Z interpolated ECI data.
        np.array
            An array of shape (n,) containing interpolated julian times.
        """
        # Get position information.
        times = self.times
        ECIdata = self.ECIdata

        # Set up ECEF interpolation functions.
        xInterp = interp1d(times, ECIdata[:, 0], kind="cubic")
        yInterp = interp1d(times, ECIdata[:, 1], kind="cubic")
        zInterp = interp1d(times, ECIdata[:, 2], kind="cubic")

        if type(indices) != type(None):
            for i in np.flip(indices):

                times, time, i = self._times_broadening(times, i, factor, dt)

                # Interpolate, insert, and reform ECEF data.
                x = np.delete(ECIdata[:, 0], i)
                y = np.delete(ECIdata[:, 1], i)
                z = np.delete(ECIdata[:, 2], i)
                xP = np.insert(x, i[0], xInterp(time)).reshape((-1, 1))
                yP = np.insert(y, i[0], yInterp(time)).reshape((-1, 1))
                zP = np.insert(z, i[0], zInterp(time)).reshape((-1, 1))
                ECIdata = np.concatenate((xP, yP, zP), axis=1)

        else:
            timesP = np.linspace(min(times), max(times), factor * times.size)
            xP, yP, zP = xInterp(timesP), yInterp(timesP), zInterp(timesP)
            ECIdata = np.concatenate((xP, yP, zP), axis=1)

        return ECIdata, times
    
    def interp_ECEF(self, factor, dt=0, indices=None):
        """Interpolate ECEF data.
        
        If `indices=None`, this method will interpolate all ECEF data,
        otherwise it will interpolate the regions defined by `indices`.
        
        Parameters
        ----------
        factor : int
            The factor increase in the number of steps in interpolated regions.
        dt : int, optional
            The number of data points adjacent to a valid encounter region to
            interpolate within. Only used if `indices!=None`.
        indices : np.array, optional
            An array of arrays of indices defining encounter regions.
            
        returns
        -------
        np.array
            An array of shape (n, 3) containing X, Y, Z interpolated ECEF data.
        np.array
            An array of shape (n,) containing interpolated julian times.
        """
        # Get position information.
        times = self.times
        ECEFdata = self.ECEFdata

        # Set up ECEF interpolation functions.
        xInterp = interp1d(times, ECEFdata[:, 0], kind="cubic")
        yInterp = interp1d(times, ECEFdata[:, 1], kind="cubic")
        zInterp = interp1d(times, ECEFdata[:, 2], kind="cubic")

        if type(indices) != type(None):
            for i in np.flip(indices):

                times, time, i = self._times_broadening(times, i, factor, dt)

                # Interpolate, insert, and reform ECEF data.
                x = np.delete(ECEFdata[:, 0], i)
                y = np.delete(ECEFdata[:, 1], i)
                z = np.delete(ECEFdata[:, 2], i)
                xP = np.insert(x, i[0], xInterp(time)).reshape((-1, 1))
                yP = np.insert(y, i[0], yInterp(time)).reshape((-1, 1))
                zP = np.insert(z, i[0], zInterp(time)).reshape((-1, 1))
                ECEFdata = np.concatenate((xP, yP, zP), axis=1)

        else:
            timesP = np.linspace(min(times), max(times), factor * times.size)
            xP, yP, zP = xInterp(timesP), yInterp(timesP), zInterp(timesP)
            ECEFdata = np.concatenate((xP, yP, zP), axis=1)

        return ECEFdata, times
    
    def interp_alt(self, groundPos, factor, dt=0, indices=None):
        """Interpolate altitude angle data.
        
        If `indices=None`, this method will interpolate all altitude data,
        otherwise it will interpolate the regions defined by `indices`.
        
        Parameters
        ----------
        groundPos : GroundPosition
            The ground location containing the altitude to interpolate.
        factor : int
            The factor increase in the number of steps in interpolated regions.
        dt : int, optional
            The number of data points adjacent to a valid encounter region to
            interpolate within. Only used if `indices!=None`.
        indices : np.array, optional
            An array of arrays of indices defining encounter regions.
            
        returns
        -------
        np.array
            An array of shape (n,) containing interpolated altitude data.
        np.array
            An array of shape (n,) containing interpolated julian times.
        """
        # Get position information.
        times = self.times
        alt = groundPos.alt

        # Set up altitude interpolation function.
        altInterp = interp1d(times, alt, kind="cubic")

        if type(indices) != type(None):
            for i in np.flip(indices):
                times, time, i = self._times_broadening(times, i, factor, dt)
                temp = np.delete(alt, i)
                altP = np.insert(temp, i[0], altInterp(time)).reshape((-1, 1))
        else:
            timesP = np.linspace(min(times), max(times), factor * times.size)
            altP = altInterp(timesP)

        return altP, times

    def interp_az(self, groundPos, factor, dt=0, indices=None):
        """Interpolate azimuth data.
        
        If `indices=None`, this method will interpolate all azimuth data,
        otherwise it will interpolate the regions defined by `indices`.
        
        Parameters
        ----------
        groundPos : GroundPosition
            The ground location containing the azimuth to interpolate.
        factor : int
            The factor increase in the number of steps in interpolated regions.
        dt : int, optional
            The number of data points adjacent to a valid encounter region to
            interpolate within. Only used if `indices!=None`.
        indices : np.array, optional
            An array of arrays of indices defining encounter regions. 
           
        returns
        -------
        np.array
            An array of shape (n,) containing interpolated azimuth data.
        np.array
            An array of shape (n,) containing interpolated julian times.
        """
        # Get position information.
        times = self.times
        az = groundPos.az

        # Set up azimuth interpolation function.
        azInterp = interp1d(times, az, kind="cubic")

        if type(indices) != type(None):
            for i in np.flip(indices):
                times, time, i = self._times_broadening(times, i, factor, dt)
                temp = np.delete(az, i)
                azP = np.insert(temp, i[0], azInterp(time)).reshape((-1, 1))
        else:
            timesP = np.linspace(min(times), max(times), factor * times.size)
            azP = azInterp(timesP)

        return azP, times

    def interp_nadir(self, groundPos, factor, dt=0, indices=None):
        """Interpolate nadir angle data.
        
        If `indices=None`, this method will interpolate all nadir angle data,
        otherwise it will interpolate the regions defined by `indices`.
        
        Parameters
        ----------
        groundPos : GroundPosition
            The ground location containing the nadir angles to interpolate.
        factor : int
            The factor increase in the number of steps in interpolated regions.
        dt : int, optional
            The number of data points adjacent to a valid encounter region to
            interpolate within. Only used if `indices!=None`.
        indices : np.array, optional
            An array of arrays of indices defining encounter regions.
            
        returns
        -------
        np.array
            An array of shape (n,) containing interpolated nadir angle data.
        np.array
            An array of shape (n,) containing interpolated julian times.
        """
        # Get position information.
        times = self.times
        nadir = groundPos.nadirAng

        # Set up nadir interpolation function.
        nadirInterp = interp1d(times, nadir, kind="cubic")

        if type(indices) != type(None):
            for i in np.flip(indices):
                times, time, i = self._times_broadening(times, i, factor, dt)
                temp = np.delete(nadir, i)
                nadirP = np.insert(temp, i[0], nadirInterp(time)).reshape((-1, 1))
        else:
            timesP = np.linspace(min(times), max(times), factor * times.size)
            nadirP = nadirInterp(timesP)

        return nadirP, times
    
    def interp_distance(self, groundPos, factor, dt=0, indices=None):
        """Interpolate distance data.
        
        If `indices=None`, this method will interpolate all distance data,
        otherwise it will interpolate the regions defined by `indices`.
        
        Parameters
        ----------
        groundPos : GroundPosition
            The ground location containing the distances to interpolate.
        factor : int
            The factor increase in the number of steps in interpolated regions.
        dt : int, optional
            The number of data points adjacent to a valid encounter region to
            interpolate within. Only used if `indices!=None`.
        indices : np.array, optional
            An array of arrays of indices defining encounter regions.
            
        returns
        -------
        np.array
            An array of shape (n,) containing interpolated distance data.
        np.array
            An array of shape (n,) containing interpolated julian times.
        """
        # Get position information.
        times = self.times
        distance = groundPos.distance

        # Set up distance interpolation function.
        distanceInterp = interp1d(times, distance, kind="cubic")

        if type(indices) != type(None):
            for i in np.flip(indices):
                times, time, i = self._times_broadening(times, i, factor, dt)
                temp = np.delete(distance, i)
                distanceP = np.insert(temp, i[0], distanceInterp(time)).reshape((-1, 1))
        else:
            timesP = np.linspace(min(times), max(times), factor * times.size)
            distanceP = distanceInterp(timesP)

        return distanceP, times
    
    def interp_elev(self, factor, dt=0, indices=None):
        """Interpolate satellite elevation data.
        
        If `indices=None`, this method will interpolate all elevation data,
        otherwise it will interpolate the regions defined by `indices`.
        
        Parameters
        ----------
        factor : int
            The factor increase in the number of steps in interpolated regions.
        dt : int, optional
            The number of data points adjacent to a valid encounter region to
            interpolate within. Only used if `indices!=None`.
        indices : np.array, optional
            An array of arrays of indices defining encounter regions.
            
        returns
        -------
        np.array
            An array of shape (n,) containing interpolated elevation data.
        np.array
            An array of shape (n,) containing interpolated julian times.
        """
        # Get position information.
        times = self.times
        altitude = self.altitude

        # Set up elevation interpolation function.
        altitudeInterp = interp1d(times, altitude, kind="cubic")

        if type(indices) != type(None):
            for i in np.flip(indices):
                times, time, i = self._times_broadening(times, i, factor, dt)
                temp = np.delete(altitude, i)
                altitudeP = np.insert(temp, i[0], altitudeInterp(time)).reshape((-1, 1))
        else:
            timesP = np.linspace(min(times), max(times), factor * times.size)
            altitudeP = altitudeInterp(timesP)

        return altitudeP, times
    
    def _GEO_to_sexagesimal(self) -> np.array:
        """Present geographical coordinates in sexagesimal format.

        This method takes geographical coordinates in decimal degrees and km
        and returns strings for each coordinate following the ISO 6709:2008
        sexagesimal format.

        Returns
        -------
        np.array
            Array of shape (n,) containing sexagesimal angles as strings.
        """
        lat = self.GEOdata[:, 0]
        long = self.GEOdata[:, 1]
        alt = self.GEOdata[:, 2]
        length = self.GEOdata.shape[0]
        sexagesimalAng = np.empty((length, 3), dtype="<U32")
        for i in range(length):
            num = lat[i]
            sign = "N" if num >= 0 else "S"
            degrees = str(int(abs(num) - num % 1)).zfill(2)
            minutes = str(int(60 * (num % 1) - (60 * (num % 1)) % 1)).zfill(2)
            seconds = "{:.3f}".format(60 * ((60 * (num % 1)) % 1)).zfill(5)
            degree_symbol = u"\u00B0"
            minute_symbol = u"\u2032"
            second_symbol = u"\u2033"
            sexagesimalAng[i][0] = f"{degrees}{degree_symbol}{minutes}{minute_symbol}{seconds}{second_symbol}{sign}"

            num = long[i]
            sign = "E" if num >= 0 else "W"
            degrees = str(int(abs(num) - num % 1)).zfill(2)
            minutes = str(int(60 * (num % 1) - (60 * (num % 1)) % 1)).zfill(2)
            seconds = "{:.3f}".format(60 * ((60 * (num % 1)) % 1)).zfill(5)
            degree_symbol = u"\u00B0"
            minute_symbol = u"\u2032"
            second_symbol = u"\u2033"
            sexagesimalAng[i][1] = f"{degrees}{degree_symbol}{minutes}{minute_symbol}{seconds}{second_symbol}{sign}"

            sexagesimalAng[i][2] = f"{alt[i]:.3f}km"

        return sexagesimalAng
    
    def _deg_to_sexagesimal(self, angles: np.array) -> np.array:
        """Convert decimal angles into sexagesimal angles.

        Parameters
        ----------
        angles : np.array
            Array of shape (n,) containing float angles.

        Returns
        -------
        np.array
            Array of shape (n,) containing sexagesimal angles as strings.
        """
        length = angles.shape[0]
        sexagesimalAng = np.empty((length,), dtype="<U32")
        for i in range(length):
            num = angles[i]
            sign = "+" if num >= 0 else "-"
            degrees = str(int(abs(num) - num % 1)).zfill(2)
            minutes = str(int(60 * (num % 1) - (60 * (num % 1)) % 1)).zfill(2)
            seconds = "{:.3f}".format(60 * ((60 * (num % 1)) % 1)).zfill(5)
            degree_symbol = u"\u00B0"
            minute_symbol = u"\u2032"
            second_symbol = u"\u2033"
            sexagesimalAng[i] = f"{sign}{degrees}{degree_symbol}{minutes}{minute_symbol}{seconds}{second_symbol}"
        return sexagesimalAng

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
            data["Time (UTC)"] = pd.Series(self.times)

        if type(self.ERAdata) != type(None):
            data["ERA (rad)"] = pd.Series(self.ERAdata)

        if type(self.ECIdata) != type(None):
            data["ECI.X"] = pd.Series(self.ECIdata[:, 0])
            data["ECI.Y"] = pd.Series(self.ECIdata[:, 1])
            data["ECI.Z"] = pd.Series(self.ECIdata[:, 2])

        if type(self.ECEFdata) != type(None):
            data["ECEF.X"] = pd.Series(self.ECEFdata[:, 0])
            data["ECEF.Y"] = pd.Series(self.ECEFdata[:, 1])
            data["ECEF.Z"] = pd.Series(self.ECEFdata[:, 2])
        
        if type(self.elevation) != type(None):
            data["Elevation"] = pd.Series(self.elevation)

        for key in self.gs:
            if type(self.gs[key].alt) != type(None):
                data[f"{key}.alt"] = self.gs[key].alt
                data[f"{key}.az"] = self.gs[key].az
            if type(self.gs[key].nadirAng) != type(None):
                data[f"{key}.NdrAng"] = self.gs[key].nadirAng
            if type(self.gs[key].distance) != type(None):
                data[f"{key}.distance"] = self.gs[key].distance

        df = pd.DataFrame(data)
        df.to_csv(fileName, sep=delimiter)
