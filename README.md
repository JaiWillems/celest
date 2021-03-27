# Celest

The Celest library was created to provide a simple interface to satellite orbital position representations at an unprecedented speed. This library is designed for small satellite applications who desire efficiency over precision in orbital calculations.

## Installation

```terminal
pip install Celest
```

## Class Overview
Celest is broken into two classes, a Satellite class and a GroundPosition class. This section details the specifics for each class.
### Satellite Class

The Satellite class represents the orbital object to store and compute orbital position representations. The class stores its data in six instance variables:

1. **self.times**: Orbital time dependencies in an (n,) shaped ndarray.
1. **self.ERAdata**: Earth rotation angles based on time dependencies in an (n,) shaped ndarray.
1. **self.ECIdata**: Earth centered inertial position data in an (n,3) shaped ndarray where the columns are the x, y, z position data.
1. **self.ECEFdata**: Earth centered earth fixed position data in an (n,3) shaped ndarray where the columns are the x, y, z position data.
1. **self. gs** : Dictionary with the GroundPosition *name* attribute as the key and the GroundPosition object as the corresponding value.
1. **self.length**: Length of data attributes represented as an int

These instance variables are interfaced by the user through eight methods:

1. **timeData**: Instantiates *times* attribute with orbital time dependency.
1. **positionData**: Instantiates *ECIdata* or *ECEFdata* attribute with orbital position data.
1. **getERA**: Instantiates *ERAdata* attribute.
1. **getECI**: Instantiates *ECIdata* attribute.
1. **getECEF**: Instantiates *ECEFdata* attribute.
1. **getAltAz**: Instantiates GroundPosition's *alt* and *az* attributes.
1. **getNdrAng**: Instantiates GroundPosition's *nadirAng* attribute.
1. **saveData**: Save class data in local directory.

The guiding principle behind the design of the Satellite class is to be user simple and efficient. The restrictions on how to use these methods are outlined in the following subsections along with detailed explanations and example usage of each of the methods.

#### .timeData(timeData)
**Description**  
The `timeData` method instantiates the *times* attribute with the orbital time dependencies (the times that each position observation was taken). The length of the input data initializes the *length* attribute.  
**Parameters**  
`timeData`: ndarray of shape (n,) containing `datetime.datetime` objects in UTC.  
**Returns**  
None  
**Usage**  
None  
**Example**  
```python
from celest import Satellite
import numpy as np

finch = Satellite()
UTCTimeData = np.array(['2020-06-01 12:00:00.0340', ..., '2020-06-01 12:01:00.0340'])
finch.timeData(timeData=UTCtimeData)
```


#### .positionData(posData, type)
**Description**  
The `posData` method interfaces with the *ECIdata* or *ECEFdata* attributes to initialize position data of the satellite. The user must specify the data as either ECI or ECEF data through the `type` parameter.  
**Parameters**  
`posData`: ndarray of shape (n,3) with columns of xyz position data.  
`type`: Specifies the type of position data as either "ECI" or "ECEF".  
**Returns**  
None  
**Usage**  
The length of `posData` must be the same as the *length* attribute.  
**Example**  
```python
from celest import Satellite
import numpy as np

finch = Satellite()
ECIvec = np.array([[-4.46e+03 -5.22e+03  1.75e-04], ..., [ 2.73e+03  2.08e+03 -6.02e+03]])
finch.positionData(posData=ECIvec, type="ECI")
```

#### .getERA(**kwargs)
**Description**  
The `getERA` method uses the *times* variable to calculate the earth rotation angles in radians. This is the angle between the <img src="https://render.githubusercontent.com/render/math?math=x_{ECI}"> and <img src="https://render.githubusercontent.com/render/math?math=x_{ECEF}"> coordinate axes. The array of these angles is stored in *ERAdata*.  
**Parameters**  
None  
**\*\*kwargs**  
`timeData`: ndarray of shape (n,) containing `datetime.datetime` objects in UTC.  
**Returns**  
*ERAdata* : ndarray of shape (n,) containing radian earth rotation angles.  
**Usage**  
The *times* attribute must be initiated to use `getERA`. If one wishes to combine the `timeData` and `getERA` methods then the time dependency must be passed in as a **kwarg.  
**Example**  
```python
from celest import Satellite
import numpy as np

finch = Satellite()
UTCTimeData = np.array(['2020-06-01 12:00:00.0340', ..., '2020-06-01 12:01:00.0340'])
finch.timeData(timeData=UTCtimeData)
ERAangles = finch.getERA()
```

The commands can be simplified as in the following,  

```python
from celest import Satellite
import numpy as np

finch = Satellite()
UTCTimeData = np.array(['2020-06-01 12:00:00.0340', ..., '2020-06-01 12:01:00.0340'])
ERAangles = finch.getERA(timeData=UTCTimeData)
```

#### .getECI(**kwargs)
**Description**  
The `getECI` method uses the earth rotation angles to convert ECEF data to ECI data while instantiating *ECIdata*.  
**Parameters**  
None  
**\*\*kwargs**  
`posData`: ndarray array of shape (n,3) with columns of X, Y, Z position data assumed ECEF.  
`timeData`: ndarray of shape (n,) containing `datetime.datetime` objects in UTC.  
**Returns**  
*ECIdata* : ndarray of shape (n,3) with columns of X, Y, Z ECI position data.  
**Usage**  
The *ECEFdata* and *times* attributes must be instantiated or passed in as **kwargs.  
**Example**  
```python
from celest import Satellite
import numpy as np

UTCTimeData = np.array(['2020-06-01 12:00:00.0340', ..., '2020-06-01 12:01:00.0340'])
ECEFvec = np.array([[-4.46e+03 -5.22e+03  1.75e-04], ..., [ 2.73e+03  2.08e+03 -6.02e+03]])

finch = Satellite()
ECIvec = finch.getECI(posData=ECEFvec, timeData=UTCTimeData)
```

#### .getECEF(**kwargs)
**Description**  
The `getECEF` method uses the earth rotation angles to convert ECI data to ECEF data while instantiating *ECEFdata*.  
**Parameters**  
None  
**\*\*kwargs**  
`posData`: ndarray array of shape (n,3) with columns of X, Y, Z position data assumed ECI.  
`timeData`: ndarray of shape (n,) containing `datetime.datetime` objects in UTC.  
**Returns**  
*ECEFdata* : ndarray of shape (n,3) with columns of X, Y, Z ECEF position data.  
**Usage**  
The *ECIdata* and *times* attributes must be initiated or passed in as **kwargs.  
**Example**  
```python
from celest import Satellite
import numpy as np

UTCTimeData = np.array(['2020-06-01 12:00:00.0340', ..., '2020-06-01 12:01:00.0340'])
ECIvec = np.array([[-4.46e+03 -5.22e+03  1.75e-04], ..., [ 2.73e+03  2.08e+03 -6.02e+03]])

finch = Satellite()
ECEFvec = finch.getECEF(posData=ECIvec, timeData=UTCTimeData)
```

#### .getAltAz(groundPos, **kwargs)
**Description**  
This method initiates the *alt* and *az* attributes of the GroundPosition object stored in the Satellites *gs* attribute. These attributes are both ndarray's of shape (n,) where rows correspond to the times in the *times* attribute. This method gains ground location information from the groundPos parameter.  
**Parameters**  
`groundPos`: GroundPosition object instantiated as per its documentation.  
**\*\*kwargs**  
`posData`: ndarray array of shape (n,3) with columns of X, Y, Z position data assumed ECI.  
`timeData`: ndarray of shape (n,) containing `datetime.datetime` objects in UTC.  
**Returns**  
*AltAz* : tuple of ndarrays of shape (n,) containing the computed altitude and azimuth data, (alt, az).  
**Usage**  
The function requires *ECEFdata* to be initiated or have the position and time dependencies passed in as **kwargs.  
If passing in a new GroundPosition object, it will become stored in the *gs* dictionary with the GroundPosition's *name* attribute as the key and the object as the value. If the object is already stored in *gs*, then its *alt* and *az* attributes will be updated.  
**Example**  
```python
from celest import Satellite, GroundPosition
import numpy as np

UTCTimeData = np.array(['2020-06-01 12:00:00.0340', ..., '2020-06-01 12:01:00.0340'])
ECIvec = np.array([[-4.46e+03 -5.22e+03  1.75e-04], ..., [ 2.73e+03  2.08e+03 -6.02e+03]])

toronto = GroundPosition(name="Toronto", coor=(43.662300, -79.394530))

finch = Satellite()
Alt, Az = finch.getAltAz(groundPos=toronto, posData=ECIvec, timeData=UTCTimeData)
```

#### .getNdrAng(groundPos, **kwargs)
**Description** 
This method initiates the *nadirAng* attribute of the GroundPosition object stored in the Satellite's *gs* attribute. The *nadirAng* attribute is a ndarray of shape (n,) with the nadir-LOS angles. This method gains ground location information from the groundPos parameters.  
**Parameters**  
`groundPos`: GroundPosition object instantiated as per its documentation.  
**\*\*kwargs**  
`posData`: ndarray array of shape (n,3) with columns of X, Y, Z position data assumed ECI.  
`timeData`: ndarray of shape (n,) containing `datetime.datetime` objects in UTC.  
**Returns**  
*nadirAng*: ndarray of shape (n,) with nadir-LOS angle data.  
**Usage**  
The function requires *ECEFdata* to be initiated or have the position and time dependencies passed in as **kwargs.  
If passing in a new GroundPosition object, it will become stored in the *gs* dictionary with the GroundPosition's *name* attribute as the key and the object as the value. If the object is already stored in *gs*, then its *nadirAng* attribute will be updated.  
**Example**  
```python
from celest import Satellite, GroundPosition
import numpy as np

UTCTimeData = np.array(['2020-06-01 12:00:00.0340', ..., '2020-06-01 12:01:00.0340'])
ECIvec = np.array([[-4.46e+03 -5.22e+03  1.75e-04], ..., [ 2.73e+03  2.08e+03 -6.02e+03]])

toronto = GroundPosition(name="Toronto", coor=(43.662300, -79.394530))

finch = Satellite()
NdrAng = finch.getNdrAng(groundPos=toronto, posData=ECIvec, timeData=UTCTimeData)
```

#### .saveData(fileName, delimiter)
**Description**  
This method saves the class data to the local working directory.  
**Parameters**  
`fileName`: str representing the output file name. The output file can be a .txt or .csv file.  
`delimiter`: str of length 1 representing the feild delimiter for the output file.  
**Returns**  
None  
**Usage**  
It is recommended to use a tab delimiter for .txt files and a comma delimiter for .csv files. Will return an error if fileName already exists in the current working directory.  
**Example**  
```python
from celest import Satellite, GroundPosition
import numpy as np

UTCTimeData = np.array(['2020-06-01 12:00:00.0340', ..., '2020-06-01 12:01:00.0340'])
ECIvec = np.array([[-4.46e+03 -5.22e+03  1.75e-04], ..., [ 2.73e+03  2.08e+03 -6.02e+03]])

toronto = GroundPosition(name="Toronto", coor=(43.662300, -79.394530))

finch = Satellite()
finch.getAltAz(groundPos=toronto, posData=ECIvec, timeData=UTCTimeData)
finch.saveData(fileName="data.csv", delimiter=",")
```

### GroundPosition Class
The GroundPosition class represents a ground location to store geographically specific orbital position representations. The class stores its data in seven instance variables:

1. **self. name**: String representing a tag for identifying the ground station location.
1. **self.coor**: Tuple of ground position coordinates as (latitude, longitude) in decimal degrees.
1. **self.radius**: Radius of earths surface at *coor* given in km.
1. **self.alt**: ndarray of shape (n,) containing altitude data instantiated by the Satellite `.getAltAz()` method.
1. **self. az**: ndarray of shape (n,) containing azimuth data instantiated by the Satellite `.getAltAz()` method.
1. **self.nadirAng**: ndarray of shape (n,) containing nadir angle data instantiated by the Satellite `.getNdrAng()` method.
1. **self.length**: Length n of data attributes.

Note that the GoundPosition object is completely instantiated through the initial class declaration. An example innitiation of a GroundPosition object is given by the following:

```python
from celest import GroundPosition

toronto = GroundPosition(name="Toronto", coor=(43.662300, -79.394530))
```

There exists only one method of the GroundPosition class that is used internally to calculate the Earths radius at the position defined by *coor*. This method approximates the earth as an ellipsoid using the World Geodetic System WGS84 to account for changing radius with latitude.

#### .getRadius(coor)
**Description**  
This method uses the World Geodetic System WGS84 to calculate the radius at the given coordinate position.  
**Parameters**  
`coor`: Tuple of ground position coordinates as (latitude, longitude) in decimal degrees.  
**Returns**  
*radius*: Float representing the radius given in km.  
**Usage**  
This method is designed to be used internally to instantiate the *radius*.  
**Example**  
None
