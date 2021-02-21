# Celest

This purpose of the Celest library is to provide a simple interface to orbital coordinate conversions at an unprecedented speed which allows users to specify the exact computations to be enacted. This library is designed for small satellite applications who desire efficiency over precision in orbital calculations.

## Satellite Class

Celest is broken into a single Satellite class representing the orbital object to store and compute orbital coordinate data. The class stores its data in six instance variables:

1. **self.times**: Orbital time dependencies in an (n,) shaped ndarray.
1. **self.ERAdata**: Earth rotation angles based on time dependencies in an (n,) shaped ndarray.
1. **self.ECIdata**: Earth centered inertial position data in an (n,3) shaped ndarray where the columns are the x, y, z position data.
1. **self.ECEFdata**: Earth centered earth fixed position data in an (n,3) shaped ndarray where the columns are the x, y, z position data.
1. **self.horizontal**: Altitude and azimuth data in an (n,2) shaped ndarray where the columns are the alt, az position data.
1. **self.length**: Length of data attributes represented as an int

These inerence variables are interfaced by the user through seven methods:

1. **timeData**: Instantiates *times* attribute with orbital time dependency.
1. **positionData**: Instantiates *ECIdata* or *ECEFdata* attribute with orbital position data.
1. **getERA**: Instantiates *ERAdata* attribute.
1. **getECI**: Instantiates *ECIdata* attribute.
1. **getECEF**: Instantiates *ECEFdata* attribute.
1. **getAltAz**: Instantiates *horizontal* attribute.
1. **saveData**: Save class data in local directory.

The guiding principle behind the design of the Satellite class is to be simple to use and efficient. The latter is gained by two means. First, the code uses simple geometric relationships within orbital mechanics to derive basic mathematical relations that can be expedited by vectored inputs and NumPy. Secondly, the methods are created such that the user has complete control over the enacted computations by only instantiating necessary instance variables. The restrictions on how to use these methods are outlined in the following subsections along with detailed explanations and example usage of each of the methods.

### .timeData(timeData)
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
from celest.satellite import *
import numpy as np

finch = Satellite()
UTCTimeData = np.array(['2020-06-01 12:00:00.0340', ..., '2020-06-01 12:01:00.0340'])
finch.timeData(timeData=UTCtimeData)
```


### .positionData(timeData, type)
**Description**  
The `positionData` method interfaces with the *ECIdata* or *ECEFdata* attributes to initialize position data of the satellite. The user must specify the data as either ECI or ECEF data through the `type` parameter.  
**Parameters**  
`posData`: ndarray of shape (n,3) with columns of xyz position data.  
`type`: Specifies the type of position data as either "ECI" or "ECEF".  
**Returns**  
None  
**Usage**  
The length of `posData` must be the same as the *length* attribute.  
**Example**  
```python
from celest.satellite import *
import numpy as np

finch = Satellite()
ECIvec = np.array([[-4.46e+03 -5.22e+03  1.75e-04], ..., [ 2.73e+03  2.08e+03 -6.02e+03]])
finch.positionData(posData=ECIvec, type="ECI")
```

### .getERA(**kwargs)
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
from celest.satellite import *
import numpy as np

finch = Satellite()
UTCTimeData = np.array(['2020-06-01 12:00:00.0340', ..., '2020-06-01 12:01:00.0340'])
finch.timeData(timeData=UTCtimeData)
ERAangles = finch.getERA()
```

The commands can be simplified as in the following,  

```python
from celest.satellite import *
import numpy as np

finch = Satellite()
UTCTimeData = np.array(['2020-06-01 12:00:00.0340', ..., '2020-06-01 12:01:00.0340'])
ERAangles = finch.getERA(timeData=UTCTimeData)
```

### .getECI(**kwargs)
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
from celest.satellite import *
import numpy as np

UTCTimeData = np.array(['2020-06-01 12:00:00.0340', ..., '2020-06-01 12:01:00.0340'])
ECEFvec = np.array([[-4.46e+03 -5.22e+03  1.75e-04], ..., [ 2.73e+03  2.08e+03 -6.02e+03]])

finch = Satellite()
ECIvec = finch.getECI(posData=ECEFvec, timeData=UTCTimeData)
```

### .getECEF(**kwargs)
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
from celest.satellite import *
import numpy as np

UTCTimeData = np.array(['2020-06-01 12:00:00.0340', ..., '2020-06-01 12:01:00.0340'])
ECIvec = np.array([[-4.46e+03 -5.22e+03  1.75e-04], ..., [ 2.73e+03  2.08e+03 -6.02e+03]])

finch = Satellite()
ECEFvec = finch.getECEF(posData=ECIvec, timeData=UTCTimeData)
```

### .getAltAz(obsCoor, radius, **kwargs)
**Description**  
This method initiates the *horizontal* attribute of the Satellite object. This attribute is a ndarray of shape (n,2) where the first column is the altitude data and the second is the azimuth data. This method takes in an observers latitude and longitude in degrees as well as the surface radius.  
**Parameters**  
`obsCoor`: tuple of floats specifying the observer position in degrees (latitude, longitude).  
`radius`: int or float representing the fixed radius of Earth.  
**\*\*kwargs**  
`posData`: ndarray array of shape (n,3) with columns of X, Y, Z position data assumed ECI.  
`timeData`: ndarray of shape (n,) containing `datetime.datetime` objects in UTC.  
**Returns**  
*horizontal*: ndarray of shape (n,2) with columns of alt, az position data.  
**Usage**  
The function requires *ECEFdata* to be initiated or have the position and time dependencies passed in as **kwargs.  
**Example**  
```python
from celest.satellite import *
import numpy as np

UTCTimeData = np.array(['2020-06-01 12:00:00.0340', ..., '2020-06-01 12:01:00.0340'])
ECIvec = np.array([[-4.46e+03 -5.22e+03  1.75e-04], ..., [ 2.73e+03  2.08e+03 -6.02e+03]])

finch = Satellite()
AltAz = finch.getAltAz(obsCoor=(43.662300, -79.394530), radius=6371, posData=ECIvec, timeData=UTCTimeData)
```

### .saveData(fileName)
**Description**  
This method saves the class data to the local working directory.  
**Parameters**  
`fileName`: str representing the file name.  
**Returns**  
None  
**Usage**  
All instance variables must be initiated.  
**Example**  
```python
from celest.satellite import *
import numpy as np

UTCTimeData = np.array(['2020-06-01 12:00:00.0340', ..., '2020-06-01 12:01:00.0340'])
ECIvec = np.array([[-4.46e+03 -5.22e+03  1.75e-04], ..., [ 2.73e+03  2.08e+03 -6.02e+03]])

finch = Satellite()
# getAltAz() in this way will initiate all instance variables
finch.getAltAz(obsCoor=(43.662300, -79.394530), radius=6371, posData=ECIvec, timeData=UTCTimeData)
finch.saveData(filename="data.txt")
```
