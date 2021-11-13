Tutorials
=========

Position and Time Conversions
-----------------------------

Celest allows for a variety of time and position representaions that may
provide use in a mission planning context. Here, we will endevour to explore
the basic interactions with the Time and Coordinate classes. For details on all
possible conversions, refer to the :ref:`Time <Time Class>` and
:ref:`Coordinate <Coordinate Class>` documentation.

.. code-block::

   import numpy as np
   from celest.satellite import Time, Coordinate, Satellite
   from celest.encounter import GroundPosition, windows, attitudes

We begin by importing the desired position and time data into the program.
Since Celest takes in data as NumPy arrays, it is recommend to use NumPy's
`loadtxt` or `genfromtxt` functions. The data should be loaded such that all
array values are of a type `float`.

.. code-block::

   # Load time and data from external sources.
   julian = np.loadtxt('julian.txt')
   gcrs = np.loadtxt('gcrs.txt')

We can then go ahead instantiating our `Time` object. The input data should be
Julian time data in the J2000 epoch. To incorporate data from other epochs, the
appropriate offset can be applied using the `offset` keyword.

.. code-block::

   # Initialize Time object using J2000 data.
   time = Time(julian)

   # Initialize Time object using Mean Julian Date data.
   time = Time(julian, offset=2400000.5)

Different time representations can then be accessed using various methods. Some
examples are shown below.

.. code-block::

   # Get Greenwich Apparent Sidereal Time.
   gast = time.gast()

   # Get Mean Hour Angle.
   mha = time.mean_solar_time(longitude=-79.3832)

   # Get datetime representaions.
   datetime_data = time.datetime()

The `Coordinate` class can be instantiated using a data array of positions, a
frame specifier, and a `Time` object associated with the position data. three
input position frames are currently supported by the `Coordinate` class:
geocentric celestial reference system (GCRS), international terrestrial
reference system (ITRS), and geodetic reference system (GRS).

.. code-block::

   # Initialize Coordinate object using GCRS data.
   position = Coordinate(gcrs, frame='gcrs', time=time)

   # Initialize Coordinate object using ITRS data.
   position = Coordinate(itrs, frame='itrs', time=time)

   # Initialize Coordinate object using GCRS data.
   position = Coordinate(geo, frame='geo', time=time)

Similar to the `Time` class, different coordinate representations can be
accessed using various methods. Some examples are shown below.

..code-block::

   from celest.encounter import GroundPosition
   toronto = GroundPosition(latitude=43.6532, longitude=-79.3832, altitude=0.0)

   # Get horizontal coordinates.
   alt, az = position.horizontal(location=toronto)

   # Get ITRS data.
   itrs_position = position.itrs()

   # Get GRS data with ISO6709 formatted output strings.
   geo_position = position.geo(iso=True)

Notice that some methods require the use of a `GroundPosition` object as
parameter to specify a ground location. The `GroundPosition` object can be
imported from the encounter module.

Window Generation Workflow
--------------------------

The basic window generation workflow can be broken down into three stages:

#. Import and prepare time and position data,
#. Specify ground locations, and
#. Generate and save winows.

The first step is to import and prepare the time and position data. This
includes setting up the `Satellite` object that holds the necessary but not
sufficient information to generate the desired windows.

.. code-block::

   import numpy as np
   from celest.satellite import Time, Coordinate, Satellite
   from celest.encounter import GroundPosition, windows, attitudes

   # Load the data.
   julian = np.loadtxt('julian.txt')
   gcrs = np.loadtxt('gcrs.txt')

   # Initialize satellite representation.
   time = Time(julian, offset=0)
   position = Coordinate(gcrs, frame='gcrs', time=time)
   satellite = Satellite(position=position)

Next, we specify the ground locations that we wish to generate windows for. To
accomplish this, we define a `GroundPosition` object for each location we wish
to encounter. If various encounter types for one location are desired, only
one `GroundPosition` object is required.

.. code-block::

   # Define ground position.
   toronto = GroundPosition(latitude=43.65, longitude=-79.38)
   saskatoon = GroundPosition(latitude=52.13, longitude=-106.67)

We are now ready to generate windows. The `windows.generate` function takes
a satellite and ground location as an input and will populate a `Window` object
with possible encounter opportunities for the encounter defined by the `enc`
and `ang` keywords.

There are two encounter types that Celest currently supports: (1) imaging
encounters where the satellite is in view of the ground location, and (2) data
transmission encounters where the ground location is in view of the satellite.
The `enc` keyword specifies the type of encounter as either and imaging
(enc="image") or data transmission (enc="data_link") type.

The `ang` keyword defines the constraint angle that boarders a
viable/non-viable encounter region. The contraint angle type used for the
imaging encounters is off-nadir angle measured in increasing degrees from the
satellite's nadir to the ground location. The transmission encounters use the
altitude angle of the satellite as measured in increasing degrees above the
horizon (as seen from the ground location).

.. code-block::

   # Generate ground location windows.
   toronto_IMG_windows = windows.generate(satellite=satellite, location=toronto, enc="image", ang=30)
   toronto_GL_windows = windows.generate(satellite=satellite, location=toronto, enc="data_link", ang=10)

   # Save satellite encounter windows.
   toronto_IMG_windows.save_encounters(fname="toronto_IMG_windows.csv", delimiter=",")
   toronto_DL_windows.save_encounters(fname="toronto_DL_windows.csv", delimiter=",")