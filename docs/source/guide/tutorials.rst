Tutorials
=========

Position and Time Conversions
-----------------------------

Celest allows for various time and position representations that may
prove useful in a mission planning context. Here, we will explore
the basic interactions with the :class:`Time` and :class:`Coordinate` classes.
For details on all possible conversions, refer to the :ref:`Time <Time Class>`
and :ref:`Coordinate <Coordinate Class>` documentation.

.. code-block::

   import numpy as np
   from celest.satellite import Time, Coordinate

We begin by importing the desired position and time data into the program.
Since Celest takes in data as NumPy arrays, it is recommended to use NumPy's
:py:func:`loadtxt` or :py:func:`genfromtxt` functions. The data should be
loaded such that all array values are of a type `float`.

.. code-block::

   # Load time and data from external sources.
   julian = np.loadtxt('julian.txt')
   gcrs = np.loadtxt('gcrs.txt')

We can then instantiate the :class:`Time` class. The input data
should be Julian time data in the J2000 epoch. The appropriate offset can be
applied using the `offset` keyword to incorporate data from other epochs. This
offset will be added to the input data when stored in the :class:`Time` object.

.. code-block::

   # Initialize Time object using J2000 data.
   time = Time(julian)

   # Initialize Time object using Mean Julian Date data.
   time = Time(julian, offset=2400000.5)

Different time representations can then be accessed through various methods.
Some examples are shown below.

.. code-block::

   # Get Greenwich Apparent Sidereal Time.
   gast = time.gast()

   # Get Mean Hour Angle.
   mha = time.mean_solar_time(longitude=-79.3832)

   # Get datetime representaions.
   datetime_data = time.datetime()

The :class:`Coordinate` class can be instantiated using a data array of
positions, a frame specifier, and time information associated with the
position data. Three input position frames are currently supported by the
:class:`Coordinate` class: the geocentric celestial reference system (GCRS),
international terrestrial reference system (ITRS), and geodetic reference
system (GRS).

.. code-block::

   # Initialize Coordinate object using GCRS data.
   position = Coordinate(gcrs, frame='gcrs', time=julian, offset=0)

   # Initialize Coordinate object using ITRS data.
   position = Coordinate(itrs, frame='itrs', time=julian, offset=0)

   # Initialize Coordinate object using GCRS data.
   position = Coordinate(geo, frame='geo', time=julian, offset=0)

The :class:`Coordinate` class inherits the :class:`Time` class; as a result,
any time methods can be accessed using the :class:`Coordinate` class.

Similar to the :class:`Time` class, different coordinate representations can be
accessed through various methods. Some examples are shown below.

.. code-block::

   from celest.encounter import GroundPosition
   toronto = GroundPosition(latitude=43.6532, longitude=-79.3832, height=0.076)

   # Get horizontal coordinates.
   alt, az = position.horizontal(location=toronto)

   # Get ITRS data.
   itrs_position = position.itrs()

   # Get GRS data with ISO6709 formatted output strings.
   geo_position = position.geo(iso=True)

Notice that some methods require using a :class:`GroundPosition` object as
a parameter to specify a ground location. The :class:`GroundPosition` object
can be imported from the encounter module.

Window Generation Workflow
--------------------------

The primary window generation workflow can be broken down into three stages:

#. Import and prepare time and position data,
#. Specify ground locations, and
#. Generate and save windows.

The first step is to import and prepare the time and position data. This
includes setting up the :class:`Satellite` object that holds the necessary but
insufficient information to generate the desired windows.

.. code-block::

   from celest.satellite import Satellite
   from celest.encounter import GroundPosition, windows
   import numpy as np

   # Load the data.
   julian = np.loadtxt('julian.txt')
   gcrs = np.loadtxt('gcrs.txt')

   # Initialize satellite representation.
   satellite = Satellite(position=gcrs_pos, velocity=gcrs_vel, frame='gcrs', time=julian, offset=0)

Next, we specify the ground locations for which we wish to generate windows. To
accomplish this, we define a :class:`GroundPosition` object for each location
we wish to encounter. If various encounter types for one location are desired,
only one :class:`GroundPosition` object is required.

.. code-block::

   # Define ground position.
   toronto = GroundPosition(latitude=43.6532, longitude=-79.3832, height=0.076)
   saskatoon = GroundPosition(latitude=52.1579, longitude=-106.6702, height=0.482)

We are now ready to generate windows. The :py:func:`windows.generate_vtw` function
takes a satellite and ground location as an input and will populate a
:class:`VTWHandling` object with visible time windows defined by the
`vis_threshold` and `lighting` keywords.


The `vis_threshold` keyword defines the minimum satellite elevation angle as
seen from the ground location that will allow for satellite-ground
interactions. The `lighting` keyword allows us to specify the lighting
conditions of the visible time windows.

.. code-block::

   # Generate ground location windows.
   toronto_IMG_windows = windows.generate_vtw(satellite=satellite, location=toronto, vis_threshold=10, lighting=1)
   toronto_GL_windows = windows.generate_vtw(satellite=satellite, location=toronto, vis_threshold=10, lighting=0)

   # Save satellite encounter windows.
   toronto_IMG_windows.save(fname="toronto_IMG_windows.csv", delimiter=",")
   toronto_DL_windows.save(fname="toronto_DL_windows.csv", delimiter=",")

Scheduling Workflow
-------------------

The scheduling workflow can be broken down into three stages:

#. Import and prepare data,
#. Specify ground locations,
#. Define encounter requests,
#. Generate and save windows.

The first step is to import and prepare the time, position, and velocity data.
This includes setting up the :class:`Satellite` object that holds necessary
information for the request scheduling.

.. code-block::

   from celest.satellite import Satellite
   import numpy as np

   # Load the data.
   julian = np.loadtxt('julian.txt')
   gcrs_pos = np.loadtxt('gcrs_pos.txt')
   gcrs_vel = np.loadtxt('gcrs_vel.txt')

   # Initialize satellite representation.
   satellite = Satellite(position=gcrs, frame='gcrs', time=julian, offset=0)

Next, we specify the ground locations for which we wish to image. To
accomplish this, we define a :class:`GroundPosition` object for each location
we wish to encounter. If various encounter types for one location are desired,
only one :class:`GroundPosition` object is required.

.. code-block::

   from celest.encounter import GroundPosition

   # Define ground position.
   toronto = GroundPosition(latitude=43.65, longitude=-79.38, height=0.076)
   north_bay = GroundPosition(latitude=46.31, longitude=-79.46, height=0.193)
   sudbury = GroundPosition(latitude=46.49, longitude=-80.99, height=0.348)
   mississauga = GroundPosition(latitude=43.59, longitude=-79.64, height=0.156)

We now need to define the requests. A request specifies an imaging encounter
that we wish to make. The scheduler will then take these requests and determine
when each request should be fulfilled. Note that there is no gaurentee that
all requests will be scheduled if insufficient opportunities or conflicts
exist.

Before defining the requests, we initialize the :class:`Schedule` class by
specifying the satellite that will participate in fulfilling the requests and
the visibility threshold that will be used to determine feasible imaging
encounters. The visibility threshold is the minimum elevation angle of the
satellite as seen from the desired ground location that allows for a viable
imaging encounter.

For each request, various parameters can be specified.

.. code-block::

   from celest.schedule import Schedule

   # Initialize the scheduling object.
   schedule = Schedule(satellite=satellite, vis_threshold=10)

   # Add imaging requests we want to schedule.
   schedule.add_request(location=toronto, deadline=2460467, duration=30, priority=1, look_ang=None)
   schedule.add_request(location=north_bay, deadline=2460467, duration=30, priority=1, look_ang=None)
   schedule.add_request(location=sudbury, deadline=2460467, duration=30, priority=4, look_ang=None)
   schedule.add_request(location=mississauga, deadline=2460467, duration=30, priority=5, look_ang=None)

The schedule can now be generated and the data exported.

.. code-block::

   # Determine a feasible schedule.
   schedule_out = schedule.generate(max_iter=100, annealing_coeff=0.8, react_factor=0.5)
   schedule_out.save(fname="ontario_imaging_schedule.csv", delimiter=",")
