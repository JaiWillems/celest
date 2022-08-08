Tutorials
=========

.. contents:: Contents
   :depth: 1
   :local:

Time Conversions
----------------

Celest provides the :class:`Time` class that can be used for various time transformations. The :class:`Time` class can
be imported from the time module by the following import statement:

.. code-block:: python

   from celest.time import Time

We begin by importing the time data into the program. Since celest takes in time series data as NumPy arrays, we
recommend using NumPy's :py:func:`loadtxt` or :py:func:`genfromtxt` functions. The time arrays should contain `float`
time values.

.. code-block:: python

   import numpy as np

   julian = np.loadtxt(`julian.txt`)

Using the time array, we can initialize the :class:`Time` class. The input time shall be some julian variant. If times
not in the JD2000 epoch are passed in, an offset that can be added to the input time to make it JD2000 must be passed
in.

.. code-block:: python

   # Initialize Time object using J2000 data.
   time = Time(julian=julian)

   # Initialize Time object using Mean Julian Date data.
   time = Time(julian=julian, offset=2400000.5)

Different time representations can then be accessed through various methods.
Some examples are shown below.

.. code-block:: python

   # Get Greenwich Apparent Sidereal Time.
   gast = time.gast()

   # Get mean solar time.
   mst = time.mean_solar_time(longitude=-79.3832)

   # Get datetime representations.
   datetime_data = time.datetime()

Position Conversions
--------------------

Celest provides a suite for converting between different coordinate representations. For conversions, the most important
classes are the frame specific classes that hold data for a single coordinate frame and the :class:`Coordinates` class
to handle frame conversions. These different classes can be imported from the coordinate module using the following
statement:

.. code-block:: python

   from celest.coordinate import AzEl, GCRS, ITRS, Coordinates

We begin by importing the initial coordinate data into the program. Since celest takes in position data as NumPy arrays,
we recommend using NumPy's :py:func:`loadtxt` or :py:func:`genfromtxt` functions. The position arrays should contain
`float` values.

.. code-block:: python

   import numpy as np

   julian_data = np.loadtxt(`julian.txt`)
   gcrs_position_data = np.loadtxt(`gcrs_position.txt`)

Using the loaded position array in the GCRS frame, we can initialize a :class:`GCRS` object.

.. code-block:: python

   from celest import units as u

   gcrs_position = GCRS(
      julian=julian_data,
      x=gcrs_position_data[:, 0],
      y=gcrs_position_data[:, 1],
      z=gcrs_position_data[:, 2],
      unit=u.km
   )

We can then create a :class:`Coordinates` object that holds the position object and use it to convert to different
frames. Most frames require no additional arguments for conversions, however, some frames require a location passed in.

.. code-block:: python

   from celest.coordinate import GroundLocation

   # Create the Coordinates object.
   coordinates = Coordinates(gcrs_position)

   # Convert to ITRS.
   itrs_position = coordinates.convert_to(ITRS)

   # Convert to AzEl (requires a ground location).
   toronto = GroundLocation(
      latitude=43.6532,
      longitude=-79.3832,
      height=0.76,
      angular_unit=u.deg,
      length_unit=u.km
   )
   azel_position = coordinates.convert_to(AzEl, location=toronto)

Once the conversions are complete, the coordinate data can be accessed through a series of attributes. Such data is
stored as :class:`Quantity` objects allowing for easy unit conversions.

.. code-block:: python

   # Get the ITRS coordinates as Quantity objects.
   x = itrs_position.x
   y = itrs_position.y
   z = itrs_position.z

   # Get the AzEl coordinates as Quantity objects in degrees.
   az = azel_position.az.to(u.deg)
   el = azel_position.el.to(u.deg)

Window Generation Workflow
--------------------------

One of the primary goals of Celest is to provide tools for satellite planning. This includes window generation and
scheduling. The window generation workflow can be broken down into three stages:

#. Processing data,
#. Specifying ground locations, and
#. Generating window data.

We begin by importing the time and position data to initialize the :class:`Satellite` class.

.. code-block:: python

   from celest.coordinates import GCRS
   from celest.satellite import Satellite
   from celest import units as u
   import numpy as np

   julian_data = np.loadtxt(`julian.txt`)
   gcrs_position_data = np.loadtxt(`gcrs_position.txt`)
   gcrs_velocity_data = np.loadtxt(`gcrs_velocity.txt`)

    gcrs_position = GCRS(
        julian=julian_data,
        x=gcrs_position_data[:, 0],
        y=gcrs_position_data[:, 1],
        z=gcrs_position_data[:, 2],
        unit=u.km
    )
    gcrs_velocity = GCRS(
        julian=julian_data,
        x=gcrs_velocity_data[:, 0],
        y=gcrs_velocity_data[:, 1],
        z=gcrs_velocity_data[:, 2],
        unit=u.m/u.s
    )

   satellite = Satellite(position=gcrs_position, velocity=gcrs_velocity)

To generate window data, we need to specify the ground locations involved in the encounters. If we want to image the
Canadian cities of Toronto, Calgary, and Vancouver, we need to initialize three :class:`GroundLocation` objects.

.. code-block:: python

   from celest.coordinates import GroundLocation

   toronto = GroundLocation(
      latitude=43.6532,
      longitude=-79.3832,
      height=0.76,
      angular_unit=u.deg,
      length_unit=u.km
   )
   calgary = GroundLocation(
      latitude=51.0486,
      longitude=-114.0708,
      height=1.045,
      angular_unit=u.deg,
      length_unit=u.km
   )
   vancouver = GroundLocation(
      latitude=49.2827,
      longitude=-123.1207,
      height=0.0,
      angular_unit=u.deg,
      length_unit=u.km
   )

We can now create window data using the :py:func:`generate_vtws` function.

.. code-block:: python

   from celest.window import generate_vtws, Lighting

   # Generate window data for Toronto.
   toronto_vtws = generate_vtws(
      satellite=satellite,
      location=toronto,
      vis_threshold=10,
      lighting=Lighting.DAYTIME
   )

   # Generate window data for Calgary.
   calgary_vtws = generate_vtws(
      satellite=satellite,
      location=calgary,
      vis_threshold=10,
      lighting=Lighting.DAYTIME
   )

   # Generate window data for Vancouver.
   vancouver_vtws = generate_vtws(
      satellite=satellite,
      location=vancouver,
      vis_threshold=10,
      lighting=Lighting.DAYTIME
   )

More About Windows
^^^^^^^^^^^^^^^^^^

Using Celest effectively requires a good understanding of the concept of windows and how to define them. In their most
basic form, a window is defined as a time where there exists a line of sight between the satellite and the ground
location. Theoretically, this line of sight should exist whenever the satellite has a positive elevation angle with
respect to the ground location. However, do to obstructions such as topography or structures, the line of sight may only
exist when the satellite's elevation angle is above some threshold: the visibility threshold. This threshold can be
passed into the :py:func:`generate_vtws` function as the `vis_threshold` argument.

The other primary parameter of a window is the lighting condition. Some types of encounters may require specific
lighting conditions. For example, an imaging encounter may require the satellite to be in the daytime, while a data
transfer encounter may occur at any time of day. The :class:`Lighting` enum provides a list of possible lighting
conditions that can be passed into the :py:func:`generate_vtws` function under the `lighting` argument.

Scheduling Workflow
-------------------

One of the primary goals of Celest is to provide tools for satellite planning. This includes window generation and
scheduling. The scheduling workflow can be broken down into four stages:

#. Processing data,
#. Specifying ground locations,
#. Defining encounter requests, and
#. Generating a satellite schedule.

We begin by importing the time and position data to initialize the :class:`Scheduler` class.

.. code-block:: python

   from celest.coordinates import GCRS
   from celest.satellite import Satellite
   from celest.schedule import Scheduler
   from celest import units as u
   import numpy as np

   julian_data = np.loadtxt(`julian.txt`)
   gcrs_position_data = np.loadtxt(`gcrs_position.txt`)
   gcrs_velocity_data = np.loadtxt(`gcrs_velocity.txt`)

    gcrs_position = GCRS(
        julian=julian_data,
        x=gcrs_position_data[:, 0],
        y=gcrs_position_data[:, 1],
        z=gcrs_position_data[:, 2],
        unit=u.km
    )
    gcrs_velocity = GCRS(
        julian=julian_data,
        x=gcrs_velocity_data[:, 0],
        y=gcrs_velocity_data[:, 1],
        z=gcrs_velocity_data[:, 2],
        unit=u.m/u.s
    )

   satellite = Satellite(position=gcrs_position, velocity=gcrs_velocity)
   scheduler = Scheduler(satellite=satellite, vis_threshold=10)

To generate window data, we need to specify the ground locations involved in the encounters. If we want to image the
Canadian cities of Toronto, Calgary, and Vancouver, we need to initialize three :class:`GroundLocation` objects.

.. code-block:: python

   from celest.coordinates import GroundLocation

    toronto = GroundLocation(
        latitude=43.6532,
        longitude=-79.3832,
        height=0.76,
        angular_unit=u.deg,
        length_unit=u.km
    )
    calgary = GroundLocation(
        latitude=51.0486,
        longitude=-114.0708,
        height=1.045,
        angular_unit=u.deg,
        length_unit=u.km
    )
    vancouver = GroundLocation(
        latitude=49.2827,
        longitude=-123.1207,
        height=0.0,
        angular_unit=u.deg,
        length_unit=u.km
    )

We now need to specify the encounters that we want to schedule; this is done through adding a request to the
`scheduler`. A request is a desired encounter between the satellite and a ground location that meets certain criteria.
The scheduler then generates the possible windows for the request and works to schedule one of the windows that meets
all criteria. There are various parameters that define a request:

#. The ground location associated with the satellite-to-ground encounter,
#. The deadline for when the encounter should occur before,
#. The duration that the encounter should be scheduled for,
#. The priority of the request over other requests,
#. The quality of the encounter which defines how nadir the encounter occurs,
#. A specific look angle for the encounter if desired, and
#. The lighting condition for the encounter.

Note that a request being added to the scheduler does not guarantee that the request will be scheduled. The scheduler
will search for the most optimal schedule without conflicts. Therefore, if significant conflicts exist (violations of
the request criteria or overlaping of scheduled tasks), the scheduler may not be able to find a schedule that both
schedules all tasks and has no conflicts, in this case, the lower priority request will be ignored.

.. code-block:: python

   from celest.encounter import Lighting

   schedule.add_request(
      location=toronto,
      deadline=2460467,
      duration=30,
      priority=5,
      quality=1,
      look_ang=None,
      lighting=Lighting.DAYTIME
   )
   schedule.add_request(
      location=calgary,
      deadline=2460467,
      duration=30,
      priority=3,
      quality=1,
      look_ang=None,
      lighting=Lighting.DAYTIME
   )
   schedule.add_request(
      location=vancouver,
      deadline=2460467,
      duration=30,
      priority=1,
      quality=1,
      look_ang=None,
      lighting=Lighting.DAYTIME
   )

The schedule can now be generated and the data exported. Note that the scheduler may take significant time to search
for an optimal solution. The parameters that govern the search are will affect this time; these parameters are:

#. The number of schedule iterations to make during the search,
#. The annealing coefficient for the simulated annealing process that governs the acceptance probability of a non-improving solution, and
#. The reactivity factor that determines the effect of past iterations on the current search for an optimal solution.

.. code-block:: python

   # Determine a feasible schedule.
   schedule = scheduler.generate(max_iter=100, annealing_coeff=0.8, react_factor=0.5)
   schedule.save_text_file(file_name="canada_imaging_schedule")
