Quickstart
==========

.. contents:: Contents
   :depth: 1
   :local:

Installation using pip
----------------------

To install `celest` using `pip`, run:

.. code-block::

   pip install celest

Getting familiar with the library
---------------------------------

Celest is broken into various modules:

* :mod:`celest.units` handles units and their conversions.
* :mod:`celest.coordinates` handles coordinate systems and their conversions.
* :mod:`celest.satellite` handles satellite specific calculations.
* :mod:`celest.encounter` generates and stores encounter windows.
* :mod:`celest.schedule` handles scheduling of satellite requests.

A typical workflow for using Celest begins by importing satellite position and velocity data.
Coordinate objects can then be initialized to be used in the rest of the library.

.. code-block:: python

   from celest.coordinates import GCRS
   from celest import units as u
   import numpy as np

   # Load the data.
   julian = np.loadtxt('julian.txt')
   position = np.loadtxt('gcrs_position.txt')
   velocity = np.loadtxt('gcrs_velocity.txt')

   # Initialize the coordinate objects.
   gcrs_position = GCRS(
      julian=julian_data,
      x=position[:, 0],
      y=position[:, 1],
      z=position[:, 2],
      unit=u.km
   )
   gcrs_velocity = GCRS(
      julian=julian_data,
      x=velocity[:, 0],
      y=velocity[:, 1],
      z=velocity[:, 2],
      unit=u.m/u.s
   )

Coordinates can then be converted into other frames using the
:class:`celest.coordinates.Coordinate` class.

.. code-block:: python

   from celest.coordinates import Coordinate, GroundLocation, ITRS

   # Create a coordinate object.
   satellite_coordinates = Coordinate(gcrs_position)

   # Convert the position to ITRS.
   itrs_position = satellite_coordinates.convert_to(ITRS)

   # Convert the position into the horizontal frame.
   toronto = GroundLocation(
      latitude=43.65,
      longitude=-79.38,
      height=0,
      angular_unit=u.deg,
      length_unit=u.km
   )
   horizontal_position = satellite_coordinates.convert_to(
      horizontal_frame,
      location=toronto
   )

The satellite coordinate representations can be used to create a :class:`celest.satellite.Satellite`
object and perform satellite specific calculations such as determining the satellite's attitude towards a ground
location.

.. code-block:: python

   from celest.satellite import Satellite

   # Create a satellite object.
   satellite = Satellite(itrs_position, itrs_velocity)

   # Get the satellite's attitude towards Toronto.
   to_toronto_attitude = satellite.attitude(location=toronto)

Using a :class:`celest.satellite.Satellite` instantiation, possible encounter
times for a satellite-ground-position pair can be generated.


.. code-block:: python

   from celest.encounter import Lighting, generate_vtw

   # Generate daylight encounters.
   toronto_daylight_encounters = generate_vtw(
      satellite=satellite,
      location=toronto,
      vis_threshold=10,
      lighting=Lighting.DAYTIME
   )

   # Night time encounters can be determined using `Lighting.NIGHTTIME`.
   # No lighting constraint is attained using `Lighting.ANYTIME`.

   # The daylight encounters can then be saved as a text file.
   toronto_daylight_encounters.save_text_file(
      file_name='toronto_daylight_encounters'
   )

A user may wish to schedule a series of different imaging encounters for a
particular satellite. The window generation for each possible imaging request
and the scheduling of such tasks is handled by the :class:`Scheduler` class.

.. code-block:: python

   from celest.schedule import Scheduler

   # Begin by defining the ground locations to image.
   toronto = GroundLocation(
      latitude=43.65,
      longitude=-79.38,
      height=0.076
      length_unit=u.km
      angular_unit=u.deg)
   north_bay = GroundLocation(
      latitude=46.31,
      longitude=-79.46,
      height=0.193
      length_unit=u.km
      angular_unit=u.deg)
   sudbury = GroundLocation(
      latitude=46.49,
      longitude=-80.99,
      height=0.348
      length_unit=u.km
      angular_unit=u.deg)
   mississauga = GroundLocation(
      latitude=43.59,
      longitude=-79.64,
      height=0.156
      length_unit=u.km
      angular_unit=u.deg)

   # Initialize the scheduling object.
   scheduler = Scheduler(satellite=satellite, vis_threshold=10)

   # Add imaging requests to schedule.
   schedule.add_request(
      location=toronto,
      deadline=2460467,
      duration=30,
      priority=1,
      look_ang=None
      lighting=Lighting.DAYLIGHT
   )
   schedule.add_request(
      location=north_bay,
      deadline=2460467,
      duration=30,
      priority=1,
      look_ang=None
      lighting=Lighting.DAYLIGHT
   )
   schedule.add_request(
      location=sudbury,
      deadline=2460467,
      duration=30,
      priority=4,
      look_ang=None
      lighting=Lighting.DAYLIGHT
   )
   schedule.add_request(
      location=mississauga,
      deadline=2460467,
      duration=30,
      priority=5,
      look_ang=None
      lighting=Lighting.DAYLIGHT
   )

   # Determine a feasible schedule.
   schedule = schedule.generate(
      max_iter=100,
      annealing_coeff=0.8,
      react_factor=0.5
   )
   schedule.save_text_file(file_name="ontario_imaging_schedule")
