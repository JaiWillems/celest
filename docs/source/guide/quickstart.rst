Quickstart
==========

Installation using pip
----------------------

If you use `pip`, you can install Celest with:

.. code-block::

   pip install celest

Getting familiar with the library
---------------------------------

Celest is broken primarily into two modules: (1) :mod:`celest.satellite` module
contains all functions related to satellite placement and dynamics such as
coordinate and time conversions. (2) :mod:`celest.encounter` module contains
all functions related to encounter calculations between a satellite and ground
location.

A typical workflow for using Celest is to begin by importing satellite time and
position data as NumPy arrays which can then be used to initialize the
:class:`Satellite` class. The :class:`Satellite` object can then be coupled
with encounter information such as a ground location and encounter parameters
to determine encounter opportunities.

.. code-block::

   from celest.satellite import Time, Coordinate, Satellite
   from celest.schedule import Schedule
   from celest.encounter import GroundPosition, windows
   import numpy as np

   # Load the data.
   julian = np.loadtxt('julian.txt')
   gcrs = np.loadtxt('gcrs.txt')

   # Initialize satellite representation.
   satellite = Satellite(position=gcrs, frame='gcrs', time=julian, offset=0)

   # Define ground position.
   toronto = GroundPosition(latitude=43.65, longitude=-79.38, height=0.076)

   # Generate ground location windows.
   toronto_IMG_windows = windows.generate_vtw(satellite=satellite, location=toronto, vis_threshold=10, lighting=1)
   toronto_GL_windows = windows.generate_vtw(satellite=satellite, location=toronto, vis_threshold=10, lighting=0)

   # Save satellite encounter windows.
   toronto_IMG_windows.save(fname="toronto_IMG_windows.csv", delimiter=",")
   toronto_DL_windows.save(fname="toronto_DL_windows.csv", delimiter=",")

We can also directly schedule a series of imaging requests.

.. code-block::

   # Start by defining ground locations.
   toronto = GroundPosition(latitude=43.65, longitude=-79.38, height=0.076)
   north_bay = GroundPosition(latitude=46.31, longitude=-79.46, height=0.193)
   sudbury = GroundPosition(latitude=46.49, longitude=-80.99, height=0.348)
   mississauga = GroundPosition(latitude=43.59, longitude=-79.64, height=0.156)

   # Initialize the scheduling object.
   schedule = Schedule(satellite=satellite, vis_threshold=10)

   # Add imaging requests we want to schedule.
   schedule.add_request(location=toronto, deadline=2460467, duration=30, priority=1, look_ang=None)
   schedule.add_request(location=north_bay, deadline=2460467, duration=30, priority=1, look_ang=None)
   schedule.add_request(location=sudbury, deadline=2460467, duration=30, priority=4, look_ang=None)
   schedule.add_request(location=mississauga, deadline=2460467, duration=30, priority=5, look_ang=None)

   # Determine a feasible schedule.
   schedule_out = schedule.generate(max_iter=100, annealing_coeff=0.8, react_factor=0.5)
   schedule_out.save(fname="ontario_imaging_schedule.csv", delimiter=",")
