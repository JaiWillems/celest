Satellite Module
================

Time Class
----------

The `Time` class allows for flexible and convenient time specifications. It
takes in Julian date data in the J2000 epoch or another epoch with the
appropriate offset to convert to the J2000 epoch. Example usage of the `Time`
class to carry out conversions can be seen on the tutorials page
`here <Position and Time Conversions>`.

.. autoclass:: celest.satellite.Time
   :members:
   :show-inheritance:
   :noindex:

Coordinate Class
----------------

The `Coordinate` class is used to represent a satellite position and can be
used to carry out coordinate transformations between different coordinate
frames. The Example usage of the `Coordinate` class to carry out conversions
can be seen in the tutorials page `here <Position and Time Conversions>`.

.. autoclass:: celest.satellite.Coordinate
   :members:
   :show-inheritance:
   :noindex:

Satellite Class
---------------

The `Satellite` class is used to hold all necessary information for a given
satellite such that it can be passed into the window generation functionality
for mission planning. Example usage of the `Satellite` in a mission planning
capacity can be seen on the tutorials page
`here <Window Generation Workflow>`.

.. autoclass:: celest.satellite.Satellite
   :members:
   :show-inheritance:
   :noindex:
