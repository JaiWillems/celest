Satellite Module
================

Time Class
----------

The :class:`Time` class allows for flexible and convenient time specifications.
It takes in Julian data in the J2000 epoch or another epoch with an appropriate
offset. This offset will be added to the input data when stored in the
:class:`Time` object. Example usage of the :class:`Time` class to carry out
conversions can be seen on the tutorials page
:ref:`here <Position and Time Conversions>`.

.. autoclass:: celest.satellite.Time
   :members:
   :show-inheritance:
   :noindex:

Coordinate Class
----------------

The :class:`Coordinate` class is used to represent a satellite's position and
can be used to carry out coordinate transformations between different
frames. The Example usage of the :class:`Coordinate` class to carry out
conversions can be seen on the tutorials page
:ref:`here <Position and Time Conversions>`.

.. autoclass:: celest.satellite.Coordinate
   :members:
   :show-inheritance:
   :noindex:

Satellite Class
---------------

The :class:`Satellite` class is used to hold all necessary information for a
given satellite such that it can be passed into the :py:func:window.generate`
function for mission planning. Example usage of the :class:`Satellite` in
a mission planning capacity can be seen on the tutorials page
:ref:`here <Window Generation Workflow>`.

.. autoclass:: celest.satellite.Satellite
   :members:
   :show-inheritance:
   :noindex:
