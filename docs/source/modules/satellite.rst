Satellite
=========

.. contents:: Contents
   :depth: 1
   :local:

Satellite Overview
------------------

The satellite module provides the :class:`Satellite` class for satellite specific computations. The following sections
will detail the functionality of the :class:`Satellite` class. Additionally, a more in depth discussion of the attitude
and look-angle calculations will be provided.

Satellite Class
---------------

The :class:`Satellite` class represents an Earth observation satellite and allows for satellite specific computations.

.. note::
   The :class:`Satellite.attitude()` method requires both satellite position and velocity to be defined.

The :class:`Satellite` class can be imported via the following:

.. code-block:: python

   from celest.satellite import Satellite

.. autoclass:: celest.satellite.Satellite
   :members:
   :noindex:

Satellite Attitude Definition
-----------------------------

Defining the satellite's attitude requires the satellite's initial and final orientations in a common frame: the
:ref:`LVLH frame <Local-Vertical-Local-Horizontal (LVLH)>`. With this knowledge, the angles required to rotate from the
initial to final orientation can be computed; these angles are the roll, pitch, and yaw angles.

Assuming an Earth observation satellite, the initial orientation is defined as having the satellite's camera pointing
in the direction of the lvlh's positive z-axis. The final orientation is the orientation in the lvlh frame
where the camera is pointing in the direction a target of interest (i.e. when the camera is pointing towards a ground
location). The attitude angles are then determined and adhere to the following definitions:

#. Roll is the angle about the lvlh's x-axis,
#. Pitch is the angle about the lvlh's y-axis, and
#. Yaw is the angle about the lvlh's z-axis.

Rotating a vector from the satellites initial orientation to the final orientation is performed by applying a Euler-213
(or pitch-roll-yaw) sequence to the vector.

Satellite Look-Angle Definition
-------------------------------

The look-angle (or off-nadir angle) is the angle between the vector pointing from the satellite to the target and the
satellite's nadir vector. A zero look-angle indicates the satellite is directly overhead the target. The look-angle is
always positive by definition.

Physically, the look-angle is a measure of a ground target's position relative to the satellite's nadir. This can be a
useful metric for determining the satellite's image quality. If a large look-angle exists for a ground target, then
imaging of that location would induce significant skew in the resultant data.
