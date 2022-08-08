Coordinates
===========

.. contents:: Contents
   :depth: 1
   :local:

Coordinates Overview
--------------------

The coordinate module provides the base classes and functions used for coordinate definitions and conversions. This
includes the ability to define ground locations and different coordinate systems, in addition to defining the framework
for conversions between these systems.

Specifying Ground Locations
---------------------------

One of the primary objectives of Celest is to generate and schedule encounter times between a satellite and a ground
location; this requires the definition of a ground location which is accomplished using the :class:`GroundLocation`
class.

The :class:`GroundLocation` class can be imported via the following:

.. code-block:: python

   from celest.coordinates import GroundLocation

.. autoclass:: celest.coordinates.ground_location.GroundLocation
   :inherited-members:
   :member-order: bysource

Different Coordinate Frames
---------------------------

Celest provides the ability to deal with different coordinate frames and their conversions. The supported frames are:

#. :ref:`Roll-Pitch-Yaw Attitude Frame<Roll-Pitch-Yaw Attitude Frame>`,
#. :ref:`Azimuth-Elevation (Horizontal frame)<Azimuth-Elevation (Horizontal Frame)>`,
#. :ref:`Geocentric Celestial Reference System (GCRS)<Geocentric Celestial Reference System (GCRS)>`,
#. :ref:`International Terrestrial Reference System (ITRS)<International Terrestrial Reference System (ITRS)>`,
#. :ref:`Local-Vertical-Local-Horizontal (LVLH)<Local-Vertical-Local-Horizontal (LVLH)>`, and
#. :ref:`World Geodetic System 84 (WGS84)<World Geodetic System 84 (WGS84)>`.

Roll-Pitch-Yaw Attitude Frame
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Roll-Pitch-Yaw attitude frame is a 3-dimensional frame that defines the satellite orientation. The roll, pitch, and
yaw angles are the angles required to rotate the satellite from the
:ref:`lvlh frame<Local-Vertical-Local-Horizontal (LVLH)>` to an orientation where the satellite's nadir intersects the
ground location. The roll of the attitude frame is to define satellite orientation for ground target tracking during
ground imaging encounters (for a nadir facing camera).

The attitude frame is calculated from the :meth:`Satellite.attitude` method and requires both the satellite position and
velocity.

The :class:`Attitude` class can be imported via the following:

.. code-block:: python

   from celest.coordinates import Attitude

.. autoclass:: celest.coordinates.frames.attitude.Attitude
   :inherited-members:
   :noindex:

Azimuth-Elevation (Horizontal Frame)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Azimuth-Elevation frame (also known as the Horizontal or Altitude-Azimuth frame) is an observer centric frame that
uses the observers local horizon to define to angular measures: azimuth and elevation. Azimuth is the angle in the plane
of the observer's local horizon between the desired target and North. Azimuth is measured as clockwise positive and
lies in the range :math:`[0, 360)`. Elevation is the angle of the desired target above and perpendicular to the
observer's local horizon. Elevation is measured as positive above the horizontal plane and lies in the range
:math:`[-90, 90]`. This frame is used extensively in the window generation workflow.

Since the Azimuth-Elevation frame is dependent on the observer's location, conversions from frames such as the GCRS or
ITRS frame requires passing in an observer's location. For more information on initializing a ground location, refer to
the :ref:`Specifying Ground Locations` section.

The :class:`AzEl` class can be imported via the following:

.. code-block:: python

   from celest.coordinates import AzEl


.. autoclass:: celest.coordinates.frames.azel.AzEl
   :inherited-members:
   :noindex:

Geocentric Celestial Reference System (GCRS)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :class:`GCRS` class can be imported via the following:

.. code-block:: python

   from celest.coordinates import GCRS

.. autoclass:: celest.coordinates.frames.gcrs.GCRS
   :inherited-members:
   :noindex:

International Terrestrial Reference System (ITRS)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :class:`ITRS` class can be imported via the following:

.. code-block:: python

   from celest.coordinates import ITRS

.. autoclass:: celest.coordinates.frames.itrs.ITRS
   :inherited-members:
   :noindex:

Local-Vertical-Local-Horizontal (LVLH)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Local-Vertical-Local-Horizontal (LVLH) frame (also known as the Hill frame) is a 3-dimensional frame oriented with
respect to the satellite's orbit. The origin is located at the satellite's center of mass. The z-axis is aligned with
the satellites geocentric radius vector and is positive towards the Earth's surface. The x-axis lies in the orbital
plane and is perpendicular to the z-axis; the axis is positive in the direction of the satellite's motion. Lastly, the
y-axis is perpendicular to the orbital plane to complete the right handed orthogonal set.

The LVLH frame is used as a basis for satellite attitude determination. That is, the satellite's attitude is defined as
the roll, pitch, and yaw angles required to rotate the satellite from the LVLH frame to it's given orientation.

The :class:`LVLH` class can be imported via the following:

.. code-block:: python

   from celest.coordinates import LVLH

.. autoclass:: celest.coordinates.frames.lvlh.LVLH
   :inherited-members:
   :noindex:

World Geodetic System 84 (WGS84)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The World Geodetic System 84 (WGS84) is a reference ellipsoid used to model the Earth as an ellipsoid rather than a
sphere. It is defined by three measures: latitude, longitude, and height above the reference ellipsoid. The WGS84 frame
is used primarily in defining ground locations such as those used in satellite-to-ground encounters.

The :class:`WGS84` class can be imported via the following:

.. code-block:: python

   from celest.coordinates import WGS84

.. autoclass:: celest.coordinates.frames.wgs84.WGS84
   :inherited-members:
   :noindex:

Converting Between Coordinate Frames
------------------------------------

Conversions between the different coordinate frames are accessible through the :class:`Coordinate` class. However, there
exist constraints on various conversions.

#. Since the :ref:`Azimuth-Elevation frame<Azimuth-Elevation (Horizontal Frame)>` is observer centric, conversions to this frane require a :ref:`defined ground location<Specifying Ground Locations>` to be passed in as a parameter to the conversion method.
#. Conversions from the :ref:`Azimuth-Elevation frame<Azimuth-Elevation (Horizontal Frame)>` are not supported.
#. Conversions of velocity data to the :ref:`WGS84 frame<World Geodetic System 84 (WGS84)>` should be avoided to to poor physical significance.

For more information on using the coordinate conversions, refer to the tutorial on converting between frames found
:ref:`here <Position Conversions>`.

.. note::
   Velocity conversions can be handled much in the same way as position conversions. However, due to physical
   insignificance, velocity coordinates should not be converted into the :ref:`WGS84 frame<World Geodetic System 84 (WGS84)>`.

The :class:`Coordinate` class can be imported via the following:

.. code-block:: python

   from celest.coordinates import Coordinate

.. autoclass:: celest.coordinates.coordinate.Coordinate
   :noindex:
