Encounter
=========

.. contents:: Contents
   :depth: 1
   :local:

Encounter Overview
------------------

One of the primary goals of Celest is mission planning of encounter based tasks such as ground imaging or data transfer;
this requires knowledge of when the satellite is in a line of sight of the ground location. Windows are the definition
of viable encounter times and can take different forms. The following sections will look more into the concept of
windows, demonstrate how to determine them using Celest, and then show the structures to handle and interface with them.

More About Windows
------------------

An encounter is any position dependent interaction between a satellite and some ground location. One type of encounter
may be an imaging encounter where a satellite wishes to image a ground location. Another instance may be downlinking
said imaging data to a ground station. In either case, a line of sight between the satellite and ground location is
required. Windows are the way of defining when this line of sight is viable.

Windows in Celest come in two flavours: visible time windows and observation windows.

Visible Time Windows
^^^^^^^^^^^^^^^^^^^^

A visible time window defines the period of time when the satellite is in direct line of sight with the ground location
and is characterized by two values: the satellite rise and set times as seen from the ground location. Between these
times, the satellite is in view of the ground location and data transmission or ground imaging is theoretically
possible.

Observation Windows
^^^^^^^^^^^^^^^^^^^

Visible time windows can be imagined as generalized to any satellite-ground encounter and define the period of time
where a line of sight exists. However, the entire visible time window is more than sufficient for imaging encounters;
additionally, imaging encounters will provide other constraints such as image quality (influencing the maximum
look-angle) or a fixed look-angle. An observation window is catered to imaging tasks and is a subset of a visible time
window which meets the addition constraints of image quality and a fixed look-angle.

Observation windows are generated using the :class:`Scheduler` class.


Generating Visible Time Windows
-------------------------------

All possible visible time windows for a satellite-ground encounter can be determined using the :py:func:`generate_vtws`
function. For more information on generating visible time windows, see the
:ref:`window generation tutorial<Window Generation Workflow>`.

The :py:func:`generate_vtws` function can be imported via the following:

.. code-block:: python

   from celest.encounter import generate_vtws

.. autofunction:: celest.encounter.window_generator.generate_vtws
   :noindex:

The :class:`Lighting` enumeration is used to specify the lighting conditions for the visible time windows and can be
imported via the following:

.. code-block:: python

   from celest.encounter import Lighting

.. autoclass:: celest.encounter.window_generator.Lighting
   :noindex:

Window Data Structures
----------------------

Celest contains various data structures to hold individual windows and collections of windows. The
:class:`VisibleTimeWindow` class holds information regarding a single visible time window. Similarly, the
:class:`ObservationWindow` class holds information regarding a single observation window. The :class:`WindowCollection`
class is used to hold collections of either visible time windows or observation windows.

.. autoclass:: celest.encounter.window_handling.VisibleTimeWindow
   :noindex:

.. autoclass:: celest.encounter.window_handling.ObservationWindow
   :noindex:

.. autoclass:: celest.encounter.window_handling.WindowCollection
   :noindex:
