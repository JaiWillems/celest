Schedule
========

.. contents:: Contents
   :depth: 1
   :local:

Schedule Overview
-----------------

The primary usage of Celest is for mission planning of Earth observation satellites. A key part of this problem is the
scheduling of observations. The following sections describe the concept of a scheduling request, the scheduling
framework apart of celest, and a discussion on the scheduling algorithm.

Concept of Requests
-------------------

A scheduling request is a formal description of a desired observation that outlines the constraints that must be met
for the observation to be viable. The parameters that govern a request include the following:

#. The ground location of the satellite-ground encounter.
#. The deadline for the request to be fulfilled by.
#. The duration required for the observation.
#. The priority of the request compared to other requests.
#. The minimum quality of the imaging data.
#. A specific look-angle for the observation.
#. The lighting condition for the request.

Schedule Class
--------------

The purpose of the scheduler is take a series of scheduling requests and find a near-optimal set of viable observations
that fulfill the highest number of requests. In generating the schedule, the scheduler will ensure all request
constraints are met and that there exist no conflicts between scheduled observations. The :class:`Scheduler` class
documentation can be found below and an example usage can be found in the :ref:`Scheduling Workflow` tutorial.

The :class:`Scheduler` class can be imported via the following:

.. code-block:: python

   from celest.schedule import Scheduler

.. autoclass:: celest.schedule.scheduler.Scheduler
   :members:
   :noindex:

Scheduling Algorithm
--------------------

The scheduling algorithm used by the :class:`Scheduler` class is an adaptive large neighborhood search (ALNS)
metaheuristic with a simulated annealing acceptance criterion and is adapted from the work of Xiaolu Liu, et al. (2017).
[Liu+17]_ The following gives a brief overview to the algorithm which can be read in greater detail from the paper by
Liu et al. [Liu+17]_

.. [Liu+17] Xiaolu Liu et al. “An adaptive large neighborhood search metaheuristic for agile satellite scheduling with
   time-dependent transition time”. en. In: Computers Operations Research 86 (Oct. 2017), pp. 41–53. issn: 03050548.
   doi: 10.1016/j.cor.2017.04.006. url: https://linkinghub.elsevier.com/retrieve/pii/S0305054817300977.

Initial Solution
^^^^^^^^^^^^^^^^

The scheduling algorithm uses the sum of scheduled requires priorities as a measure of the quality of the solution.
That is, the goal of the algorithm is to determine the highest sum of priorities of all scheduled observations. As a
result, a good initial solution is greedy one. The :class:`Scheduler`'s ALNS algorithm sorts all requests in decreasing
order of their priorities and then in increasing order of start time; the scheduler then attempts to schedule all tasks
in order whilst adhering to the request constraints and takes the resulting solution to be the initial solution.

Note that if all requests are scheduled in the initial solution, the algorithm stops and returns the initial solution.

Scheduler Iteration
^^^^^^^^^^^^^^^^^^^

The general methodology of a neighborhood search algorithm is to create a new solution in the neighborhood of the
current solution by making small modifications through the use of destroy and repair operations. In each iteration,
the scheduler will apply destroy and repair operation to the current solution to change it slightly. The cost of the
resulting solution (the sum of the priorities of scheduled requests) is determined and compared to the cost of the
current solution. The new solution is accepted by the :ref:`simulated annealing acceptance criterion<Acceptance Criterion>`.
If the new solution is better than the best solution to date, the best solution is updated.

Acceptance Criterion
^^^^^^^^^^^^^^^^^^^^

A new solution is accepted by the simulated annealing acceptance criterion. If the cost of the new solution is better
than the cost of the current solution, the new solution is accepted. If the cost of the new solution is worse than the
cost of the current solution, the new solution is accepted with the following probability:

.. math:: \rho = \exp\left(\frac{100}{T}\frac{f(s)' - f(s)}{f(s)}\right)

where :math:`T` is the current temperature, :math:`f(s)` is the cost of the current solution and :math:`f(s)'` is the
cost of the new solution.

Allowing non-improving solutions to be occasionally accepted builds robustness against getting stuck in a local optimum.

Adaptive Strategy
^^^^^^^^^^^^^^^^^

The adaptive strategy of the algorithm is incorporated into the selection of the destroy and repair operations using
weights. Initially the weights are all set to be equal but as the algorithm progresses, the weights are adjusted based
on their impact in improving the cost of the solution. The weight for the :math:`i^{th}` destroy or repair function,
:math:`w_i` from a total of :math:`H` destroy or repair functions is adjusted using the following formula:

    .. math:: w_i = (1 - \lambda)w_i + \lambda\frac{\pi_i}{\sum_{j=1}^H\pi_j}

where :math:`\lambda` is the reactivity number and :math:`\pi_i` is the weight of the :math:`i^{th}` destroy operator.

The reactivity number governs the effect that previous iterations have on the current wight and can be passed into the
scheduler as a parameter. A reactivity number of 1 will cause the algorithm is disregard the success of the operator on
previous iterations. A zero reactivity factor will prevent the wights from being updated and will retain the initial
and equal weights.

The both the destroy and repair operators are chosen using a probability distribution based on the weights where the
probability of each operator is defined as:

.. math:: p_i = w_i/\sum_{j=1}^Hw_j
