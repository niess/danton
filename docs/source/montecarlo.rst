Monte Carlo workflow
====================

In order to simulate the :math:`\nu_\tau-\tau` transport through a section of
the :doc:`Earth <geometry>`, Danton employs a Monte Carlo technique. In the
traditional approach, primary neutrinos are injected at the top of the
:ref:`atmosphere <geometry:Atmosphere>` and the resulting tau decays are
collected.

In addition to the traditional forward method, Danton also allows for the
backward sampling of tau decays. In practice, tau leptons are injected at decay
points of interest, and potential primary neutrinos are collected at the top of
the :ref:`atmosphere <geometry:Atmosphere>`. In the backward case, the candidate
primary neutrinos have associated Monte Carlo weights indicating their
respective likelihoods.

.. topic:: Backward or Forward

   The default setting for Danton is a backward Monte Carlo simulation. To
   perform a conventional forward simulation, one must set the :py:attr:`mode
   <danton.Simulation.mode>` attribute of the :py:class:`simulation
   <danton.Simulation>` to :python:`"forward"`.

   Backward sampling can enhance CPU efficiency by orders of magnitude in
   certain use cases (see, e.g., [Niess18b]_). However, this approach is more
   more intricate, given that nature progresses forwards.

.. topic:: Decays or fluxes

   The default configuration of Danton is to sample decay densities. However, it
   is also possible to sample fluxes of particles instead, should this be
   required. This can be achieved by setting the :py:attr:`tau_decays
   <danton.Simulation.tau_decays>` attribute of the :py:class:`simulation
   <danton.Simulation>` to :python:`False`.

.. topic:: Longitudinal approximation

   By default, Danton performs a comprehensive three-dimensional simulation of
   particle transport. However, at ultra-high energies, the deflection angles
   w.r.t. the initial direction of the neutrinos are typically minimal, and thus
   often disregarded. This approximation is activated by setting the
   :py:attr:`longitudinal <danton.Simulation.longitudinal>` attribute of the
   :py:class:`simulation <danton.Simulation>` to :python:`True`.


Target point
------------

Danton's approach is to assume that primary neutrinos always originate from the
top of the atmosphere. In the forward Monte Carlo case, for the sake of
convenience, the input position is considered to be a target point for the
incoming neutrino, in lieu of its actual injection point. Prior to running the
Monte Carlo, the neutrino is repositioned to the top of the atmosphere by
extrapolating its track backwards from the target point.

.. note::

   The relocation of primary neutrinos is purely geometrical. Thus, the
   simulated neutrinos may not intersect their target point, for example, if the
   longitudinal approximation is disabled.


Point estimate
--------------

A salient property of backward Monte Carlo is its ability to perform a point
estimate for a final state of interest. To illustrate, Danton may compute the
density of tau decays at a specific location on Earth, with a given direction
and tau energy at decay. This is obtained through a Monte Carlo weighted
average of the fluxes of backward sampled neutrinos. For instance,

.. literalinclude:: montecarlo/point-estimate.py
   :language: python
   :lines: 9-16


Monte Carlo integration
-----------------------

In a typical scenario, the objective is to collect particles within a specified
region of interest, e.g., one that matches the field of view of a detector. This
requires performing a Monte Carlo integration over initial (final) states in the
forward (backward) case. To elaborate, one would generate incoming neutrinos
over the top of the atmosphere, in the forward case (or tau decay vertices over
the region of interest, in the backward case).

In order to facilitate these generation procedures, Danton provides a
:py:class:`ParticlesGenerator <danton.ParticlesGenerator>` object, which is
instantiated using the :py:func:`particles <danton.Simulation.particles>` method
of the :py:class:`Simulation <danton.Simulation>` class. The particles generator
is configured using a `builder`_ pattern, and then triggered with the
:py:func:`generate <danton.ParticlesGenerator.generate>` method. For instance,
the following

.. literalinclude:: montecarlo/generator.py
   :language: python
   :lines: 6-9

would result in the generation of N particles with energies distributed
according to a :math:`1 / E^2` power-law (between 1 |nbsp| PeV and 1 |nbsp|
EeV), and across a solid angle encompassing elevation values between -1 and 1
|nbsp| deg.

.. topic:: Generation weight

   The default setting is to weight particles by the inverse of their generation
   likelihood (:math:`\omega = 1 / \text{pdf}(\text{S})`, for a Monte Carlo
   state :math:`\text{S}`). The :python:`weight` argument of the
   :py:func:`particles <danton.Simulation.particles>` method can be set to
   :python:`False` to disable weighting. It should be noted that weighting can
   also be enabled or disabled individually at the level of each
   :py:class:`generator <danton.ParticlesGenerator>` method.

.. topic:: Sample size

   The :py:class:`generator <danton.ParticlesGenerator>` may employ rejection
   sampling techniques (depending on its configuration), resulting in a Monte
   Carlo sample that is larger than the number of returned particles (N). In
   this case, along with the selected particles, the :py:class:`generate
   <danton.ParticlesGenerator.generate>` method also returns the sample size as
   a second argument.

.. topic:: Box region

   The :py:class:`Box <danton.Box>` object allows one to specify a region of
   interest by means of a bounding box. This region may subsequently be
   :py:func:`targeted <danton.ParticlesGenerator.target>` by a
   :py:class:`generator <danton.ParticlesGenerator>`, e.g. when employing
   forward Monte Carlo. Conversely, in the backward case, the
   :py:class:`generator <danton.ParticlesGenerator>` can be configured to
   sample decay vertices :py:func:`inside <danton.ParticlesGenerator.inside>`
   the box.


Monte Carlo history
-------------------

Danton simulations are pseudo-random, relying on a Pseudo-Random Number
Generator (`PRNG`_) as the source of randomness (a Permuted Congruential
Generator (`PCG`_) namely). The generator is, by default, seeded from the
operating system's `entropy`_. The :py:attr:`seed <danton.Random.seed>` value,
in conjunction with the `PRNG`_ stream :py:attr:`index <danton.Random.index>`,
determines the subsequent fate of a Monte Carlo particle. This can be exploited
to re-simulate a specific event, as Danton systematically logs the `PRNG`_
stream :py:attr:`index <danton.Random.index>` in addition to the Monte Carlo
particles data.

A typical use case is to re-simulate selected events by enabling a detailed log
of Monte Carlo steps (by setting the :py:attr:`record_steps
<danton.Simulation.record_steps>` attribute of the :py:class:`simulation
<danton.Simulation>` to :python:`True`).

.. note::

   It is inadvisable to systematically record Monte Carlo steps for large
   samples due to the resulting memory overhead.


.. ============================================================================
.. 
.. URL links.
.. 
.. ============================================================================

.. _builder: https://en.wikipedia.org/wiki/Builder_pattern
.. _entropy: https://en.wikipedia.org/wiki/Entropy_(computing)
.. _PRNG: https://en.wikipedia.org/wiki/Pseudorandom_number_generator
.. _PCG: https://en.wikipedia.org/wiki/Permuted_congruential_generator
