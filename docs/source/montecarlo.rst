Monte Carlo workflow
====================

In order to simulate the transport of :math:`\nu_\tau-\tau` leptons through a
section of the Earth, Danton employs a Monte Carlo technique. In the traditional
approach, primary neutrinos are injected at the top of the :ref:`atmosphere
<geometry:Atmosphere>` and the resulting tau decays are collected.

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

   Backward sampling can enhance CPU efficiency in certain Monte Carlo
   simulations (see, e.g., [Niess18b]_). However, this approach entails a
   more intricate workflow, given that nature progresses forwards.

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

Danton assumes that primary neutrinos always originate from the top of the
atmosphere. Therefore, in the forward Monte Carlo case, as a conveniency, the
input position defines a target point for the incoming neutrino, not its actual
injection point (though, the latter might be set as target). Before running the
Monte Carlo, the neutrino is first relocated to the top of the atmosphere by
extrapolating its track backwards from the target point.

.. note::

   The relocation of primary neutrinos is purely geometrical. The simulated
   neutrinos might not intersect their target point, e.g. if the longitudinal
   approximation is disabled.


Point estimate
--------------

A unique feature of backward Monte Carlo is that it let us perform a point
estimate for a final state of interest. For instance, Danton can compute the
density of tau decays at a given location on the Earth (and for a given
direction and tau energy at decay). This is done simply by performing the Monte
Carlo weighted average of the fluxes of backward sampled neutrinos. For example,

.. literalinclude:: montecarlo/point-estimate.py
   :language: python
   :lines: 9-16


Monte Carlo integration
-----------------------

Usually, one is interested in collecting particles over a region of interest,
typically matching the field of view of a detector, rather than at a specific
location. This requires performing a Monte Carlo integration over initial
(final) states in the forward (backward) case. For instance, one would generate
incoming neutrinos over the top of the :ref:`atmosphere <geometry:Atmosphere>`
(in the forward case), or tau decay points over a region of interest (in the
backward case).

In order to facilitate these generation procedures, Danton provides a
:py:class:`ParticlesGenerator <danton.ParticlesGenerator>` object, which is
instanciated using the :py:func:`particles <danton.Simulation.particles>` method
of the :py:class:`simulation <danton.Simulation>`. The particles generator is
configured using a `builder`_ pattern, and then triggered with the
:py:func:`generate <danton.ParticlesGenerator.generate>` method. For instance,
the following

.. literalinclude:: montecarlo/generator.py
   :language: python
   :lines: 6-9

would generate N particles using a :math:`1 / E^2` power-law distribution for
the energy (between 1 |nbsp| PeV and 1 |nbsp| EeV), and over a solid angle
spanning -1 to 1 |nbsp| deg of elevation values.

.. topic:: Generation weight

   By default, the generated particles are weighted by the inverse of their
   likelihood (i.e. :math:`\omega = 1 / \text{pdf}(\text{S})`, for a Monte Carlo
   state :math:`\text{S}`). This is consistent with the backward Monte Carlo
   workflow. Set to :python:`False` the :python:`weight` argument of the
   :py:func:`particles <danton.Simulation.particles>` method in order to disable
   weighting. Note that the weighting can also be enabled or disabled at the
   level of individual :py:class:`generator <danton.ParticlesGenerator>`
   methods.

.. topic:: Sample size

   Depending on its configuration, the :py:class:`generator
   <danton.ParticlesGenerator>` might rely on rejection sampling methods,
   resulting in the Monte Carlo samples size being larger than the number of
   returned particles (N). In such cases, in addition to the selected particles,
   the :py:func:`generate <danton.ParticlesGenerator.generate>` method also
   returns the samples size as second argument.


.. topic:: Box region

   The :py:class:`Box <danton.Box>` object let us define a region of interest
   using a bounding-box. This region can be :py:func:`targeted
   <danton.ParticlesGenerator.target>` by a :py:class:`generator
   <danton.ParticlesGenerator>`, typically in the case of a forward Monte Carlo.
   In the backward case, one might instead require the :py:class:`generator
   <danton.ParticlesGenerator>` to sample decays :py:func:`inside
   <danton.ParticlesGenerator.inside>` the box volume.


Monte Carlo history
-------------------


.. ============================================================================
.. 
.. URL links.
.. 
.. ============================================================================

.. _builder: https://en.wikipedia.org/wiki/Builder_pattern
