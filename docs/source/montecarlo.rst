Monte Carlo workflow
====================

Danton simulates, by Monte Carlo, the coupled transport of :math:`\nu_\tau-\tau`
leptons through a section of the Earth. Traditionaly, this would be done by
injecting primary neutrinos (on top of the :ref:`atmosphere
<geometry:Atmosphere>` layers), and by collecting the resulting tau decays
(within any :ref:`atmosphere <geometry:Atmosphere>` layer).

In addition to the traditional (forward) method, Danton also allows for the
backward sampling of tau decays. In practice, one injects tau leptons at decay
points of interest, and one collects potential primary neutrinos at the top of
the :ref:`atmosphere <geometry:Atmosphere>` layers. In the backward case, the
candidate primary neutrinos have associated Monte Carlo weights indicating their
respective likelihoods.

.. topic:: Backward or Forward

   By default, Danton performs a backward Monte Carlo simulation. In order to
   instead perform a conventional forward simulation, set to :python:`"forward"`
   the :py:attr:`mode <danton.Simulation.mode>` attribute of the
   :py:class:`simulation <danton.Simulation>`.

   Backward sampling can greatly increase the CPU efficiency of a Monte Carlo
   simulation under some circumstances (see e.g. [Niess18b]_). However, this
   usually implies a (slightly) more complex workflow, since nature actually
   goes forward.

.. topic:: Decays or fluxes

   By default, Danton is configured to sample decay densities. However, Danton
   can also sample fluxes of particles instead, if desired. This is done by
   setting to :python:`False` the :py:attr:`tau_decays
   <danton.Simulation.tau_decays>` attribute of the :py:class:`simulation
   <danton.Simulation>`.

.. topic:: Longitudinal approximation

   By default, Danton performs a full 3D simulation of particles transport.
   However, at ultra-high energies deflection angles w.r.t. the
   primary neutrino direction are usually small, and thus frequently neglected.
   This `longitudinal` approximation is enabled by setting to :python:`True` the
   :py:attr:`longitudinal <danton.Simulation.longitudinal>` attribute of the
   :py:class:`simulation <danton.Simulation>`.


Target point
------------


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
:py:class:`ParticlesGenerator <danton.ParticlesGenerator>` object, which can be
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

   By default, ...

.. topic:: Sample size

   By default, ...


Local box
---------

The :py:class:`Box <danton.Box>` object let us define a region of interest using
a bounding-box. Besides, it also allows us to use :ref:`local coordinates
<coordinates:Local coordinates>` instead of :ref:`geographic
<coordinates:Geographic coordinates>` ones. This :py:class:`Box <danton.Box>`
region can then be set as a target for a :py:class:`ParticlesGenerator
<danton.ParticlesGenerator>`, in the forward case typically. In the backward
case, the box interior can be used to specify potential decay points for the
generator.



Monte Carlo history
-------------------


.. ============================================================================
.. 
.. URL links.
.. 
.. ============================================================================

.. _builder: https://en.wikipedia.org/wiki/Builder_pattern
