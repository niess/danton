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

Danton assumes that primary neutrinos always originate from the top of the
atmosphere. Therefore, in the forward Monte Carlo case, as a conveniency, the
input position actually defines a target point for the incoming neutrino, not
necessarily its injection point. Before running the Monte Carlo, the neutrino is
first relocated to the top of the atmosphere by extrapolating its track
backwards from the target point.

.. note::

   The relocation of primary neutrinos is purely geometrical. Therefore,
   simulated neutrinos might not intersect their target point, if the
   longitudinal approximation is disabled.


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
   likelihood (i.e. :math:`\omega = 1 / \text{pdf}(\text{S})`, for state
   :math:`\text{S}`). This is consistent with the backward Monte Carlo workflow.
   Set to :python:`False` the :python:`weight` argument of the
   :py:func:`particles <danton.Simulation.particles>` method in order to disable
   weighting. Note that the weighting can also be enabled or disabled at the
   level of individual :py:class:`generator <danton.ParticlesGenerator>`
   methods.

.. topic:: Sample size

   Depending on its configuration, the :py:class:`generator
   <danton.ParticlesGenerator>` might rely on rejection sampling methods,
   resulting in the Monte Carlo samples size being actually larger than the
   number of returned particles (N). In such cases, in addition to the selected
   particles, the :py:func:`generate <danton.ParticlesGenerator.generate>`
   method also returns the actual samples size as second argument.


Local box
---------

The :py:class:`Box <danton.Box>` object let us define a region of interest using
a bounding-box. Besides, it also allows us to use :ref:`local coordinates
<coordinates:Local coordinates>` instead of :ref:`geographic
<coordinates:Geographic coordinates>` ones. The :py:class:`Box <danton.Box>`
region can then be set as a :py:func:`target <danton.ParticlesGenerator.target>`
for a :py:class:`generator <danton.ParticlesGenerator>`, typically in the case
of a forward Monte Carlo. In the backward case, one might instead configure the
:py:class:`generator <danton.ParticlesGenerator>` to sample decays
:py:func:`inside <danton.ParticlesGenerator.inside>` the box.


Monte Carlo history
-------------------


.. ============================================================================
.. 
.. URL links.
.. 
.. ============================================================================

.. _builder: https://en.wikipedia.org/wiki/Builder_pattern
