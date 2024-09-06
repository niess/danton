Physics description
===================

Danton implements an original backward Monte-Carlo technique (see [Niess18a]_,
[Niess18b]_ and [Niess23]_), requiring dedicated (backward-compliant) Monte
Carlo engines, as described below.

.. topic:: Validation

   Danton physics implementation has been validated, e.g. by comparisons to
   `NuTauSim`_ (see [Niess18b]_) as well as to `NuPropEarth`_ (see [Abraham22]_,
   sec. |nbsp| 6.6, fig. |nbsp| 63).


Neutrino interactions
---------------------

High energy neutrino interactions are simulated with `Ent`_. The dominant
process is Deep Inelastic Scattering (`DIS`_) on nuclei partons, through the
exchange of a virtual W or Z boson. Collisions with atomic electrons are also
considered, including the Glashow Resonant (`GR`_) process.

.. topic:: DIS collisions

   DIS processes are randomised from the Doubly Differential Cross-Section
   (DDCS) in Bjorken-:math:`x` and :math:`Q^2`, using Leading Order (LO)
   expressions with configurable Parton Distribution Functions (PDFs). The total
   DIS cross-sections are however rescaled to more detailed computations (e.g.
   the [CMS11]_ or [BGR18]_ cross-sections).


Tau interactions
----------------

Tau interactions are simulated with `Pumas`_, a dedicated heavy leptons
transport engine, initially developed for `muography`_ applications. `Pumas`_ has
been extensivelly validated, by comparisons to other Monte Carlo transport
engines, including `Geant4`_ and `PROPOSAL`_ (see [Niess22]_).

High energy tau interactions are usually separated into two components,
ionisation (by inelastic collisions with atomic electrons) and radiative
processes (dominant above ~10 |nbsp| TeV). Radiative processes include
`bremsstrahlung`_, :math:`e^+e^-` `pair-production`_ and photonuclear
interactions (i.e. `DIS`_ through the exchange of a virtual photon). In
addition, `Pumas`_ can also account for the multiple scattering off nuclei
through Coulomb collisions.

.. topic:: Cross-sections

   `Pumas`_ implements various models of cross-sections for the radiative
   processes (see e.g. the :py:class:`Physics <danton.Physics>` class). By
   default, the `PROPOSAL`_ cross-sections (Sandrock, Soedingresko and Rhode
   (`SSR`_)) are used for the `bremsstrahlung`_ and `pair-production`_
   processes. Ionisation and multiple scattering processes are modelled
   according to [Salvat13]_, which is similar to `Penelope`_ implementation.
   See [Niess22]_ for more details on `Pumas`_.


Tau decays
----------

Tau decays are simulated with `Alouette`_, which is a thin wrapper around
`TAUOLA`_ (C++ release), allowing for backward decays. This is essentially
achieved by a relativistic change of frame (see [Niess23]_), starting from
`TAUOLA`_ generated Centre-of-Mass (CM) decays.

.. topic:: Polarisation

   The angular distribution of decays products (in the CM) depends on the tau
   polarisation at decay. Danton assumes that tau leptons are 100%
   anti-polarised at production (through DIS processes, essentially), and that
   they retain their polarisation until their decay.


.. ============================================================================
.. 
.. URL links.
.. 
.. ============================================================================

.. _Alouette: https://github.com/niess/alouette
.. _Bremsstrahlung: https://en.wikipedia.org/wiki/Bremsstrahlung
.. _DIS: https://en.wikipedia.org/wiki/Deep_inelastic_scattering
.. _Ent: https://github.com/niess/ent
.. _Geant4: https://geant4.web.cern.ch
.. _GR: https://en.wikipedia.org/wiki/Glashow_resonance
.. _muography: https://en.wikipedia.org/wiki/Muon_tomography
.. _NuPropEarth: https://github.com/pochoarus/NuPropEarth
.. _NuTauSim: https://github.com/harmscho/NuTauSim
.. _Pair-Production: https://en.wikipedia.org/wiki/Pair_production
.. _Pumas: https://github.com/niess/pumas
.. _Penelope: https://www.oecd-nea.org/upload/docs/application/pdf/2020-10/penelope-2018__a_code_system_for_monte_carlo_simulation_of_electron_and_photon_transport.pdf
.. _PROPOSAL: https://github.com/tudo-astroparticlephysics/PROPOSAL
.. _SSR: https://arxiv.org/abs/1910.07050
.. _TAUOLA: https://tauolapp.web.cern.ch/tauolapp
