Physics description
===================

Danton implements an original backward Monte Carlo technique (see [NBCL18]_,
[NM18]_ and [Nie23]_), which requires special (backward-compliant) Monte
Carlo engines, as described below.

.. topic:: Validation

   The Danton physics implementation has been validated, e.g. by comparisons
   with `NuTauSim`_ (see [NM18]_) as well as with `NuPropEarth`_ (see [AAA+22]_,
   sec. |nbsp| 6.6, fig. |nbsp| 63).


Neutrino interactions
---------------------

The simulation of high-energy neutrino interactions is carried out with `Ent`_.
The dominant process is Deep Inelastic Scattering (`DIS`_) on nuclei partons,
through the exchange of a virtual W or Z boson. Additionally, collisions with
atomic electrons are also considered, including the Glashow Resonant (`GR`_)
process.

.. topic:: DIS collisions

   DIS collisions are randomised from the Doubly Differential Cross-Section
   (DDCS) in Bjorken-:math:`x` and :math:`Q^2`, using Leading Order (LO)
   expressions with configurable Parton Distribution Functions (PDFs). However,
   the total DIS cross-sections are rescaled to more detailed computations (e.g.
   the [CMS11]_ or [BGR18]_ cross-sections).


Tau interactions
----------------

Tau interactions are simulated with `Pumas`_, a dedicated heavy leptons
transport engine, initially developed for `muography`_ applications. `Pumas`_ has
been extensively validated, by comparisons to other Monte Carlo transport
engines, including `Geant4`_ and `PROPOSAL`_ (see [Nie22]_).

High-energy tau interactions are typically divided into two categories:
ionisation (resulting from inelastic collisions with atomic electrons) and
radiative processes (dominant above ~10 |nbsp| TeV). Radiative processes
encompasses `bremsstrahlung`_, :math:`e^+e^-` `pair-production`_ and
photonuclear interactions (i.e. `DIS`_ through the exchange of a virtual
photon). Additionally, `Pumas`_ can also account for the multiple scattering off
nuclei through Coulomb collisions.

.. topic:: Cross-sections

   `Pumas`_ implements a variety of cross-section models for radiative processes
   (see e.g. the :py:class:`Physics <danton.Physics>` class). By default, the
   `PROPOSAL`_ cross-sections  are utilised for the `bremsstrahlung`_ and
   `pair-production`_ processes (see e.g. [SSR19]_). Ionisation and multiple
   scattering processes are modelled following [Sal13]_, which is comparable to
   the `Penelope`_ implementation. For further information on `Pumas`_, refer to
   [Nie22]_.


Tau decays
----------

Tau decays are simulated with `Alouette`_, which is a thin wrapper around
`TAUOLA`_ (C++ release), allowing for backward decays. This is essentially
achieved by a relativistic change of frame (see [Nie23]_), starting from
`TAUOLA`_-generated Centre-of-Mass (CM) decays.

.. topic:: Polarisation

   The angular distribution of the decays products (in the CM) is dependent on
   the tau polarisation at decay. Danton's assumption is that tau leptons are
   100% anti-polarised at production, and that this polarisation is retained
   until the decay.


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
.. _TAUOLA: https://tauolapp.web.cern.ch/tauolapp
