Materials description
=====================

Default materials
-----------------

By default, Danton models the Earth with only three materials: :python:`"Rock"`,
:python:`"Water"` and :python:`"Air"`. These materials are defined within
a Materials Definition File (MDF) in `TOML`_ format, as outlined below.

.. literalinclude:: materials/default.toml
   :language: toml

.. topic:: Geometry composition

    The :ref:`interior <geometry:Earth structure>` layers of the Earth, in
    addition to the :ref:`topography <geometry:Topography>` layer, are composed
    of :python:`"Rock"` with varying bulk density. The `PREM`_ ocean, or oceanic
    regions, are constituted by :python:`"Water"`. The :ref:`atmospheric
    <geometry:Atmosphere>` layers are constituted by :python:`"Air"` with
    varying density.


Custom materials
----------------

It is possible to use custom materials by providing one's own MDF, with the
section titles corresponding to the material names (capitalised, by convention).
It is required that each material section indicates the material
:python:`density` (in kg/m\ :sup:`3`) and its :python:`composition`.
Additionally, the material mean excitation energy (`MEE`_) may be provided using
the :python:`I` attribute (expressed in GeV), if desired.

In the event that the :python:`"Rock"` :python:`"Water"` or :python:`"Air"`
material is not defined, Danton will add the corresponding default entry to the
MDF. For example, the following file extends the default materials by defining
two additional types of rocks, :python:`"Limestone"` and :python:`"Sandstone"`.

.. literalinclude:: materials/rocks.toml
   :language: toml

.. note::

   It is only possible to assign a different material to the :ref:`topography
   <geometry:Topography>` layer. However, the :python:`"Rock"`,
   :python:`"Water"` and :python:`"Air"` materials can be modified by redefining
   the corresponding MDF entries.

.. topic:: Composition

   The material :python:`composition` can be expressed in the form of a
   molecular formula (e.g. :python:`"H2O"`) or by providing a mass composition
   of atomic elements, molecules, or other materials (see the preceding
   examples). It should be noted that in the event of a reference to other
   materials, these must have been defined *before* within the MDF.


.. ============================================================================
.. 
.. URL links.
.. 
.. ============================================================================

.. _MEE: https://en.wikipedia.org/wiki/Bethe_formula#The_mean_excitation_energy
.. _PREM: https://en.wikipedia.org/wiki/Preliminary_reference_Earth_model
.. _TOML: https://toml.io/en/
