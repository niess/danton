Materials description
=====================

Default materials
-----------------

By default, Danton models the Earth with only three materials: :python:`"Rock"`,
:python:`"Water"` and :python:`"Air"`. These default materials are defined from
a Materials Definition File (MDF) in `TOML`_ format, as following.

.. literalinclude:: materials/default.toml
   :language: toml

.. topic:: Geometry composition

    The Earth :ref:`interior <geometry:Earth structure>` layers, as well as the
    :ref:`topography <geometry:Topography>` layer, consist of :python:`"Rock"`
    (with a varying bulk density). The `PREM`_ ocean, or oceanic regions,
    consist of :python:`"Water"`. The :ref:`atmosphere <geometry:Atmosphere>`
    layers are made of :python:`"Air"` (with a varying density).


Custom materials
----------------

Custom materials sets can be used by providing your own MDF, where the sections
title correspond to the materials names (capitalised, by convention). Each
material section must indicate the material :python:`density` (in kg/m\
:sup:`3`) and its :python:`composition`. Optionally, the material Mean
Excitation Energy (`MEE`_) can also be provided, using the :python:`I` attribute
(expressed in GeV).

If not defined, Danton will add the default :python:`"Rock"`, :python:`"Water"`
and :python:`"Air"` entries to the material set. For instance, the following
file extends the available materials by defining two additional types of rocks,
:python:`"Limestone"` and :python:`"Sandstone"`.

.. literalinclude:: materials/rocks.toml
   :language: toml

.. note::

   Only the :ref:`topography <geometry:Topography>` layer can be set to a
   different material than its default, i.e. :python:`"Rock"`. However, the
   :python:`"Rock"`, :python:`"Water"` and :python:`"Air"` materials can be
   modified by redefining the corresponding MDF entries.

.. topic:: Composition

   The material :python:`composition` can be expressed as a molecular formula
   (e.g. :python:`"H2O"`), or by giving a mass composition of atomic elements,
   molecules, or other materials (see the previous examples). Note that when
   refering to other materials, those must have been defined above in the MDF.


.. ============================================================================
.. 
.. URL links.
.. 
.. ============================================================================

.. _MEE: https://en.wikipedia.org/wiki/Bethe_formula#The_mean_excitation_energy
.. _PREM: https://en.wikipedia.org/wiki/Preliminary_reference_Earth_model
.. _TOML: https://toml.io/en/
