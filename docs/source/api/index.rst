API Reference
=============

This section contains the Python API documentation for ``pyfraglib``.

Core Classes
------------

.. autosummary::
   :toctree: generated
   :nosignatures:

   pyfraglib.Fragment
   pyfraglib.FragmentList
   pyfraglib.FragFile

Module List
-----------

.. autosummary::
   :toctree: generated
   :nosignatures:

   pyfraglib.core
   pyfraglib.fragment
   pyfraglib.fragfile
   pyfraglib.math
   pyfraglib.stats
   pyfraglib.lengths
   pyfraglib.scores
   pyfraglib.simulator

Core Simulation Classes
-----------------------

The simulation code lives in its own module, but the most important classes get
re-exported into the ``pyfraglib`` namespace:

.. autosummary::
   :toctree: generated
   :nosignatures:

   pyfraglib.FragmentSimulator
   pyfraglib.TissueMixtureSimulator
   pyfraglib.NucleaseProfile

