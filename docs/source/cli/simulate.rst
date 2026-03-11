Simulate Command
================

The ``simulate`` command generates synthetic cfDNA fragment data based on biological parameters and real genomic sequences.

Syntax
------

.. code-block:: bash

   pyfrag.py simulate -o <OUT_DIR> [OPTIONS]

Options
-------

.. code-block:: bash

   --config PATH          Configuration file (JSON format) (required)

Examples
--------

.. code-block:: bash

   pyfrag.py --out-dir synthetic/ simulate --config configs/simulation_example.json


Configuration File Format
--------------------------

Please take a look at the example configurations in ``configs/``. Most options should be self-explanatory.

Simulation Modes
----------------

Basic Mode
~~~~~~~~~~

- Single tissue type simulation
- Uniform fragment characteristics
- Customizable length distribution
- Simple nuclease profile

Tissue Mixture Mode
~~~~~~~~~~~~~~~~~~~

- Multi-tissue cfDNA simulation
- Tissue-specific fragment patterns
- Configurable tissue fractions

Cancer Progression Mode
~~~~~~~~~~~~~~~~~~~~~~~

- Simulates tumor progression
- Multiple timepoints
- Increasing tumor fractions
- Mutation accumulation


Available Tissue Profiles
-------------------------

- **hematopoietic**: Blood cell-derived cfDNA
- **liver**: Hepatic cfDNA patterns
- **placenta**: Placental cfDNA (fetal)
- **tumor**: Cancer-derived cfDNA

Output
------

Fragment Files
~~~~~~~~~~~~~~

- ``{output_name}.frag`` - Main synthetic fragment file
- ``{output_name}_timepoint_{time}.frag`` - Timepoint-specific files (cancer progression)
- ``{output_name}_tissue_{tissue}.frag`` - Tissue-specific files (tissue mixture)

