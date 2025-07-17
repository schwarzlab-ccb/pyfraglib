Simulate Command
================

The ``simulate`` command generates synthetic cfDNA fragment data based on biological parameters and real genomic sequences.

Syntax
------

.. code-block:: bash

   pyfrag.py simulate [OPTIONS]

Options
-------

.. code-block:: bash

   --config PATH          Configuration file (JSON format) (required)
   --out-dir PATH         Output directory for synthetic fragments (required)
   --mode MODE           Simulation mode (basic, tissue_mixture, cancer_progression, fetal_fraction)

Examples
--------

Basic Simulation
~~~~~~~~~~~~~~~~

.. code-block:: bash

   pyfrag.py simulate --config configs/simulation_example.json --out-dir synthetic/

Tissue Mixture Simulation
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   pyfrag.py simulate --config configs/simulation_tissue_mixture.json --mode tissue_mixture --out-dir synthetic/

Cancer Progression Simulation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   pyfrag.py simulate --config configs/simulation_cancer_progression.json --mode cancer_progression --out-dir synthetic/

Configuration File Format
--------------------------

Basic Configuration
~~~~~~~~~~~~~~~~~~~

.. code-block:: json

   {
       "fasta_path": "/path/to/reference.fasta",
       "output_name": "synthetic_sample",
       "regions": [
           {"chr": "chr1", "start": 1000000, "end": 2000000},
           {"chr": "chr2", "start": 500000, "end": 1500000}
       ],
       "n_fragments": 50000,
       "fragment_length_distribution": {
           "mean": 167,
           "std": 15
       },
       "simulation_mode": "basic"
   }

Tissue Mixture Configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: json

   {
       "fasta_path": "/path/to/reference.fasta",
       "output_name": "tissue_mixture",
       "regions": [
           {"chr": "chr1", "start": 1000000, "end": 2000000}
       ],
       "n_fragments": 100000,
       "simulation_mode": "tissue_mixture",
       "tissue_fractions": {
           "hematopoietic": 0.6,
           "liver": 0.2,
           "placenta": 0.15,
           "tumor": 0.05
       }
   }

Cancer Progression Configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: json

   {
       "fasta_path": "/path/to/reference.fasta",
       "output_name": "cancer_progression",
       "regions": [
           {"chr": "chr1", "start": 1000000, "end": 2000000}
       ],
       "n_fragments": 75000,
       "simulation_mode": "cancer_progression",
       "timepoints": [
           {"time": "baseline", "tumor_fraction": 0.0},
           {"time": "early", "tumor_fraction": 0.02},
           {"time": "progression", "tumor_fraction": 0.1},
           {"time": "advanced", "tumor_fraction": 0.3}
       ]
   }

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
- Realistic cfDNA composition

Cancer Progression Mode
~~~~~~~~~~~~~~~~~~~~~~~

- Simulates tumor evolution
- Multiple timepoints
- Increasing tumor fractions
- Mutation accumulation

Fetal Fraction Mode
~~~~~~~~~~~~~~~~~~~

- Maternal-fetal cfDNA mixture
- Configurable fetal fraction
- Fetal-specific fragment patterns
- NIPT applications

Available Tissue Profiles
-------------------------

- **hematopoietic**: Blood cell-derived cfDNA
- **liver**: Hepatic cfDNA patterns
- **placenta**: Placental cfDNA (fetal)
- **tumor**: Cancer-derived cfDNA
- **fetal**: Fetal cfDNA patterns

Output
------

Fragment Files
~~~~~~~~~~~~~~

- ``{output_name}.frag`` - Main synthetic fragment file
- ``{output_name}_timepoint_{time}.frag`` - Timepoint-specific files (cancer progression)
- ``{output_name}_tissue_{tissue}.frag`` - Tissue-specific files (tissue mixture)

Metadata Files
~~~~~~~~~~~~~~

- ``{output_name}_metadata.json`` - Simulation parameters and statistics
- ``{output_name}_summary.txt`` - Human-readable summary
- ``{output_name}_validation.csv`` - Quality control metrics

Requirements
------------

- Reference FASTA file must be indexed (samtools faidx)
- Sufficient disk space for output files
- Valid genomic regions in configuration
- Appropriate tissue fractions (must sum to 1.0)

Performance Notes
-----------------

- Simulation time scales with fragment count
- Memory usage depends on region size
- FASTA indexing improves performance
- Large-scale simulations may require significant time

Validation
----------

Quality Control Metrics
~~~~~~~~~~~~~~~~~~~~~~~

- Fragment length distribution validation
- End motif frequency validation
- Tissue fraction verification
- Genomic region coverage

Comparison with Real Data
~~~~~~~~~~~~~~~~~~~~~~~~~

- Length distribution similarity
- Motif diversity comparison
- Fragmentomics score validation
- Tissue-specific pattern verification

Troubleshooting
---------------

**"FASTA file not indexed" error**
   Create an index with: ``samtools faidx reference.fasta``

**"Invalid tissue fractions" error**
   Ensure tissue fractions sum to 1.0

**Memory errors**
   Reduce the number of fragments or region size

**Poor simulation quality**
   Adjust nuclease profiles or tissue parameters