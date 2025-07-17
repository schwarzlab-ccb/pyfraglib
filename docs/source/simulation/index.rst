Simulation
==========

pyfraglib includes a comprehensive simulation module for generating synthetic cfDNA data. This is useful for testing algorithms, validating methods, and creating benchmark datasets.

.. toctree::
   :maxdepth: 2

   overview
   fragment_simulator
   tissue_mixture
   configuration

Overview
--------

The simulation module provides several simulation modes:

1. **Basic Simulation**: Single-tissue cfDNA simulation
2. **Tissue Mixture**: Multi-tissue mixture simulation  
3. **Cancer Progression**: Tumor fraction progression over time
4. **Fetal Fraction**: NIPT simulation with fetal DNA contribution

All simulations use real genomic sequences from FASTA files to generate realistic fragment end motifs and maintain biological accuracy.

Quick Start
-----------

.. code-block:: python

   from pyfraglib import FragmentSimulator
   
   # Create simulator
   simulator = FragmentSimulator(
       fasta_path="reference.fasta",
       output_name="synthetic_sample"
   )
   
   # Define regions
   regions = [("chr1", 1000000, 2000000)]
   
   # Generate fragments
   fragments = simulator.simulate_fragments(regions, n_fragments=10000)
   
   # Save results
   fragments.to_frag_file("synthetic_sample", "output/")

Command Line Usage
------------------

.. code-block:: bash

   pyfrag.py simulate --config simulation_config.json --out-dir synthetic/

Key Features
------------

* **Realistic Fragment Properties**: Length distributions, end motifs, and genomic coordinates
* **Tissue-Specific Profiles**: Different fragmentation patterns for various tissue types
* **Mutation Simulation**: Optional variant incorporation
* **Scalable**: Generate datasets of any size
* **Configurable**: Extensive configuration options through JSON files

Simulation Modes
----------------

Basic Simulation
~~~~~~~~~~~~~~~~

Generate fragments from a single tissue type with customizable parameters.

Tissue Mixture Simulation
~~~~~~~~~~~~~~~~~~~~~~~~~~

Simulate mixtures of multiple tissue types with specified fractions:

- Hematopoietic (blood cells)
- Liver
- Placenta  
- Tumor
- Fetal

Cancer Progression
~~~~~~~~~~~~~~~~~~

Simulate tumor fraction changes over time for longitudinal studies.

Fetal Fraction
~~~~~~~~~~~~~~

Specialized simulation for non-invasive prenatal testing (NIPT) applications.

Configuration
-------------

Simulations are configured through JSON files with the following structure:

.. code-block:: json

   {
       "fasta_path": "/path/to/reference.fasta",
       "output_name": "simulation_output",
       "regions": [
           {"chr": "chr1", "start": 1000000, "end": 2000000}
       ],
       "n_fragments": 50000,
       "simulation_mode": "tissue_mixture",
       "tissue_fractions": {
           "hematopoietic": 0.7,
           "liver": 0.2,
           "tumor": 0.1
       }
   }

Applications
------------

* **Method Validation**: Test fragmentomics algorithms with known ground truth
* **Benchmark Creation**: Generate standardized datasets for comparison
* **Parameter Optimization**: Explore algorithm sensitivity to various conditions
* **Protocol Development**: Design and validate new analysis workflows
* **Training Data**: Create datasets for machine learning applications