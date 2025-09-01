Simulation
==========

``pyfraglib`` includes a simulation module for generating synthetic cfDNA data. This is useful for testing algorithms, validating methods, and creating benchmarking datasets.

Overview
--------

Simulations are based on the biological processes that underlie DNA fragmentation and shedding into the blood stream. The API reference explains the assumptions (and the scientific references) in greater detail. The simulation module provides different simulation modes (see below). All simulations use genomic sequences from FASTA files to generate realistic fragment end motifs and maintain biological accuracy to a certain degree.

The API integrates with the rest of the ``pyfraglib`` ecosystem:

.. code-block:: python

   from pyfraglib import FragmentSimulator

   simulator = FragmentSimulator(fasta_path="reference.fasta")

   fragments = simulator.simulate_fragments(
       chrom="chr1", start=1_000_000, end=1_100_000,
       num_fragments=10000
   )
   fragments.to_frag_file("synthetic_sample", "output/")


For use in pipelines, we provide a CLI, too:

.. code-block:: bash

   pyfrag.py simulate --config simulation_config.json --out-dir synthetic/


Simulation Modes
----------------

1. Basic Simulation: Generate fragments from a single tissue type with customizable parameters
2. Tissue Mixture Simulation: Simulate mixtures of multiple tissue types with specified fractions
3. Cancer Progression: Simulate tumor fraction changes over time for longitudinal studies

Please refer to the API documentation to learn more about the specific parameters available.

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
       "num_fragments": 50000,
       "simulation_mode": "tissue_mixture",
       "tissue_fractions": {
           "hematopoietic": 0.7,
           "liver": 0.2,
           "tumor": 0.1
       }
   }

Since the API is still under development, this specification might slightly change in the future. Please have a look at the example configuration files to learn more about what is available.
