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

Probabilistic Cleavage Model
----------------------------

The simulator uses a biologically-grounded probabilistic model to determine where cfDNA fragments are cleaved. The per-base cleavage probability integrates multiple biological factors through multiplicative interactions.

The cleavage probability at genomic position :math:`i` is given by:

.. math::
   P'_{\text{cleavage}}(i) = P_{\text{base}}(i) \times F_{\text{nucleosome}}(i) \times F_{\text{tissue}} \times F_{\text{nuclease}}(i) \times F_{\text{tf}}(i)

where each factor represents a different biological process influencing DNA accessibility and nuclease activity. The base accessibility depends on chromatin state:

.. math::
   P_{\text{base}}(i) = \begin{cases}
   0.6 & \text{if position } i \in \text{ open chromatin regions} \\
   0.1 & \text{if position } i \notin \text{ open chromatin regions}
   \end{cases}

Nucleosome positioning provides distance-dependent protection:

.. math::
   F_{\text{nucleosome}}(i) = \begin{cases}
   0.05 + (1 - O_j) \times 0.2 & \text{if } d_{ij} \leq 73 \text{ bp (core)} \\
   0.3 + (1 - O_j) \times 0.4 & \text{if } 73 < d_{ij} \leq 120 \text{ bp (edge)} \\
   0.8 + (1 - O_j) \times 0.2 & \text{if } d_{ij} > 120 \text{ bp (linker)}
   \end{cases}

where :math:`d_{ij}` is the distance from position :math:`i` to the nearest nucleosome center :math:`j`, :math:`O_j` is the occupancy score of nucleosome :math:`j` (range: 0-1). Different tissues exhibit characteristic chromatin accessibility patterns:

.. math::
   F_{\text{tissue}} = \begin{cases}
   1.0 & \text{healthy baseline} \\
   1.2 & \text{hematopoietic (open chromatin)} \\
   0.9 & \text{liver (compact chromatin)} \\
   1.4 & \text{tumor (disrupted chromatin)}
   \end{cases}

The nuclease factor integrates activity levels and sequence preferences across all active nucleases:

.. math::
   F_{\text{nuclease}}(i) = \frac{\sum_{k} A_k \times P_k(i)}{\sum_{k} A_k}

where :math:`k \in \{\text{DNase1}, \text{DNase1L3}, \text{DFFB}\}`, :math:`A_k` is the activity level of nuclease :math:`k`, :math:`P_k(i)` is the sequence preference factor for nuclease :math:`k` at position :math:`i`.

**DNase I and DNase1L3** sequence preferences:

.. math::
   P_k(i) = \prod_{m} \left[1 + (\text{pref}_m - 1) \times \text{freq}_m(s_i)\right]

**DFFB** combines sequence preferences with nucleosome linker preference:

.. math::
   P_{\text{DFFB}}(i) = \left[\prod_{m} \left(1 + (\text{pref}_m - 1) \times \text{freq}_m(s_i)\right)\right] \times \left(2.0 - F_{\text{nucleosome}}(i)\right)

where :math:`m` represents different motifs (e.g., "CC", "AT", "A", "T"), :math:`\text{pref}_m` is the preference value for motif :math:`m` (>1.0 favored, <1.0 disfavored), :math:`\text{freq}_m(s_i)` is the frequency/occurrence of motif :math:`m` in the local sequence context :math:`s_i`, :math:`s_i` is the 20bp sequence window centered at position :math:`i`. Transcription factor binding sites provide protection from nuclease cleavage:

.. math::
   F_{\text{tf}}(i) = \begin{cases}
   0.3 & \text{if position } i \in \text{ TF binding sites} \\
   1.0 & \text{if position } i \notin \text{ TF binding sites}
   \end{cases}

The final cleavage probability is bounded:

.. math::
   P_{\text{cleavage}}(i) = \min(1.0, \max(0.001, P'_{\text{cleavage}}(i)))
