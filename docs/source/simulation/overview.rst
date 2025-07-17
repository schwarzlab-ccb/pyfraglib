Simulation Overview
===================

The pyfraglib simulation module provides comprehensive tools for generating synthetic cell-free DNA (cfDNA) data with realistic biological characteristics. This is essential for method validation, algorithm testing, and benchmark creation.

Key Features
------------

* **Realistic Fragment Properties**: Generates fragments with biologically accurate length distributions, end motifs, and genomic coordinates
* **Tissue-Specific Profiles**: Different fragmentation patterns for various tissue types (blood, liver, placenta, tumor, fetal)
* **Multiple Simulation Modes**: Basic, tissue mixture, cancer progression, and fetal fraction simulations
* **Genomic Sequence Integration**: Uses real FASTA files to generate authentic fragment end motifs
* **Scalable Design**: Generate datasets of any size from small test sets to large cohorts
* **Configurable Parameters**: Extensive customization through JSON configuration files

Simulation Modes
----------------

Basic Simulation
~~~~~~~~~~~~~~~~

The fundamental simulation mode generates cfDNA fragments from a single tissue type:

.. code-block:: python

   from pyfraglib import FragmentSimulator

   simulator = FragmentSimulator(
       fasta_path="reference.fasta",
       output_name="basic_simulation"
   )

   regions = [("chr1", 1000000, 2000000)]
   fragments = simulator.simulate_fragments(regions, n_fragments=10000)

**Use cases:**
- Method development and testing
- Algorithm validation
- Basic benchmark creation
- Educational purposes

Tissue Mixture Simulation
~~~~~~~~~~~~~~~~~~~~~~~~~

Simulates cfDNA from multiple tissue types with specified fractions:

.. code-block:: python

   from pyfraglib import TissueMixtureSimulator

   simulator = TissueMixtureSimulator(
       fasta_path="reference.fasta",
       output_name="tissue_mixture"
   )

   tissue_fractions = {
       "hematopoietic": 0.7,  # Blood cells
       "liver": 0.2,          # Liver tissue
       "tumor": 0.1           # Tumor tissue
   }

   fragments = simulator.simulate_mixture(regions, tissue_fractions, n_fragments=50000)

**Available tissue types:**
- ``hematopoietic``: Blood cells (most abundant in healthy individuals)
- ``liver``: Liver tissue
- ``placenta``: Placental tissue
- ``tumor``: Tumor tissue
- ``fetal``: Fetal tissue (for NIPT applications)

**Use cases:**
- Liquid biopsy method development
- Tumor fraction detection algorithms
- Tissue-of-origin analysis
- Multi-tissue biomarker discovery

Cancer Progression Simulation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Models tumor fraction changes over time for longitudinal studies:

.. code-block:: python

   # Configuration for cancer progression
   timepoints = [
       {"name": "baseline", "tumor_fraction": 0.01},
       {"name": "month_3", "tumor_fraction": 0.05},
       {"name": "month_6", "tumor_fraction": 0.15},
       {"name": "month_12", "tumor_fraction": 0.30}
   ]

   for timepoint in timepoints:
       tissue_fractions = {
           "hematopoietic": 1.0 - timepoint["tumor_fraction"],
           "tumor": timepoint["tumor_fraction"]
       }

       fragments = simulator.simulate_mixture(regions, tissue_fractions, n_fragments=25000)

**Use cases:**
- Longitudinal study design
- Minimal residual disease detection
- Treatment response monitoring
- Disease progression modeling

Fetal Fraction Simulation
~~~~~~~~~~~~~~~~~~~~~~~~~

Specialized simulation for non-invasive prenatal testing (NIPT):

.. code-block:: python

   fetal_fractions = [0.05, 0.10, 0.15, 0.20, 0.25]

   for fetal_fraction in fetal_fractions:
       tissue_fractions = {
           "hematopoietic": 0.9 - fetal_fraction,  # Maternal blood
           "placenta": 0.1,                        # Maternal placenta
           "fetal": fetal_fraction                 # Fetal contribution
       }

       fragments = simulator.simulate_mixture(regions, tissue_fractions, n_fragments=100000)

**Use cases:**
- NIPT algorithm validation
- Fetal fraction estimation methods
- Aneuploidy detection algorithms
- Quality control for prenatal screening

Biological Realism
------------------

Fragment Length Distributions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Each tissue type has characteristic fragment length distributions modeled as Gaussian mixtures:

.. code-block:: python

   # Example tissue-specific parameters
   tissue_profiles = {
       "hematopoietic": {
           "mean_length": 167,
           "std_length": 25,
           "components": [
               {"mean": 145, "std": 20, "weight": 0.3},
               {"mean": 167, "std": 15, "weight": 0.5},
               {"mean": 190, "std": 30, "weight": 0.2}
           ]
       },
       "liver": {
           "mean_length": 175,
           "std_length": 30,
           "components": [
               {"mean": 150, "std": 25, "weight": 0.2},
               {"mean": 175, "std": 20, "weight": 0.6},
               {"mean": 200, "std": 35, "weight": 0.2}
           ]
       }
   }

End Motif Patterns
~~~~~~~~~~~~~~~~~~

Fragment end motifs (5' and 3' ends) are generated from real genomic sequences:

- **Sequence-based generation**: Reads actual FASTA sequences to create realistic end motifs
- **Tissue-specific preferences**: Different tissues may have subtle preferences for certain motifs
- **Nuclease activity modeling**: Incorporates known nuclease cleavage preferences

Genomic Distribution
~~~~~~~~~~~~~~~~~~~~

Fragments are distributed across the genome based on:

- **Accessibility**: Open chromatin regions are more fragmented
- **Nucleosome positioning**: Periodic fragmentation patterns around nucleosomes
- **Tissue-specific chromatin states**: Different tissues have different chromatin accessibility

Quality Control
---------------

Validation Metrics
~~~~~~~~~~~~~~~~~~

The simulator includes built-in quality control metrics:

.. code-block:: python

   from pyfraglib.scores import motif_diversity
   from pyfraglib.math import goodness_of_fit_stats

   # Validate fragment length distribution
   lengths = [f.length for f in fragments]

   # Compare with expected distribution
   expected_params = [167, 25, 1.0]  # mean, std, weight
   gof_stats = goodness_of_fit_stats(lengths, expected_params, 1)

   print(f"KS statistic: {gof_stats['kolmogorov_smirnov_statistic']:.4f}")
   print(f"Wasserstein distance: {gof_stats['wasserstein_distance']:.4f}")

   # Validate motif diversity
   diversity = motif_diversity(fragments, kmer_len=4, index="shannon")
   print(f"Shannon entropy: {diversity:.4f}")

Comparison with Real Data
~~~~~~~~~~~~~~~~~~~~~~~~~

Best practices for validating simulated data:

.. code-block:: python

   from pyfraglib import Fragment
   import numpy as np
   from scipy import stats

   # Load real data for comparison
   real_fragments = Fragment.from_bam("real_sample.bam", "real_variants.vcf")
   real_lengths = [f.length for f in real_fragments if not f.is_bogus]

   # Compare with simulated data
   sim_lengths = [f.length for f in simulated_fragments]

   # Statistical comparison
   ks_stat, ks_p = stats.ks_2samp(real_lengths, sim_lengths)
   print(f"KS test: D={ks_stat:.4f}, p={ks_p:.2e}")

   # Compare summary statistics
   print(f"Real data: mean={np.mean(real_lengths):.1f}, std={np.std(real_lengths):.1f}")
   print(f"Simulated: mean={np.mean(sim_lengths):.1f}, std={np.std(sim_lengths):.1f}")

Performance Considerations
--------------------------

Memory Management
~~~~~~~~~~~~~~~~~

For large-scale simulations:

.. code-block:: python

   # Process in chunks for memory efficiency
   def simulate_large_dataset(regions, n_fragments_total, chunk_size=10000):
       all_fragments = []

       for i in range(0, n_fragments_total, chunk_size):
           chunk_size_actual = min(chunk_size, n_fragments_total - i)

           fragments = simulator.simulate_fragments(regions, n_fragments=chunk_size_actual)
           all_fragments.extend(fragments)

           # Optional: save intermediate results
           if i % 50000 == 0:
               print(f"Processed {i + chunk_size_actual} fragments...")

       return all_fragments

Parallel Processing
~~~~~~~~~~~~~~~~~~~

For multiple samples or parameter sweeps:

.. code-block:: python

   from multiprocessing import Pool
   from functools import partial

   def simulate_single_sample(params):
       sample_name, fasta_path, tissue_fractions, n_fragments = params

       simulator = TissueMixtureSimulator(
           fasta_path=fasta_path,
           output_name=sample_name
       )

       regions = [("chr1", 1000000, 2000000)]
       fragments = simulator.simulate_mixture(regions, tissue_fractions, n_fragments)

       return sample_name, fragments

   # Define parameter sets
   param_sets = [
       ("sample_1", "reference.fasta", {"hematopoietic": 0.9, "tumor": 0.1}, 10000),
       ("sample_2", "reference.fasta", {"hematopoietic": 0.8, "tumor": 0.2}, 10000),
       ("sample_3", "reference.fasta", {"hematopoietic": 0.7, "tumor": 0.3}, 10000),
   ]

   # Process in parallel
   with Pool(processes=4) as pool:
       results = pool.map(simulate_single_sample, param_sets)

Applications
------------

Method Development
~~~~~~~~~~~~~~~~~~

- **Algorithm validation**: Test new fragmentomics methods with known ground truth
- **Parameter optimization**: Explore algorithm sensitivity to various conditions
- **Benchmark creation**: Generate standardized datasets for method comparison
- **Proof of concept**: Demonstrate feasibility of new analytical approaches

Clinical Applications
~~~~~~~~~~~~~~~~~~~~~

- **Assay development**: Design and validate clinical tests
- **Quality control**: Create control samples for routine testing
- **Training datasets**: Generate data for machine learning applications
- **Sensitivity analysis**: Determine detection limits for clinical markers

Research Applications
~~~~~~~~~~~~~~~~~~~~~

- **Hypothesis testing**: Generate data to test specific biological hypotheses
- **Study design**: Power analysis and sample size calculations
- **Protocol optimization**: Optimize experimental protocols
- **Educational resources**: Create datasets for training and teaching

See Also
--------

* :doc:`fragment_simulator` - Basic fragment simulation
* :doc:`tissue_mixture` - Tissue mixture simulation details
* :doc:`configuration` - Configuration file format
* :doc:`../examples/simulation_examples` - Practical simulation examples
