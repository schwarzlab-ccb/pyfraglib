Tissue Mixture Simulation
=========================

The :class:`pyfraglib.TissueMixtureSimulator` class enables simulation of cfDNA from multiple tissue types with specified fractions, modeling the complex mixture of cell-free DNA found in liquid biopsies.

Overview
--------

Tissue mixture simulation is essential for:

* **Liquid biopsy development**: Model realistic cfDNA compositions
* **Tumor fraction detection**: Test sensitivity of ctDNA detection methods
* **Multi-tissue analysis**: Study tissue-of-origin algorithms
* **Clinical validation**: Create control samples with known compositions

Available Tissue Types
----------------------

The simulator supports predefined tissue profiles with distinct fragmentomics characteristics:

.. code-block:: python

   # Available tissue types
   tissue_types = [
       "hematopoietic",  # Blood cells (most abundant in healthy individuals)
       "liver",          # Liver tissue
       "placenta",       # Placental tissue
       "tumor",          # Tumor tissue
       "fetal"           # Fetal tissue (for NIPT applications)
   ]

Each tissue type has characteristic:

* **Fragment length distributions**: Different mean lengths and variances
* **End motif patterns**: Tissue-specific nuclease cleavage preferences
* **Chromatin accessibility**: Varying fragmentation patterns across genomic regions

Basic Usage
-----------

.. code-block:: python

   from pyfraglib import TissueMixtureSimulator
   
   # Create simulator instance
   simulator = TissueMixtureSimulator(
       fasta_path="reference.fasta",
       output_name="tissue_mixture"
   )
   
   # Define tissue composition
   tissue_fractions = {
       "hematopoietic": 0.7,  # 70% blood cells
       "liver": 0.2,          # 20% liver
       "tumor": 0.1           # 10% tumor
   }
   
   # Define genomic regions
   regions = [
       ("chr1", 1000000, 2000000),
       ("chr2", 3000000, 4000000)
   ]
   
   # Generate tissue mixture
   fragments = simulator.simulate_mixture(
       regions, tissue_fractions, n_fragments=50000
   )
   
   # Save results
   fragments.to_frag_file("tissue_mixture", "output/")

Tissue-Specific Characteristics
-------------------------------

Each tissue type has distinct fragmentomics properties:

Hematopoietic (Blood Cells)
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   hematopoietic_profile = {
       "length_distribution": {
           "mean": 167,
           "std": 25,
           "components": [
               {"mean": 145, "std": 20, "weight": 0.3},
               {"mean": 167, "std": 15, "weight": 0.5},
               {"mean": 190, "std": 30, "weight": 0.2}
           ]
       },
       "nuclease_preferences": {
           "CCCA": 1.2,
           "TTTC": 1.1,
           "AAAG": 0.9
       },
       "characteristics": [
           "Most abundant in healthy individuals",
           "Tight nucleosomal packaging",
           "Regular fragmentation patterns"
       ]
   }

Liver Tissue
~~~~~~~~~~~~

.. code-block:: python

   liver_profile = {
       "length_distribution": {
           "mean": 175,
           "std": 30,
           "components": [
               {"mean": 150, "std": 25, "weight": 0.2},
               {"mean": 175, "std": 20, "weight": 0.6},
               {"mean": 200, "std": 35, "weight": 0.2}
           ]
       },
       "nuclease_preferences": {
           "CCCA": 1.3,
           "GGGA": 1.2,
           "TTTT": 0.8
       },
       "characteristics": [
           "Slightly longer mean length",
           "Higher variance in lengths",
           "Metabolically active tissue"
       ]
   }

Tumor Tissue
~~~~~~~~~~~~

.. code-block:: python

   tumor_profile = {
       "length_distribution": {
           "mean": 155,
           "std": 35,
           "components": [
               {"mean": 140, "std": 30, "weight": 0.4},
               {"mean": 165, "std": 25, "weight": 0.4},
               {"mean": 200, "std": 40, "weight": 0.2}
           ]
       },
       "nuclease_preferences": {
           "CCCA": 0.8,
           "TTTC": 1.4,
           "GGGA": 1.1
       },
       "characteristics": [
           "Shorter mean length",
           "Higher fragmentation variability",
           "Altered chromatin structure",
           "Increased DNA damage"
       ]
   }

Placental Tissue
~~~~~~~~~~~~~~~~

.. code-block:: python

   placenta_profile = {
       "length_distribution": {
           "mean": 160,
           "std": 20,
           "components": [
               {"mean": 143, "std": 18, "weight": 0.4},
               {"mean": 160, "std": 15, "weight": 0.5},
               {"mean": 180, "std": 25, "weight": 0.1}
           ]
       },
       "nuclease_preferences": {
           "CCCA": 1.1,
           "TTTC": 1.0,
           "AAAG": 1.2
       },
       "characteristics": [
           "Tight length distribution",
           "Maternal contribution",
           "Unique fragmentation patterns"
       ]
   }

Fetal Tissue
~~~~~~~~~~~~

.. code-block:: python

   fetal_profile = {
       "length_distribution": {
           "mean": 143,
           "std": 18,
           "components": [
               {"mean": 143, "std": 15, "weight": 0.7},
               {"mean": 165, "std": 20, "weight": 0.3}
           ]
       },
       "nuclease_preferences": {
           "CCCA": 1.0,
           "TTTC": 0.9,
           "AAAG": 1.3
       },
       "characteristics": [
           "Shorter fragments",
           "Distinct from maternal DNA",
           "Lower variance",
           "NIPT applications"
       ]
   }

Advanced Configuration
----------------------

Custom tissue profiles can be defined:

.. code-block:: python

   from pyfraglib import TissueMixtureSimulator, TissueProfile
   
   # Create custom tissue profile
   custom_tissue = TissueProfile(
       name="custom_tissue",
       length_params={
           "components": [
               {"mean": 170, "std": 22, "weight": 0.6},
               {"mean": 200, "std": 30, "weight": 0.4}
           ],
           "bounds": {"min_length": 100, "max_length": 400}
       },
       nuclease_preferences={
           "CCCA": 1.5,
           "TTTC": 1.2,
           "GGGA": 0.8,
           "AAAG": 1.1
       }
   )
   
   # Use custom tissue in simulation
   simulator = TissueMixtureSimulator(
       fasta_path="reference.fasta",
       output_name="custom_mixture",
       custom_tissues={"custom_tissue": custom_tissue}
   )
   
   tissue_fractions = {
       "hematopoietic": 0.6,
       "custom_tissue": 0.4
   }
   
   fragments = simulator.simulate_mixture(
       regions, tissue_fractions, n_fragments=30000
   )

Liquid Biopsy Scenarios
-----------------------

Cancer Detection
~~~~~~~~~~~~~~~~

Model different stages of cancer progression:

.. code-block:: python

   # Early stage cancer
   early_cancer = {
       "hematopoietic": 0.95,
       "liver": 0.04,
       "tumor": 0.01
   }
   
   # Advanced cancer
   advanced_cancer = {
       "hematopoietic": 0.70,
       "liver": 0.10,
       "tumor": 0.20
   }
   
   # Metastatic cancer
   metastatic_cancer = {
       "hematopoietic": 0.60,
       "liver": 0.05,
       "tumor": 0.35
   }
   
   scenarios = [
       ("early", early_cancer),
       ("advanced", advanced_cancer),
       ("metastatic", metastatic_cancer)
   ]
   
   for stage, fractions in scenarios:
       simulator = TissueMixtureSimulator(
           fasta_path="reference.fasta",
           output_name=f"cancer_{stage}"
       )
       
       fragments = simulator.simulate_mixture(
           regions, fractions, n_fragments=100000
       )
       
       # Analyze tumor fraction
       tumor_count = sum(1 for f in fragments if f.tissue_origin == "tumor")
       actual_fraction = tumor_count / len(fragments)
       
       print(f"{stage}: expected {fractions['tumor']:.1%}, actual {actual_fraction:.1%}")

Treatment Monitoring
~~~~~~~~~~~~~~~~~~~~

Model response to therapy:

.. code-block:: python

   # Treatment response timeline
   timepoints = [
       {"name": "baseline", "tumor_fraction": 0.20},
       {"name": "week_2", "tumor_fraction": 0.15},
       {"name": "week_4", "tumor_fraction": 0.10},
       {"name": "week_8", "tumor_fraction": 0.05},
       {"name": "week_12", "tumor_fraction": 0.02}
   ]
   
   for timepoint in timepoints:
       tissue_fractions = {
           "hematopoietic": 0.75,
           "liver": 0.05,
           "tumor": timepoint["tumor_fraction"]
       }
       
       # Normalize fractions
       total = sum(tissue_fractions.values())
       tissue_fractions = {k: v/total for k, v in tissue_fractions.items()}
       
       simulator = TissueMixtureSimulator(
           fasta_path="reference.fasta",
           output_name=f"treatment_{timepoint['name']}"
       )
       
       fragments = simulator.simulate_mixture(
           regions, tissue_fractions, n_fragments=50000
       )
       
       print(f"{timepoint['name']}: {len(fragments)} fragments generated")

Validation and Quality Control
------------------------------

Validate tissue mixture composition:

.. code-block:: python

   from pyfraglib import TissueMixtureSimulator
   import numpy as np
   
   def validate_tissue_mixture(fragments, expected_fractions):
       """Validate actual vs expected tissue fractions"""
       
       # Count fragments by tissue origin
       tissue_counts = {}
       for fragment in fragments:
           tissue = fragment.tissue_origin
           tissue_counts[tissue] = tissue_counts.get(tissue, 0) + 1
       
       # Calculate actual fractions
       total_fragments = len(fragments)
       actual_fractions = {
           tissue: count / total_fragments 
           for tissue, count in tissue_counts.items()
       }
       
       # Compare with expected
       print("Tissue fraction validation:")
       for tissue, expected in expected_fractions.items():
           actual = actual_fractions.get(tissue, 0)
           difference = abs(actual - expected)
           print(f"  {tissue}: expected {expected:.1%}, actual {actual:.1%}, diff {difference:.1%}")
       
       return actual_fractions
   
   # Generate mixture
   simulator = TissueMixtureSimulator(
       fasta_path="reference.fasta",
       output_name="validation_test"
   )
   
   expected_fractions = {
       "hematopoietic": 0.7,
       "liver": 0.2,
       "tumor": 0.1
   }
   
   fragments = simulator.simulate_mixture(
       regions, expected_fractions, n_fragments=100000
   )
   
   # Validate composition
   actual_fractions = validate_tissue_mixture(fragments, expected_fractions)

Length Distribution Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Compare length distributions between tissues:

.. code-block:: python

   from pyfraglib import TissueMixtureSimulator
   import matplotlib.pyplot as plt
   import numpy as np
   
   # Generate pure tissue samples
   tissues = ["hematopoietic", "liver", "tumor", "placenta", "fetal"]
   tissue_lengths = {}
   
   for tissue in tissues:
       simulator = TissueMixtureSimulator(
           fasta_path="reference.fasta",
           output_name=f"pure_{tissue}"
       )
       
       # Pure tissue (100% of one type)
       tissue_fractions = {tissue: 1.0}
       
       fragments = simulator.simulate_mixture(
           regions, tissue_fractions, n_fragments=20000
       )
       
       lengths = [f.length for f in fragments if not f.is_bogus]
       tissue_lengths[tissue] = lengths
   
   # Visualize length distributions
   fig, axes = plt.subplots(2, 3, figsize=(15, 10))
   axes = axes.flatten()
   
   for i, (tissue, lengths) in enumerate(tissue_lengths.items()):
       axes[i].hist(lengths, bins=50, alpha=0.7, density=True)
       axes[i].set_title(f"{tissue.capitalize()} Length Distribution")
       axes[i].set_xlabel("Fragment Length (bp)")
       axes[i].set_ylabel("Density")
       axes[i].grid(True, alpha=0.3)
       
       # Add statistics
       mean_length = np.mean(lengths)
       std_length = np.std(lengths)
       axes[i].axvline(mean_length, color='red', linestyle='--', 
                      label=f'Mean: {mean_length:.1f}±{std_length:.1f}')
       axes[i].legend()
   
   # Combined comparison
   for tissue, lengths in tissue_lengths.items():
       axes[5].hist(lengths, bins=50, alpha=0.6, density=True, label=tissue)
   
   axes[5].set_title("All Tissues Comparison")
   axes[5].set_xlabel("Fragment Length (bp)")
   axes[5].set_ylabel("Density")
   axes[5].legend()
   axes[5].grid(True, alpha=0.3)
   
   plt.tight_layout()
   plt.savefig("tissue_length_distributions.png", dpi=300, bbox_inches='tight')
   plt.show()

Performance Optimization
------------------------

Efficient simulation for large datasets:

.. code-block:: python

   from pyfraglib import TissueMixtureSimulator
   import multiprocessing as mp
   from functools import partial
   
   def simulate_sample(params):
       """Simulate single sample with tissue mixture"""
       sample_name, tissue_fractions, n_fragments = params
       
       simulator = TissueMixtureSimulator(
           fasta_path="reference.fasta",
           output_name=sample_name
       )
       
       regions = [("chr1", 1000000, 2000000)]
       fragments = simulator.simulate_mixture(
           regions, tissue_fractions, n_fragments=n_fragments
       )
       
       return sample_name, len(fragments)
   
   # Define sample parameters
   sample_params = [
       ("low_tumor", {"hematopoietic": 0.95, "tumor": 0.05}, 50000),
       ("medium_tumor", {"hematopoietic": 0.85, "tumor": 0.15}, 50000),
       ("high_tumor", {"hematopoietic": 0.70, "tumor": 0.30}, 50000),
   ]
   
   # Process samples in parallel
   with mp.Pool(processes=3) as pool:
       results = pool.map(simulate_sample, sample_params)
   
   # Print results
   for sample_name, fragment_count in results:
       print(f"{sample_name}: {fragment_count} fragments generated")

Error Handling
--------------

Robust tissue mixture simulation:

.. code-block:: python

   from pyfraglib import TissueMixtureSimulator
   import logging
   
   def robust_tissue_simulation(fasta_path, output_name, tissue_fractions, 
                               regions, n_fragments):
       """Tissue mixture simulation with error handling"""
       
       try:
           # Validate tissue fractions
           total_fraction = sum(tissue_fractions.values())
           if abs(total_fraction - 1.0) > 1e-6:
               raise ValueError(f"Tissue fractions must sum to 1.0, got {total_fraction}")
           
           # Check for valid tissue types
           valid_tissues = {"hematopoietic", "liver", "tumor", "placenta", "fetal"}
           invalid_tissues = set(tissue_fractions.keys()) - valid_tissues
           if invalid_tissues:
               raise ValueError(f"Invalid tissue types: {invalid_tissues}")
           
           # Create simulator
           simulator = TissueMixtureSimulator(
               fasta_path=fasta_path,
               output_name=output_name
           )
           
           # Generate mixture
           fragments = simulator.simulate_mixture(
               regions, tissue_fractions, n_fragments=n_fragments
           )
           
           # Validate results
           if not fragments:
               raise RuntimeError("No fragments generated")
           
           # Check tissue distribution
           tissue_counts = {}
           for fragment in fragments:
               tissue = fragment.tissue_origin
               tissue_counts[tissue] = tissue_counts.get(tissue, 0) + 1
           
           for tissue, expected_fraction in tissue_fractions.items():
               actual_count = tissue_counts.get(tissue, 0)
               actual_fraction = actual_count / len(fragments)
               
               if abs(actual_fraction - expected_fraction) > 0.05:
                   logging.warning(f"Tissue fraction deviation: {tissue} "
                                 f"expected {expected_fraction:.1%}, "
                                 f"actual {actual_fraction:.1%}")
           
           return fragments
           
       except Exception as e:
           logging.error(f"Tissue mixture simulation failed: {str(e)}")
           raise

See Also
--------

* :doc:`overview` - Simulation overview
* :doc:`fragment_simulator` - Basic fragment simulation
* :doc:`configuration` - Configuration file format
* :class:`pyfraglib.TissueMixtureSimulator` - API reference
* :class:`pyfraglib.TissueProfile` - Custom tissue profiles