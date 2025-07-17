Simulation Examples
===================

Generate and analyze synthetic cfDNA data using the pyfraglib simulation module:

.. code-block:: python

   from pyfraglib import FragmentSimulator, TissueMixtureSimulator
   from pyfraglib.lengths import fragment_length_plot
   import numpy as np
   import matplotlib.pyplot as plt
   
   # Basic fragment simulation
   simulator = FragmentSimulator(
       fasta_path="reference.fasta",
       output_name="synthetic_basic"
   )
   
   regions = [("chr1", 1000000, 2000000)]
   fragments = simulator.simulate_fragments(regions, n_fragments=10000)
   
   # Analyze simulated data
   lengths = [f.length for f in fragments]
   print(f"Simulated {len(fragments)} fragments")
   print(f"Mean length: {np.mean(lengths):.1f} bp")
   print(f"Length range: {np.min(lengths)} - {np.max(lengths)} bp")

Tissue Mixture Simulation
-------------------------

Simulate cfDNA from multiple tissue types:

.. code-block:: python

   from pyfraglib import TissueMixtureSimulator
   import json
   
   # Create tissue mixture configuration
   config = {
       "fasta_path": "reference.fasta",
       "output_name": "tissue_mixture",
       "regions": [
           {"chr": "chr1", "start": 1000000, "end": 2000000},
           {"chr": "chr2", "start": 3000000, "end": 4000000}
       ],
       "n_fragments": 50000,
       "simulation_mode": "tissue_mixture",
       "tissue_fractions": {
           "hematopoietic": 0.7,
           "liver": 0.2,
           "tumor": 0.1
       }
   }
   
   # Save configuration
   with open("tissue_mixture_config.json", "w") as f:
       json.dump(config, f, indent=2)
   
   # Run simulation
   simulator = TissueMixtureSimulator(
       fasta_path="reference.fasta",
       output_name="tissue_mixture"
   )
   
   tissue_fractions = {
       "hematopoietic": 0.7,
       "liver": 0.2,
       "tumor": 0.1
   }
   
   fragments = simulator.simulate_mixture(
       regions, tissue_fractions, n_fragments=50000
   )
   
   # Analyze tissue-specific characteristics
   print(f"Total fragments: {len(fragments)}")
   print(f"Expected hematopoietic: {50000 * 0.7}")
   print(f"Expected liver: {50000 * 0.2}")
   print(f"Expected tumor: {50000 * 0.1}")

Cancer Progression Simulation
-----------------------------

Simulate increasing tumor fraction over time:

.. code-block:: python

   from pyfraglib import TissueMixtureSimulator
   import json
   import matplotlib.pyplot as plt
   
   # Create cancer progression configuration
   config = {
       "fasta_path": "reference.fasta",
       "output_name": "cancer_progression",
       "regions": [
           {"chr": "chr1", "start": 1000000, "end": 2000000}
       ],
       "n_fragments": 25000,
       "simulation_mode": "cancer_progression",
       "timepoints": [
           {"name": "baseline", "tumor_fraction": 0.01},
           {"name": "month_3", "tumor_fraction": 0.05},
           {"name": "month_6", "tumor_fraction": 0.15},
           {"name": "month_12", "tumor_fraction": 0.30}
       ]
   }
   
   # Save configuration
   with open("cancer_progression_config.json", "w") as f:
       json.dump(config, f, indent=2)
   
   # Simulate each timepoint
   timepoints = []
   tumor_fractions = []
   fragment_counts = []
   
   for timepoint in config["timepoints"]:
       simulator = TissueMixtureSimulator(
           fasta_path="reference.fasta",
           output_name=f"cancer_{timepoint['name']}"
       )
       
       tissue_fractions = {
           "hematopoietic": 1.0 - timepoint["tumor_fraction"],
           "tumor": timepoint["tumor_fraction"]
       }
       
       fragments = simulator.simulate_mixture(
           regions, tissue_fractions, n_fragments=25000
       )
       
       timepoints.append(timepoint["name"])
       tumor_fractions.append(timepoint["tumor_fraction"])
       fragment_counts.append(len(fragments))
       
       print(f"{timepoint['name']}: {len(fragments)} fragments, {timepoint['tumor_fraction']:.1%} tumor")
   
   # Visualize progression
   fig, axes = plt.subplots(1, 2, figsize=(12, 5))
   
   axes[0].plot(range(len(timepoints)), tumor_fractions, 'o-')
   axes[0].set_xlabel('Timepoint')
   axes[0].set_ylabel('Tumor Fraction')
   axes[0].set_title('Tumor Fraction Over Time')
   axes[0].set_xticks(range(len(timepoints)))
   axes[0].set_xticklabels(timepoints, rotation=45)
   axes[0].grid(True)
   
   axes[1].bar(timepoints, fragment_counts)
   axes[1].set_xlabel('Timepoint')
   axes[1].set_ylabel('Fragment Count')
   axes[1].set_title('Fragment Count by Timepoint')
   axes[1].tick_params(axis='x', rotation=45)
   
   plt.tight_layout()
   plt.savefig('cancer_progression.png', dpi=300, bbox_inches='tight')
   plt.show()

Fetal Fraction Simulation
--------------------------

Simulate NIPT samples with varying fetal fractions:

.. code-block:: python

   from pyfraglib import TissueMixtureSimulator
   import json
   
   # Create fetal fraction configuration
   config = {
       "fasta_path": "reference.fasta",
       "output_name": "fetal_fraction",
       "regions": [
           {"chr": "chr21", "start": 1000000, "end": 2000000}
       ],
       "n_fragments": 100000,
       "simulation_mode": "fetal_fraction",
       "fetal_fractions": [0.05, 0.10, 0.15, 0.20, 0.25]
   }
   
   # Save configuration
   with open("fetal_fraction_config.json", "w") as f:
       json.dump(config, f, indent=2)
   
   # Simulate different fetal fractions
   results = []
   
   for fetal_fraction in config["fetal_fractions"]:
       simulator = TissueMixtureSimulator(
           fasta_path="reference.fasta",
           output_name=f"fetal_{fetal_fraction:.0%}"
       )
       
       tissue_fractions = {
           "hematopoietic": 0.9 - fetal_fraction,  # Maternal blood
           "placenta": 0.1,  # Maternal placenta
           "fetal": fetal_fraction  # Fetal contribution
       }
       
       fragments = simulator.simulate_mixture(
           regions, tissue_fractions, n_fragments=100000
       )
       
       results.append({
           'fetal_fraction': fetal_fraction,
           'total_fragments': len(fragments),
           'expected_fetal': int(100000 * fetal_fraction),
           'expected_maternal': int(100000 * (1.0 - fetal_fraction))
       })
       
       print(f"Fetal fraction {fetal_fraction:.0%}: {len(fragments)} fragments")
   
   # Analyze results
   import pandas as pd
   
   df = pd.DataFrame(results)
   print("\nFetal Fraction Simulation Results:")
   print(df)

Command Line Simulation
-----------------------

Use the command line interface for large-scale simulations:

.. code-block:: bash

   # Basic simulation
   pyfrag.py simulate --config simulation_basic.json --out-dir synthetic_basic/
   
   # Tissue mixture simulation
   pyfrag.py simulate --config tissue_mixture_config.json --out-dir synthetic_mixture/
   
   # Cancer progression simulation
   pyfrag.py simulate --config cancer_progression_config.json --out-dir synthetic_cancer/
   
   # Fetal fraction simulation
   pyfrag.py simulate --config fetal_fraction_config.json --out-dir synthetic_fetal/

Configuration Files
-------------------

Create comprehensive configuration files for different simulation scenarios:

.. code-block:: python

   import json
   
   # Basic simulation configuration
   basic_config = {
       "fasta_path": "/path/to/reference.fasta",
       "output_name": "synthetic_basic",
       "regions": [
           {"chr": "chr1", "start": 1000000, "end": 2000000},
           {"chr": "chr2", "start": 3000000, "end": 4000000}
       ],
       "n_fragments": 25000,
       "simulation_mode": "basic"
   }
   
   # Advanced tissue mixture configuration
   advanced_config = {
       "fasta_path": "/path/to/reference.fasta",
       "output_name": "synthetic_advanced",
       "regions": [
           {"chr": "chr1", "start": 1000000, "end": 5000000},
           {"chr": "chr2", "start": 1000000, "end": 5000000},
           {"chr": "chr3", "start": 1000000, "end": 5000000}
       ],
       "n_fragments": 100000,
       "simulation_mode": "tissue_mixture",
       "tissue_fractions": {
           "hematopoietic": 0.60,
           "liver": 0.25,
           "placenta": 0.10,
           "tumor": 0.05
       },
       "nuclease_profiles": {
           "hematopoietic": {"mean_length": 167, "std_length": 25},
           "liver": {"mean_length": 175, "std_length": 30},
           "placenta": {"mean_length": 160, "std_length": 20},
           "tumor": {"mean_length": 155, "std_length": 35}
       }
   }
   
   # Save configurations
   with open("simulation_basic.json", "w") as f:
       json.dump(basic_config, f, indent=2)
   
   with open("simulation_advanced.json", "w") as f:
       json.dump(advanced_config, f, indent=2)
   
   print("Configuration files created:")
   print("- simulation_basic.json")
   print("- simulation_advanced.json")

Validation and Analysis
-----------------------

Validate simulated data against real data:

.. code-block:: python

   from pyfraglib import Fragment, FragmentSimulator
   from pyfraglib.lengths import fragment_length_plot
   from pyfraglib.scores import motif_diversity
   import numpy as np
   import matplotlib.pyplot as plt
   from scipy import stats
   
   # Load real data
   real_fragments = Fragment.from_bam("real_sample.bam", "real_variants.vcf")
   real_lengths = [f.length for f in real_fragments if not f.is_bogus]
   
   # Generate simulated data
   simulator = FragmentSimulator(
       fasta_path="reference.fasta",
       output_name="validation"
   )
   
   regions = [("chr1", 1000000, 2000000)]
   sim_fragments = simulator.simulate_fragments(regions, n_fragments=len(real_lengths))
   sim_lengths = [f.length for f in sim_fragments]
   
   # Compare distributions
   print("Real data statistics:")
   print(f"  Mean: {np.mean(real_lengths):.1f} bp")
   print(f"  Std: {np.std(real_lengths):.1f} bp")
   print(f"  Median: {np.median(real_lengths):.1f} bp")
   
   print("\nSimulated data statistics:")
   print(f"  Mean: {np.mean(sim_lengths):.1f} bp")
   print(f"  Std: {np.std(sim_lengths):.1f} bp")
   print(f"  Median: {np.median(sim_lengths):.1f} bp")
   
   # Statistical tests
   ks_stat, ks_p = stats.ks_2samp(real_lengths, sim_lengths)
   print(f"\nKolmogorov-Smirnov test: D={ks_stat:.4f}, p={ks_p:.2e}")
   
   # Compare motif diversity
   real_diversity = motif_diversity(real_fragments, kmer_len=4, index="shannon")
   sim_diversity = motif_diversity(sim_fragments, kmer_len=4, index="shannon")
   
   print(f"\nMotif diversity (Shannon):")
   print(f"  Real: {real_diversity:.4f}")
   print(f"  Simulated: {sim_diversity:.4f}")
   print(f"  Difference: {abs(real_diversity - sim_diversity):.4f}")
   
   # Visualization
   fig, axes = plt.subplots(2, 2, figsize=(12, 10))
   
   # Length distributions
   axes[0, 0].hist(real_lengths, bins=50, alpha=0.7, label='Real', density=True)
   axes[0, 0].hist(sim_lengths, bins=50, alpha=0.7, label='Simulated', density=True)
   axes[0, 0].set_xlabel('Fragment Length (bp)')
   axes[0, 0].set_ylabel('Density')
   axes[0, 0].set_title('Length Distribution Comparison')
   axes[0, 0].legend()
   
   # Q-Q plot
   stats.probplot(real_lengths, dist="norm", plot=axes[0, 1])
   axes[0, 1].set_title('Q-Q Plot: Real Data')
   
   # Cumulative distributions
   real_sorted = np.sort(real_lengths)
   sim_sorted = np.sort(sim_lengths)
   real_cdf = np.arange(1, len(real_sorted) + 1) / len(real_sorted)
   sim_cdf = np.arange(1, len(sim_sorted) + 1) / len(sim_sorted)
   
   axes[1, 0].plot(real_sorted, real_cdf, label='Real')
   axes[1, 0].plot(sim_sorted, sim_cdf, label='Simulated')
   axes[1, 0].set_xlabel('Fragment Length (bp)')
   axes[1, 0].set_ylabel('Cumulative Probability')
   axes[1, 0].set_title('Cumulative Distribution Comparison')
   axes[1, 0].legend()
   
   # Difference plot
   axes[1, 1].scatter(real_lengths[:1000], sim_lengths[:1000], alpha=0.5)
   axes[1, 1].plot([min(real_lengths), max(real_lengths)], 
                   [min(real_lengths), max(real_lengths)], 'r--')
   axes[1, 1].set_xlabel('Real Fragment Length (bp)')
   axes[1, 1].set_ylabel('Simulated Fragment Length (bp)')
   axes[1, 1].set_title('Length Correlation')
   
   plt.tight_layout()
   plt.savefig('simulation_validation.png', dpi=300, bbox_inches='tight')
   plt.show()

See Also
--------

* :doc:`basic_workflow` - Basic fragmentomics workflow
* :doc:`length_analysis` - Fragment length analysis
* :doc:`batch_processing` - Processing multiple samples
* :class:`pyfraglib.FragmentSimulator` - Basic fragment simulator
* :class:`pyfraglib.TissueMixtureSimulator` - Tissue mixture simulator