Fragment Simulator
==================

The :class:`pyfraglib.FragmentSimulator` class provides basic cfDNA fragment simulation capabilities for single-tissue scenarios.

Basic Usage
-----------

.. code-block:: python

   from pyfraglib import FragmentSimulator
   
   # Create simulator instance
   simulator = FragmentSimulator(
       fasta_path="reference.fasta",
       output_name="basic_simulation"
   )
   
   # Define genomic regions
   regions = [
       ("chr1", 1000000, 2000000),
       ("chr2", 3000000, 4000000)
   ]
   
   # Generate fragments
   fragments = simulator.simulate_fragments(regions, n_fragments=10000)
   
   # Save results
   fragments.to_frag_file("basic_simulation", "output/")

Constructor Parameters
----------------------

.. autoclass:: pyfraglib.FragmentSimulator
   :members: __init__
   :no-index:

The constructor accepts the following parameters:

* **fasta_path** (str): Path to reference FASTA file (must be indexed with samtools faidx)
* **output_name** (str): Name for output files
* **length_distribution** (dict, optional): Custom length distribution parameters
* **nuclease_profile** (NucleaseProfile, optional): Custom nuclease cleavage profile

Fragment Generation
-------------------

.. automethod:: pyfraglib.FragmentSimulator.simulate_fragments
   :no-index:

The :meth:`simulate_fragments` method generates fragments with the following characteristics:

* **Realistic lengths**: Based on Gaussian mixture models fitted to real cfDNA data
* **Authentic end motifs**: Extracted from actual genomic sequences
* **Proper coordinates**: 0-based genomic coordinates following BAM conventions
* **Quality flags**: Includes bogus fragment detection

Length Distribution Modeling
----------------------------

Default length distribution parameters:

.. code-block:: python

   default_length_params = {
       "components": [
           {"mean": 145, "std": 20, "weight": 0.3},  # Short fragments
           {"mean": 167, "std": 15, "weight": 0.5},  # Nucleosomal fragments
           {"mean": 190, "std": 30, "weight": 0.2}   # Long fragments
       ],
       "bounds": {
           "min_length": 50,
           "max_length": 500
       }
   }

Custom length distributions:

.. code-block:: python

   custom_length_params = {
       "components": [
           {"mean": 160, "std": 25, "weight": 0.7},
           {"mean": 200, "std": 35, "weight": 0.3}
       ],
       "bounds": {
           "min_length": 100,
           "max_length": 300
       }
   }
   
   simulator = FragmentSimulator(
       fasta_path="reference.fasta",
       output_name="custom_simulation",
       length_distribution=custom_length_params
   )

Nuclease Profile Customization
------------------------------

The :class:`pyfraglib.NucleaseProfile` class allows customization of fragmentation patterns:

.. code-block:: python

   from pyfraglib import NucleaseProfile
   
   # Create custom nuclease profile
   custom_profile = NucleaseProfile(
       name="custom_dnase",
       cleavage_bias={
           "A": 0.3,  # Preference for cleaving after A
           "T": 0.3,  # Preference for cleaving after T
           "G": 0.2,  # Preference for cleaving after G
           "C": 0.2   # Preference for cleaving after C
       },
       motif_preferences={
           "CCCA": 1.5,  # 1.5x preference for this motif
           "TTTC": 1.2,  # 1.2x preference for this motif
           "GGGA": 0.8   # 0.8x preference (disfavored)
       }
   )
   
   simulator = FragmentSimulator(
       fasta_path="reference.fasta",
       output_name="custom_nuclease",
       nuclease_profile=custom_profile
   )

Advanced Configuration
----------------------

Batch simulation with different parameters:

.. code-block:: python

   from pyfraglib import FragmentSimulator
   
   # Define parameter sets
   simulation_params = [
       {
           "name": "short_fragments",
           "length_params": {
               "components": [{"mean": 140, "std": 20, "weight": 1.0}],
               "bounds": {"min_length": 100, "max_length": 200}
           }
       },
       {
           "name": "long_fragments", 
           "length_params": {
               "components": [{"mean": 200, "std": 30, "weight": 1.0}],
               "bounds": {"min_length": 150, "max_length": 300}
           }
       },
       {
           "name": "bimodal_fragments",
           "length_params": {
               "components": [
                   {"mean": 140, "std": 15, "weight": 0.6},
                   {"mean": 200, "std": 25, "weight": 0.4}
               ],
               "bounds": {"min_length": 100, "max_length": 300}
           }
       }
   ]
   
   regions = [("chr1", 1000000, 2000000)]
   
   # Generate different fragment populations
   for params in simulation_params:
       simulator = FragmentSimulator(
           fasta_path="reference.fasta",
           output_name=params["name"],
           length_distribution=params["length_params"]
       )
       
       fragments = simulator.simulate_fragments(regions, n_fragments=10000)
       fragments.to_frag_file(params["name"], "output/")
       
       # Analyze results
       lengths = [f.length for f in fragments]
       print(f"{params['name']}: mean={np.mean(lengths):.1f}, std={np.std(lengths):.1f}")

Quality Control
---------------

Validate simulated fragments:

.. code-block:: python

   from pyfraglib import FragmentSimulator
   from pyfraglib.scores import motif_diversity
   import numpy as np
   
   simulator = FragmentSimulator(
       fasta_path="reference.fasta",
       output_name="validation_test"
   )
   
   regions = [("chr1", 1000000, 2000000)]
   fragments = simulator.simulate_fragments(regions, n_fragments=10000)
   
   # Basic statistics
   print(f"Total fragments: {len(fragments)}")
   print(f"Bogus fragments: {sum(1 for f in fragments if f.is_bogus)}")
   
   # Length distribution validation
   lengths = [f.length for f in fragments if not f.is_bogus]
   print(f"Length stats: mean={np.mean(lengths):.1f}, std={np.std(lengths):.1f}")
   print(f"Length range: {np.min(lengths)}-{np.max(lengths)}")
   
   # End motif validation
   end5p_motifs = [f.end5p[:4] for f in fragments if not f.is_bogus]
   end3p_motifs = [f.end3p[:4] for f in fragments if not f.is_bogus]
   
   print(f"Unique 5' motifs: {len(set(end5p_motifs))}")
   print(f"Unique 3' motifs: {len(set(end3p_motifs))}")
   
   # Motif diversity
   diversity = motif_diversity(fragments, kmer_len=4, index="shannon")
   print(f"Shannon entropy: {diversity:.4f}")
   
   # Chromosome distribution
   chr_distribution = {}
   for fragment in fragments:
       if not fragment.is_bogus:
           chr_distribution[fragment.chrom] = chr_distribution.get(fragment.chrom, 0) + 1
   
   print("Chromosome distribution:")
   for chrom, count in sorted(chr_distribution.items()):
       print(f"  {chrom}: {count}")

Performance Optimization
-------------------------

Memory-efficient simulation for large datasets:

.. code-block:: python

   from pyfraglib import FragmentSimulator
   
   def simulate_large_dataset(fasta_path, regions, total_fragments, chunk_size=10000):
       """Simulate large dataset in chunks to manage memory"""
       
       simulator = FragmentSimulator(
           fasta_path=fasta_path,
           output_name="large_simulation"
       )
       
       all_fragments = []
       processed = 0
       
       while processed < total_fragments:
           current_chunk = min(chunk_size, total_fragments - processed)
           
           fragments = simulator.simulate_fragments(regions, n_fragments=current_chunk)
           all_fragments.extend(fragments)
           processed += current_chunk
           
           print(f"Processed {processed}/{total_fragments} fragments...")
           
           # Optional: save intermediate results
           if processed % 50000 == 0:
               temp_collection = FragmentList(all_fragments)
               temp_collection.to_frag_file(f"temp_{processed}", "temp/")
       
       return all_fragments
   
   # Generate 100,000 fragments in chunks
   regions = [("chr1", 1000000, 10000000)]
   large_fragments = simulate_large_dataset(
       "reference.fasta", regions, 100000, chunk_size=5000
   )

Integration with Analysis Pipeline
----------------------------------

Complete workflow from simulation to analysis:

.. code-block:: python

   from pyfraglib import FragmentSimulator
   from pyfraglib.lengths import fragment_length_plot, fragment_length_gmm
   from pyfraglib.scores import windowed_protection_score
   from pyfraglib.stats import fragments_per_chromosome_barplot
   
   # 1. Generate fragments
   simulator = FragmentSimulator(
       fasta_path="reference.fasta",
       output_name="pipeline_test"
   )
   
   regions = [("chr1", 1000000, 2000000)]
   fragments = simulator.simulate_fragments(regions, n_fragments=25000)
   
   # 2. Save fragments
   fragments.to_frag_file("pipeline_test", "output/")
   
   # 3. Length analysis
   fragment_length_plot(fragments, "output/", "pipeline_test")
   
   # 4. GMM fitting
   with open("gmm_config.json", "w") as f:
       json.dump({
           "number_of_gaussians": 3,
           "subsample_percentage": 0.1,
           "means_lower_bounds": [100, 150, 200],
           "means_upper_bounds": [150, 200, 250],
           "std_lower_bounds": [10, 10, 10],
           "std_upper_bounds": [40, 40, 40],
           "initial_means": [140, 167, 210],
           "initial_stds": [20, 20, 30],
           "initial_weights": [0.3, 0.5, 0.2]
       }, f, indent=2)
   
   fragment_length_gmm(fragments, "gmm_config.json", "output/", "pipeline_test")
   
   # 5. Chromosome distribution
   fragments_per_chromosome_barplot(fragments, "output/", "pipeline_test")
   
   # 6. Calculate scores (if regions file available)
   # wps_scores = windowed_protection_score(fragments, regions_file)
   
   print("Complete analysis pipeline finished!")

Error Handling
--------------

Robust simulation with error handling:

.. code-block:: python

   from pyfraglib import FragmentSimulator
   import logging
   
   def robust_simulation(fasta_path, output_name, regions, n_fragments):
       """Simulation with comprehensive error handling"""
       
       try:
           # Validate inputs
           if not os.path.exists(fasta_path):
               raise FileNotFoundError(f"FASTA file not found: {fasta_path}")
           
           if not os.path.exists(f"{fasta_path}.fai"):
               raise FileNotFoundError(f"FASTA index not found: {fasta_path}.fai")
           
           if not regions:
               raise ValueError("No regions specified")
           
           if n_fragments <= 0:
               raise ValueError("Number of fragments must be positive")
           
           # Create simulator
           simulator = FragmentSimulator(
               fasta_path=fasta_path,
               output_name=output_name
           )
           
           # Generate fragments
           fragments = simulator.simulate_fragments(regions, n_fragments=n_fragments)
           
           # Validate results
           if not fragments:
               raise RuntimeError("No fragments generated")
           
           valid_fragments = [f for f in fragments if not f.is_bogus]
           if len(valid_fragments) < n_fragments * 0.8:
               logging.warning(f"High proportion of bogus fragments: {len(valid_fragments)}/{len(fragments)}")
           
           return fragments
           
       except Exception as e:
           logging.error(f"Simulation failed: {str(e)}")
           raise

Troubleshooting
---------------

Common issues and solutions:

**FASTA file issues:**

.. code-block:: python

   # Check FASTA file and index
   import os
   
   fasta_path = "reference.fasta"
   
   if not os.path.exists(fasta_path):
       print(f"ERROR: FASTA file not found: {fasta_path}")
   elif not os.path.exists(f"{fasta_path}.fai"):
       print(f"ERROR: FASTA index not found. Run: samtools faidx {fasta_path}")
   else:
       print("FASTA file and index OK")

**Memory issues:**

.. code-block:: python

   # Monitor memory usage
   import psutil
   
   process = psutil.Process()
   initial_memory = process.memory_info().rss / 1024 / 1024  # MB
   
   # Simulate with memory monitoring
   fragments = simulator.simulate_fragments(regions, n_fragments=10000)
   
   final_memory = process.memory_info().rss / 1024 / 1024  # MB
   print(f"Memory usage: {initial_memory:.1f} -> {final_memory:.1f} MB")

**Region specification issues:**

.. code-block:: python

   # Validate regions
   def validate_regions(regions):
       for region in regions:
           if len(region) != 3:
               raise ValueError(f"Invalid region format: {region}")
           
           chrom, start, end = region
           if not isinstance(chrom, str):
               raise ValueError(f"Chromosome must be string: {chrom}")
           if not isinstance(start, int) or not isinstance(end, int):
               raise ValueError(f"Coordinates must be integers: {start}, {end}")
           if start >= end:
               raise ValueError(f"Start must be less than end: {start} >= {end}")
           if start < 0:
               raise ValueError(f"Start must be non-negative: {start}")
       
       return True

See Also
--------

* :doc:`overview` - Simulation overview
* :doc:`tissue_mixture` - Tissue mixture simulation
* :doc:`configuration` - Configuration file format
* :class:`pyfraglib.FragmentSimulator` - API reference
* :class:`pyfraglib.NucleaseProfile` - Nuclease profile customization