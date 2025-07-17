Examples
========

This section provides comprehensive examples of using pyfraglib for various fragmentomics analyses.

.. toctree::
   :maxdepth: 2

   basic_workflow
   mutation_analysis
   fragmentomics_scores
   length_analysis
   simulation_examples
   batch_processing

Basic Workflow
--------------

Complete example of processing BAM files through fragmentomics analysis:

.. code-block:: python

   from pyfraglib import Fragment, FragFile, fit_gmm, plot_gmm
   import matplotlib.pyplot as plt
   
   # 1. Extract fragments from BAM file
   fragments = Fragment.from_bam("sample.bam", "variants.vcf")
   
   # 2. Save to fragment file
   fragments.to_frag_file("sample", "output/")
   
   # 3. Load fragment file
   frag_file = FragFile("output/sample.frag")
   fragments = frag_file.get_fragment_list()
   
   # 4. Basic statistics
   print(f"Total fragments: {fragments.length()}")
   print(f"Bogus fragments: {fragments.count_bogus_fragments()}")
   print(f"Mutated fragments: {fragments.count_mutated_fragments()}")
   
   # 5. Length distribution analysis
   lengths = [f.length for f in fragments if not f.is_bogus]
   gmm_result = fit_gmm(lengths, n_components=3)
   plot_gmm(lengths, gmm_result, "length_distribution.png")
   
   frag_file.close()

Mutation Analysis
-----------------

Analyze fragment mutation patterns:

.. code-block:: python

   from pyfraglib import Fragment, FragFile
   
   # Extract fragments with mutation annotation
   fragments = Fragment.from_bam("tumor.bam", "somatic_variants.vcf")
   
   # Analyze mutation patterns
   total_fragments = fragments.length()
   mutated_fragments = fragments.count_mutated_fragments()
   mutation_rate = mutated_fragments / total_fragments
   
   print(f"Mutation rate: {mutation_rate:.2%}")
   
   # Analyze mutation by chromosome
   mutation_by_chr = {}
   for fragment in fragments:
       if fragment.is_mutated:
           chr_name = fragment.chrom
           mutation_by_chr[chr_name] = mutation_by_chr.get(chr_name, 0) + 1
   
   print("Mutations by chromosome:")
   for chr_name, count in sorted(mutation_by_chr.items()):
       print(f"  {chr_name}: {count}")

Fragmentomics Scores
--------------------

Calculate windowed protection scores (WPS) and motif diversity:

.. code-block:: python

   from pyfraglib.scores import windowed_protection_score, calculate_motif_diversity
   
   # Define genomic regions of interest
   bed_regions = [
       ("chr1", 1000000, 1001000),
       ("chr2", 2000000, 2001000),
       ("chr3", 3000000, 3001000)
   ]
   
   # Calculate WPS scores
   wps_scores = windowed_protection_score(fragments, bed_regions)
   
   print("WPS scores:")
   for region, score in zip(bed_regions, wps_scores):
       chr_name, start, end = region
       print(f"  {chr_name}:{start}-{end}: {score:.3f}")
   
   # Calculate motif diversity
   diversity = calculate_motif_diversity(fragments, kmer_length=4)
   print(f"Shannon entropy: {diversity['shannon']:.3f}")
   print(f"Simpson index: {diversity['simpson']:.3f}")

Batch Processing
----------------

Process multiple samples efficiently:

.. code-block:: python

   from pyfraglib import Fragment
   import os
   
   # Process multiple BAM files
   bam_files = ["sample1.bam", "sample2.bam", "sample3.bam"]
   vcf_files = ["sample1.vcf", "sample2.vcf", "sample3.vcf"]
   
   # Use parallel processing
   Fragment.bams_to_frags(bam_files, vcf_files, "output/", is_nanopore=False)
   
   # Analyze all fragment files
   results = {}
   for bam_file in bam_files:
       sample_name = os.path.splitext(os.path.basename(bam_file))[0]
       frag_file = FragFile(f"output/{sample_name}.frag")
       fragments = frag_file.get_fragment_list()
       
       results[sample_name] = {
           'total_fragments': fragments.length(),
           'bogus_fragments': fragments.count_bogus_fragments(),
           'mutated_fragments': fragments.count_mutated_fragments()
       }
       
       frag_file.close()
   
   # Summary report
   print("Sample Summary:")
   for sample, stats in results.items():
       print(f"  {sample}:")
       print(f"    Total: {stats['total_fragments']}")
       print(f"    Bogus: {stats['bogus_fragments']}")
       print(f"    Mutated: {stats['mutated_fragments']}")

Simulation Example
------------------

Generate and analyze synthetic data:

.. code-block:: python

   from pyfraglib import FragmentSimulator, TissueMixtureSimulator
   
   # Basic simulation
   simulator = FragmentSimulator(
       fasta_path="reference.fasta",
       output_name="synthetic_basic"
   )
   
   regions = [("chr1", 1000000, 2000000)]
   fragments = simulator.simulate_fragments(regions, n_fragments=10000)
   
   # Tissue mixture simulation
   mixture_simulator = TissueMixtureSimulator(
       fasta_path="reference.fasta",
       output_name="synthetic_mixture"
   )
   
   tissue_fractions = {
       "hematopoietic": 0.7,
       "liver": 0.2,
       "tumor": 0.1
   }
   
   mixture_fragments = mixture_simulator.simulate_mixture(
       regions, tissue_fractions, n_fragments=50000
   )
   
   # Compare fragment characteristics
   print("Basic simulation:")
   print(f"  Mean length: {np.mean([f.length for f in fragments]):.1f}")
   
   print("Mixture simulation:")
   print(f"  Mean length: {np.mean([f.length for f in mixture_fragments]):.1f}")

Advanced Analysis
-----------------

Combine multiple analysis types:

.. code-block:: python

   from pyfraglib import Fragment, fit_gmm, fragments_per_chromosome_barplot
   from pyfraglib.scores import windowed_protection_score
   import numpy as np
   
   # Load data
   fragments = Fragment.from_bam("sample.bam", "variants.vcf")
   
   # 1. Length distribution analysis
   lengths = [f.length for f in fragments if not f.is_bogus]
   gmm_result = fit_gmm(lengths, n_components=3)
   
   # Extract component statistics
   means = gmm_result.means_.flatten()
   weights = gmm_result.weights_
   
   print("GMM Components:")
   for i, (mean, weight) in enumerate(zip(means, weights)):
       print(f"  Component {i+1}: mean={mean:.1f}, weight={weight:.3f}")
   
   # 2. Chromosome distribution
   fragments_per_chromosome_barplot(fragments, "chr_distribution.png")
   
   # 3. Fragment characteristics by chromosome
   chr_stats = {}
   for fragment in fragments:
       if fragment.is_bogus:
           continue
       
       chr_name = fragment.chrom
       if chr_name not in chr_stats:
           chr_stats[chr_name] = []
       chr_stats[chr_name].append(fragment.length)
   
   print("Length statistics by chromosome:")
   for chr_name, lengths in sorted(chr_stats.items()):
       mean_len = np.mean(lengths)
       std_len = np.std(lengths)
       print(f"  {chr_name}: mean={mean_len:.1f}, std={std_len:.1f}")
   
   # 4. End motif analysis
   motifs_5p, motifs_3p, num_frags = fragments.count_endmotifs(4)
   
   print(f"Top 5' motifs:")
   for motif, count in sorted(motifs_5p.items(), key=lambda x: x[1], reverse=True)[:5]:
       print(f"  {motif}: {count} ({count/num_frags:.3f})")

Integration with External Tools
-------------------------------

Using pyfraglib with other bioinformatics tools:

.. code-block:: python

   import pysam
   from pyfraglib import Fragment
   
   # Extract specific regions using pysam
   bam = pysam.AlignmentFile("sample.bam", "rb")
   vcf = pysam.VariantFile("variants.vcf")
   
   # Get reads from specific region
   region_reads = list(bam.fetch("chr1", 1000000, 2000000))
   print(f"Reads in region: {len(region_reads)}")
   
   # Process with pyfraglib
   fragments = Fragment.from_bam("sample.bam", "variants.vcf")
   
   # Filter fragments to same region
   region_fragments = [
       f for f in fragments 
       if f.chrom == "chr1" and f.start_pos >= 1000000 and f.end_pos <= 2000000
   ]
   
   print(f"Fragments in region: {len(region_fragments)}")
   
   bam.close()
   vcf.close()

These examples demonstrate the flexibility and power of pyfraglib for comprehensive fragmentomics analysis. For more specific use cases, refer to the API documentation and the notebooks in the repository.