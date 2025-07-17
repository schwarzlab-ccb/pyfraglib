Fragmentomics Scores
====================

Calculate windowed protection scores (WPS) and motif diversity for fragmentomics analysis:

.. code-block:: python

   from pyfraglib import Fragment, FragFile
   from pyfraglib.scores import windowed_protection_score, motif_diversity
   import pandas as pd
   import matplotlib.pyplot as plt
   
   # Load fragments
   fragments = Fragment.from_bam("sample.bam", "variants.vcf")
   
   # Calculate motif diversity
   shannon_diversity = motif_diversity(fragments, kmer_len=4, index="shannon")
   simpson_diversity = motif_diversity(fragments, kmer_len=4, index="simpson")
   
   print(f"Shannon entropy: {shannon_diversity:.3f}")
   print(f"Simpson index: {simpson_diversity:.3f}")

Windowed Protection Score (WPS)
-------------------------------

Calculate WPS for genomic regions of interest:

.. code-block:: python

   from pyfraglib import Fragment
   from pyfraglib.scores import windowed_protection_score
   import pysam
   
   # Load fragments
   fragments = Fragment.from_bam("sample.bam", "variants.vcf")
   
   # Create BED file with regions of interest
   bed_content = """chr1	1000000	1001000	gene1
   chr1	2000000	2001000	gene2
   chr2	3000000	3001000	gene3"""
   
   with open("regions.bed", "w") as f:
       f.write(bed_content)
   
   # Create tabix file
   regions = pysam.TabixFile("regions.bed")
   
   # Calculate WPS scores
   wps_df = windowed_protection_score(fragments, regions, win_size=120, genome="hg19")
   
   print(wps_df.head())
   
   # Plot WPS scores
   plt.figure(figsize=(12, 6))
   plt.plot(wps_df['abs_pos'], wps_df['wps'])
   plt.xlabel('Absolute Genomic Position')
   plt.ylabel('WPS Score')
   plt.title('Windowed Protection Score')
   plt.savefig('wps_scores.png', dpi=300, bbox_inches='tight')
   plt.show()

Advanced WPS Analysis
---------------------

Analyze WPS patterns across different genomic regions:

.. code-block:: python

   from pyfraglib import Fragment
   from pyfraglib.scores import windowed_protection_score, score_line_plot
   import pysam
   import pandas as pd
   
   # Load fragments
   fragments = Fragment.from_bam("sample.bam", "variants.vcf")
   
   # Create regions around TSS sites
   tss_regions = pd.DataFrame({
       'chrom': ['chr1', 'chr1', 'chr2', 'chr2'],
       'start': [1000000, 2000000, 3000000, 4000000],
       'end': [1001000, 2001000, 3001000, 4001000],
       'gene': ['GENE1', 'GENE2', 'GENE3', 'GENE4']
   })
   
   # Save as BED file
   tss_regions.to_csv("tss_regions.bed", sep='\t', index=False, header=False)
   
   # Calculate WPS
   regions = pysam.TabixFile("tss_regions.bed")
   wps_df = windowed_protection_score(fragments, regions, win_size=120)
   
   # Generate WPS line plot
   score_line_plot(wps_df, "sample", "output/", score="wps", genome="hg19")
   
   # Analyze WPS statistics by region
   region_stats = wps_df.groupby('chrom').agg({
       'wps': ['mean', 'std', 'min', 'max'],
       'depth': ['mean', 'std']
   }).round(3)
   
   print("WPS statistics by chromosome:")
   print(region_stats)

Motif Diversity Analysis
------------------------

Detailed analysis of fragment end motifs:

.. code-block:: python

   from pyfraglib import Fragment
   from pyfraglib.scores import motif_diversity
   from pyfraglib.core import shannon_entropy, simpson_index
   import matplotlib.pyplot as plt
   import seaborn as sns
   from collections import defaultdict
   
   # Load fragments
   fragments = Fragment.from_bam("sample.bam", "variants.vcf")
   
   # Calculate diversity for different k-mer lengths
   kmer_lengths = [2, 3, 4, 5, 6]
   shannon_scores = []
   simpson_scores = []
   
   for k in kmer_lengths:
       shannon = motif_diversity(fragments, kmer_len=k, index="shannon")
       simpson = motif_diversity(fragments, kmer_len=k, index="simpson")
       
       shannon_scores.append(shannon)
       simpson_scores.append(simpson)
   
   # Plot diversity vs k-mer length
   fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
   
   ax1.plot(kmer_lengths, shannon_scores, 'o-', label='Shannon entropy')
   ax1.set_xlabel('K-mer length')
   ax1.set_ylabel('Shannon entropy')
   ax1.set_title('Shannon Entropy vs K-mer Length')
   ax1.grid(True)
   
   ax2.plot(kmer_lengths, simpson_scores, 'o-', label='Simpson index', color='orange')
   ax2.set_xlabel('K-mer length')
   ax2.set_ylabel('Simpson index')
   ax2.set_title('Simpson Index vs K-mer Length')
   ax2.grid(True)
   
   plt.tight_layout()
   plt.savefig('motif_diversity_kmer.png', dpi=300, bbox_inches='tight')
   plt.show()

Comparative Analysis
--------------------

Compare fragmentomics scores between different samples:

.. code-block:: python

   from pyfraglib import Fragment
   from pyfraglib.scores import motif_diversity, windowed_protection_score
   import pandas as pd
   import matplotlib.pyplot as plt
   import seaborn as sns
   
   # Define samples
   samples = [
       ("healthy", "healthy.bam", "healthy.vcf"),
       ("cancer", "cancer.bam", "cancer.vcf"),
       ("treated", "treated.bam", "treated.vcf")
   ]
   
   # Calculate scores for each sample
   results = []
   for sample_name, bam_file, vcf_file in samples:
       fragments = Fragment.from_bam(bam_file, vcf_file)
       
       # Calculate motif diversity
       shannon = motif_diversity(fragments, kmer_len=4, index="shannon")
       simpson = motif_diversity(fragments, kmer_len=4, index="simpson")
       
       results.append({
           'sample': sample_name,
           'shannon_entropy': shannon,
           'simpson_index': simpson,
           'total_fragments': fragments.length(),
           'mutation_rate': fragments.count_mutated_fragments() / fragments.length()
       })
   
   # Create results dataframe
   df = pd.DataFrame(results)
   print(df)
   
   # Create comparison plots
   fig, axes = plt.subplots(2, 2, figsize=(12, 10))
   
   # Shannon entropy
   axes[0, 0].bar(df['sample'], df['shannon_entropy'])
   axes[0, 0].set_title('Shannon Entropy')
   axes[0, 0].set_ylabel('Entropy')
   
   # Simpson index
   axes[0, 1].bar(df['sample'], df['simpson_index'])
   axes[0, 1].set_title('Simpson Index')
   axes[0, 1].set_ylabel('Index')
   
   # Total fragments
   axes[1, 0].bar(df['sample'], df['total_fragments'])
   axes[1, 0].set_title('Total Fragments')
   axes[1, 0].set_ylabel('Count')
   
   # Mutation rate
   axes[1, 1].bar(df['sample'], df['mutation_rate'])
   axes[1, 1].set_title('Mutation Rate')
   axes[1, 1].set_ylabel('Rate')
   
   plt.tight_layout()
   plt.savefig('scores_comparison.png', dpi=300, bbox_inches='tight')
   plt.show()

Performance Considerations
--------------------------

For large datasets, consider these optimization strategies:

.. code-block:: python

   from pyfraglib import Fragment
   from pyfraglib.scores import windowed_protection_score_fast
   import pysam
   
   # Use the fast WPS implementation for large datasets
   fragments = Fragment.from_bam("large_sample.bam", "variants.vcf")
   regions = pysam.TabixFile("regions.bed")
   
   # Fast WPS calculation
   wps_df = windowed_protection_score_fast(fragments, regions, win_size=120)
   
   # For very large datasets, consider processing in chunks
   from pyfraglib.fragment import FragmentCollection
   
   # Process multiple BAM files in parallel
   bam_files = ["sample1.bam", "sample2.bam", "sample3.bam"]
   vcf_files = ["sample1.vcf", "sample2.vcf", "sample3.vcf"]
   
   Fragment.bams_to_frags(bam_files, vcf_files, "output/", is_nanopore=False)

See Also
--------

* :doc:`basic_workflow` - Basic fragmentomics workflow
* :doc:`mutation_analysis` - Mutation analysis techniques
* :doc:`length_analysis` - Fragment length analysis
* :func:`pyfraglib.scores.windowed_protection_score` - WPS calculation
* :func:`pyfraglib.scores.motif_diversity` - Motif diversity calculation