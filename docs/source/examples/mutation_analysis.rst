Mutation Analysis
=================

Analyze fragment mutation patterns using VCF files for variant annotation:

.. code-block:: python

   from pyfraglib import Fragment, FragFile
   import pandas as pd
   import matplotlib.pyplot as plt
   
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

Advanced Mutation Analysis
--------------------------

Compare fragment lengths between mutated and wildtype fragments:

.. code-block:: python

   from pyfraglib import Fragment, FragFile
   import numpy as np
   import matplotlib.pyplot as plt
   import seaborn as sns
   
   # Load fragments
   fragments = Fragment.from_bam("tumor.bam", "somatic_variants.vcf")
   
   # Separate mutated and wildtype fragments
   mutated_lengths = []
   wildtype_lengths = []
   
   for fragment in fragments:
       if fragment.is_bogus:
           continue
       
       if fragment.is_mutated:
           mutated_lengths.append(fragment.length)
       else:
           wildtype_lengths.append(fragment.length)
   
   # Statistical comparison
   print(f"Mutated fragments: {len(mutated_lengths)}")
   print(f"Wildtype fragments: {len(wildtype_lengths)}")
   print(f"Mean length (mutated): {np.mean(mutated_lengths):.1f}")
   print(f"Mean length (wildtype): {np.mean(wildtype_lengths):.1f}")
   
   # Visualization
   plt.figure(figsize=(10, 6))
   plt.hist(wildtype_lengths, bins=50, alpha=0.7, label='Wildtype', density=True)
   plt.hist(mutated_lengths, bins=50, alpha=0.7, label='Mutated', density=True)
   plt.xlabel('Fragment Length (bp)')
   plt.ylabel('Density')
   plt.legend()
   plt.title('Fragment Length Distribution: Mutated vs Wildtype')
   plt.savefig('mutation_length_comparison.png', dpi=300, bbox_inches='tight')
   plt.show()

Mutation Context Analysis
-------------------------

Analyze the sequence context around mutations:

.. code-block:: python

   from pyfraglib import Fragment
   from collections import defaultdict
   
   # Extract fragments with mutation annotation
   fragments = Fragment.from_bam("tumor.bam", "somatic_variants.vcf")
   
   # Analyze end motifs for mutated vs wildtype fragments
   mutated_5p_motifs = defaultdict(int)
   wildtype_5p_motifs = defaultdict(int)
   
   for fragment in fragments:
       if fragment.is_bogus:
           continue
       
       motif_5p = fragment.end5p[:4]  # First 4 bases of 5' end
       
       if fragment.is_mutated:
           mutated_5p_motifs[motif_5p] += 1
       else:
           wildtype_5p_motifs[motif_5p] += 1
   
   # Compare top motifs
   print("Top 5' motifs in mutated fragments:")
   for motif, count in sorted(mutated_5p_motifs.items(), key=lambda x: x[1], reverse=True)[:10]:
       print(f"  {motif}: {count}")
   
   print("\nTop 5' motifs in wildtype fragments:")
   for motif, count in sorted(wildtype_5p_motifs.items(), key=lambda x: x[1], reverse=True)[:10]:
       print(f"  {motif}: {count}")

Working with Multiple Samples
------------------------------

Compare mutation patterns across different samples:

.. code-block:: python

   from pyfraglib import Fragment
   import pandas as pd
   
   # Define samples
   samples = [
       ("patient1_baseline", "patient1_baseline.bam", "patient1_baseline.vcf"),
       ("patient1_treatment", "patient1_treatment.bam", "patient1_treatment.vcf"),
       ("patient2_baseline", "patient2_baseline.bam", "patient2_baseline.vcf"),
       ("patient2_treatment", "patient2_treatment.bam", "patient2_treatment.vcf"),
   ]
   
   # Analyze each sample
   results = []
   for sample_name, bam_file, vcf_file in samples:
       fragments = Fragment.from_bam(bam_file, vcf_file)
       
       total = fragments.length()
       mutated = fragments.count_mutated_fragments()
       mutation_rate = mutated / total if total > 0 else 0
       
       results.append({
           'sample': sample_name,
           'total_fragments': total,
           'mutated_fragments': mutated,
           'mutation_rate': mutation_rate
       })
   
   # Create summary dataframe
   df = pd.DataFrame(results)
   print(df)
   
   # Visualize mutation rates
   import matplotlib.pyplot as plt
   
   plt.figure(figsize=(10, 6))
   bars = plt.bar(df['sample'], df['mutation_rate'])
   plt.xlabel('Sample')
   plt.ylabel('Mutation Rate')
   plt.title('Mutation Rates Across Samples')
   plt.xticks(rotation=45)
   
   # Add value labels on bars
   for bar, rate in zip(bars, df['mutation_rate']):
       plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.001, 
                f'{rate:.3f}', ha='center', va='bottom')
   
   plt.tight_layout()
   plt.savefig('mutation_rates_comparison.png', dpi=300, bbox_inches='tight')
   plt.show()

Important Notes
---------------

* VCF files must be properly formatted and indexed
* Only SNVs (single nucleotide variants) are currently supported
* Mutation annotation requires that reads overlap variant positions
* Use :func:`pyfraglib.Fragment.build_mutated_reads_set` to understand mutation detection process

See Also
--------

* :doc:`basic_workflow` - Basic fragmentomics workflow
* :doc:`fragmentomics_scores` - Computing fragmentomics scores
* :doc:`batch_processing` - Processing multiple samples