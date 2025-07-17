Batch Processing
================

Efficiently process multiple samples using parallel processing and batch operations:

.. code-block:: python

   from pyfraglib import Fragment, FragFile
   from pyfraglib.stats import fragments_per_chromosome_barplot
   import os
   import pandas as pd
   
   # Process multiple BAM files in parallel
   bam_files = ["sample1.bam", "sample2.bam", "sample3.bam"]
   vcf_files = ["sample1.vcf", "sample2.vcf", "sample3.vcf"]
   
   # Use parallel processing to convert BAM to FRAG files
   Fragment.bams_to_frags(bam_files, vcf_files, "output/", is_nanopore=False)
   
   # Analyze all generated fragment files
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
   
   # Print summary
   print("Batch Processing Results:")
   for sample, stats in results.items():
       mutation_rate = stats['mutated_fragments'] / stats['total_fragments'] * 100
       print(f"  {sample}:")
       print(f"    Total fragments: {stats['total_fragments']:,}")
       print(f"    Bogus fragments: {stats['bogus_fragments']:,}")
       print(f"    Mutated fragments: {stats['mutated_fragments']:,}")
       print(f"    Mutation rate: {mutation_rate:.2f}%")

Large-Scale Cohort Analysis
---------------------------

Process and analyze large cohorts efficiently:

.. code-block:: python

   from pyfraglib import Fragment, FragFile
   from pyfraglib.scores import motif_diversity
   from pyfraglib.lengths import fragment_length_gmm
   import pandas as pd
   import numpy as np
   import matplotlib.pyplot as plt
   import seaborn as sns
   from pathlib import Path
   
   # Define cohort structure
   cohort_data = {
       'healthy': [
           ('healthy_01', 'healthy_01.bam', 'healthy_01.vcf'),
           ('healthy_02', 'healthy_02.bam', 'healthy_02.vcf'),
           ('healthy_03', 'healthy_03.bam', 'healthy_03.vcf'),
       ],
       'cancer': [
           ('cancer_01', 'cancer_01.bam', 'cancer_01.vcf'),
           ('cancer_02', 'cancer_02.bam', 'cancer_02.vcf'),
           ('cancer_03', 'cancer_03.bam', 'cancer_03.vcf'),
       ],
       'treated': [
           ('treated_01', 'treated_01.bam', 'treated_01.vcf'),
           ('treated_02', 'treated_02.bam', 'treated_02.vcf'),
           ('treated_03', 'treated_03.bam', 'treated_03.vcf'),
       ]
   }
   
   # Create output directory structure
   output_dir = Path("cohort_analysis")
   output_dir.mkdir(exist_ok=True)
   
   # Process each group
   cohort_results = []
   
   for group, samples in cohort_data.items():
       print(f"\nProcessing {group} group...")
       
       # Extract BAM and VCF files
       bam_files = [sample[1] for sample in samples]
       vcf_files = [sample[2] for sample in samples]
       
       # Parallel processing
       Fragment.bams_to_frags(bam_files, vcf_files, f"{output_dir}/{group}/", is_nanopore=False)
       
       # Analyze each sample in the group
       for sample_name, bam_file, vcf_file in samples:
           frag_file = FragFile(f"{output_dir}/{group}/{sample_name}.frag")
           fragments = frag_file.get_fragment_list()
           
           # Calculate comprehensive statistics
           total_fragments = fragments.length()
           bogus_fragments = fragments.count_bogus_fragments()
           mutated_fragments = fragments.count_mutated_fragments()
           
           # Fragment length statistics
           lengths = [f.length for f in fragments if not f.is_bogus]
           mean_length = np.mean(lengths)
           std_length = np.std(lengths)
           median_length = np.median(lengths)
           
           # Motif diversity
           shannon_diversity = motif_diversity(fragments, kmer_len=4, index="shannon")
           simpson_diversity = motif_diversity(fragments, kmer_len=4, index="simpson")
           
           # Store results
           cohort_results.append({
               'sample': sample_name,
               'group': group,
               'total_fragments': total_fragments,
               'bogus_fragments': bogus_fragments,
               'mutated_fragments': mutated_fragments,
               'mutation_rate': mutated_fragments / total_fragments if total_fragments > 0 else 0,
               'mean_length': mean_length,
               'std_length': std_length,
               'median_length': median_length,
               'shannon_diversity': shannon_diversity,
               'simpson_diversity': simpson_diversity
           })
           
           frag_file.close()
           print(f"  {sample_name}: {total_fragments:,} fragments processed")
   
   # Create results dataframe
   df = pd.DataFrame(cohort_results)
   
   # Save results
   df.to_csv(f"{output_dir}/cohort_analysis_results.csv", index=False)
   print(f"\nResults saved to {output_dir}/cohort_analysis_results.csv")

Statistical Analysis and Visualization
--------------------------------------

Perform comprehensive statistical analysis on batch results:

.. code-block:: python

   import pandas as pd
   import numpy as np
   import matplotlib.pyplot as plt
   import seaborn as sns
   from scipy import stats
   
   # Load results (assuming df from previous example)
   # df = pd.read_csv("cohort_analysis/cohort_analysis_results.csv")
   
   # Group statistics
   group_stats = df.groupby('group').agg({
       'total_fragments': ['mean', 'std', 'min', 'max'],
       'mutation_rate': ['mean', 'std', 'min', 'max'],
       'mean_length': ['mean', 'std', 'min', 'max'],
       'shannon_diversity': ['mean', 'std', 'min', 'max'],
       'simpson_diversity': ['mean', 'std', 'min', 'max']
   }).round(4)
   
   print("Group Statistics:")
   print(group_stats)
   
   # Statistical tests between groups
   groups = df['group'].unique()
   print(f"\nStatistical comparisons between groups:")
   
   for i, group1 in enumerate(groups):
       for j, group2 in enumerate(groups):
           if i < j:  # Avoid duplicate comparisons
               data1 = df[df['group'] == group1]
               data2 = df[df['group'] == group2]
               
               # T-test for mutation rate
               t_stat, p_val = stats.ttest_ind(data1['mutation_rate'], data2['mutation_rate'])
               print(f"  {group1} vs {group2} (mutation rate): t={t_stat:.3f}, p={p_val:.3e}")
               
               # T-test for mean length
               t_stat, p_val = stats.ttest_ind(data1['mean_length'], data2['mean_length'])
               print(f"  {group1} vs {group2} (mean length): t={t_stat:.3f}, p={p_val:.3e}")
               
               # T-test for Shannon diversity
               t_stat, p_val = stats.ttest_ind(data1['shannon_diversity'], data2['shannon_diversity'])
               print(f"  {group1} vs {group2} (Shannon diversity): t={t_stat:.3f}, p={p_val:.3e}")
   
   # Create comprehensive visualization
   fig, axes = plt.subplots(2, 3, figsize=(18, 12))
   
   # Mutation rate boxplot
   sns.boxplot(data=df, x='group', y='mutation_rate', ax=axes[0, 0])
   axes[0, 0].set_title('Mutation Rate by Group')
   axes[0, 0].set_ylabel('Mutation Rate')
   
   # Mean length boxplot
   sns.boxplot(data=df, x='group', y='mean_length', ax=axes[0, 1])
   axes[0, 1].set_title('Mean Fragment Length by Group')
   axes[0, 1].set_ylabel('Mean Length (bp)')
   
   # Shannon diversity boxplot
   sns.boxplot(data=df, x='group', y='shannon_diversity', ax=axes[0, 2])
   axes[0, 2].set_title('Shannon Diversity by Group')
   axes[0, 2].set_ylabel('Shannon Entropy')
   
   # Total fragments barplot
   sns.barplot(data=df, x='group', y='total_fragments', ax=axes[1, 0])
   axes[1, 0].set_title('Total Fragments by Group')
   axes[1, 0].set_ylabel('Fragment Count')
   
   # Correlation matrix
   corr_data = df[['mutation_rate', 'mean_length', 'shannon_diversity', 'simpson_diversity']]
   corr_matrix = corr_data.corr()
   sns.heatmap(corr_matrix, annot=True, cmap='coolwarm', center=0, ax=axes[1, 1])
   axes[1, 1].set_title('Correlation Matrix')
   
   # PCA plot
   from sklearn.decomposition import PCA
   from sklearn.preprocessing import StandardScaler
   
   features = ['mutation_rate', 'mean_length', 'shannon_diversity', 'simpson_diversity']
   X = df[features].values
   
   scaler = StandardScaler()
   X_scaled = scaler.fit_transform(X)
   
   pca = PCA(n_components=2)
   X_pca = pca.fit_transform(X_scaled)
   
   for group in groups:
       mask = df['group'] == group
       axes[1, 2].scatter(X_pca[mask, 0], X_pca[mask, 1], label=group, alpha=0.7)
   
   axes[1, 2].set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.2%} variance)')
   axes[1, 2].set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.2%} variance)')
   axes[1, 2].set_title('PCA of Fragmentomics Features')
   axes[1, 2].legend()
   
   plt.tight_layout()
   plt.savefig('cohort_analysis_visualization.png', dpi=300, bbox_inches='tight')
   plt.show()

Automated Report Generation
---------------------------

Generate comprehensive reports for batch processing results:

.. code-block:: python

   import pandas as pd
   import matplotlib.pyplot as plt
   import seaborn as sns
   from pathlib import Path
   
   def generate_cohort_report(results_df, output_dir):
       """Generate comprehensive cohort analysis report"""
       
       output_path = Path(output_dir)
       output_path.mkdir(exist_ok=True)
       
       # Summary statistics
       summary_stats = results_df.groupby('group').agg({
           'total_fragments': ['count', 'mean', 'std'],
           'mutation_rate': ['mean', 'std'],
           'mean_length': ['mean', 'std'],
           'shannon_diversity': ['mean', 'std']
       }).round(4)
       
       # Save summary table
       summary_stats.to_csv(f"{output_path}/summary_statistics.csv")
       
       # Generate plots
       fig, axes = plt.subplots(2, 2, figsize=(15, 12))
       
       # Plot 1: Mutation rate distribution
       sns.violinplot(data=results_df, x='group', y='mutation_rate', ax=axes[0, 0])
       axes[0, 0].set_title('Mutation Rate Distribution')
       axes[0, 0].set_ylabel('Mutation Rate')
       
       # Plot 2: Fragment length distribution
       sns.violinplot(data=results_df, x='group', y='mean_length', ax=axes[0, 1])
       axes[0, 1].set_title('Fragment Length Distribution')
       axes[0, 1].set_ylabel('Mean Length (bp)')
       
       # Plot 3: Diversity comparison
       sns.scatterplot(data=results_df, x='shannon_diversity', y='simpson_diversity', 
                      hue='group', ax=axes[1, 0])
       axes[1, 0].set_title('Motif Diversity Comparison')
       axes[1, 0].set_xlabel('Shannon Entropy')
       axes[1, 0].set_ylabel('Simpson Index')
       
       # Plot 4: Sample overview
       sample_metrics = results_df.melt(id_vars=['sample', 'group'], 
                                       value_vars=['mutation_rate', 'mean_length', 'shannon_diversity'],
                                       var_name='metric', value_name='value')
       
       sns.boxplot(data=sample_metrics, x='metric', y='value', hue='group', ax=axes[1, 1])
       axes[1, 1].set_title('Multi-metric Comparison')
       axes[1, 1].set_xlabel('Metric')
       axes[1, 1].set_ylabel('Value')
       axes[1, 1].tick_params(axis='x', rotation=45)
       
       plt.tight_layout()
       plt.savefig(f"{output_path}/cohort_analysis_report.png", dpi=300, bbox_inches='tight')
       plt.close()
       
       # Generate HTML report
       html_content = f"""
       <!DOCTYPE html>
       <html>
       <head>
           <title>Cohort Analysis Report</title>
           <style>
               body {{ font-family: Arial, sans-serif; margin: 40px; }}
               table {{ border-collapse: collapse; width: 100%; }}
               th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
               th {{ background-color: #f2f2f2; }}
               .section {{ margin: 20px 0; }}
           </style>
       </head>
       <body>
           <h1>Cohort Analysis Report</h1>
           
           <div class="section">
               <h2>Summary Statistics</h2>
               {summary_stats.to_html()}
           </div>
           
           <div class="section">
               <h2>Detailed Results</h2>
               {results_df.to_html(index=False)}
           </div>
           
           <div class="section">
               <h2>Visualization</h2>
               <img src="cohort_analysis_report.png" alt="Cohort Analysis Plots" style="max-width: 100%;">
           </div>
           
           <div class="section">
               <h2>Key Findings</h2>
               <ul>
                   <li>Total samples analyzed: {len(results_df)}</li>
                   <li>Groups compared: {', '.join(results_df['group'].unique())}</li>
                   <li>Average fragments per sample: {results_df['total_fragments'].mean():.0f}</li>
                   <li>Overall mutation rate: {results_df['mutation_rate'].mean():.3f}</li>
               </ul>
           </div>
       </body>
       </html>
       """
       
       with open(f"{output_path}/cohort_report.html", "w") as f:
           f.write(html_content)
       
       print(f"Report generated in {output_path}/")
       print("Files created:")
       print(f"  - summary_statistics.csv")
       print(f"  - cohort_analysis_report.png")
       print(f"  - cohort_report.html")
   
   # Generate report (using df from previous example)
   generate_cohort_report(df, "cohort_analysis/report")

Distributed Processing
----------------------

Scale to very large datasets using distributed processing:

.. code-block:: python

   from pyfraglib import Fragment
   from pathlib import Path
   import multiprocessing as mp
   from functools import partial
   
   def process_sample_batch(sample_batch, output_dir):
       """Process a batch of samples"""
       results = []
       
       for sample_name, bam_file, vcf_file in sample_batch:
           try:
               # Process individual sample
               fragments = Fragment.from_bam(bam_file, vcf_file)
               
               # Save to fragment file
               fragments.to_frag_file(sample_name, output_dir)
               
               # Calculate basic statistics
               total_fragments = fragments.length()
               bogus_fragments = fragments.count_bogus_fragments()
               mutated_fragments = fragments.count_mutated_fragments()
               
               results.append({
                   'sample': sample_name,
                   'status': 'success',
                   'total_fragments': total_fragments,
                   'bogus_fragments': bogus_fragments,
                   'mutated_fragments': mutated_fragments
               })
               
           except Exception as e:
               results.append({
                   'sample': sample_name,
                   'status': 'error',
                   'error': str(e)
               })
       
       return results
   
   def process_large_cohort(sample_list, output_dir, batch_size=10, n_processes=None):
       """Process large cohort with batching and multiprocessing"""
       
       if n_processes is None:
           n_processes = mp.cpu_count()
       
       # Create output directory
       Path(output_dir).mkdir(parents=True, exist_ok=True)
       
       # Split samples into batches
       batches = []
       for i in range(0, len(sample_list), batch_size):
           batch = sample_list[i:i + batch_size]
           batches.append(batch)
       
       print(f"Processing {len(sample_list)} samples in {len(batches)} batches using {n_processes} processes")
       
       # Process batches in parallel
       with mp.Pool(processes=n_processes) as pool:
           process_func = partial(process_sample_batch, output_dir=output_dir)
           batch_results = pool.map(process_func, batches)
       
       # Combine results
       all_results = []
       for batch_result in batch_results:
           all_results.extend(batch_result)
       
       # Summary
       successful = sum(1 for r in all_results if r['status'] == 'success')
       failed = sum(1 for r in all_results if r['status'] == 'error')
       
       print(f"\nProcessing complete:")
       print(f"  Successful: {successful}")
       print(f"  Failed: {failed}")
       
       return all_results
   
   # Example usage
   large_sample_list = [
       (f"sample_{i:03d}", f"data/sample_{i:03d}.bam", f"data/sample_{i:03d}.vcf")
       for i in range(1, 1001)  # 1000 samples
   ]
   
   # Process large cohort
   results = process_large_cohort(large_sample_list, "large_cohort_output/", 
                                 batch_size=50, n_processes=8)

Performance Monitoring
-----------------------

Monitor processing performance and resource usage:

.. code-block:: python

   import time
   import psutil
   from pyfraglib import Fragment
   
   def monitor_processing_performance(bam_files, vcf_files, output_dir):
       """Monitor performance during batch processing"""
       
       start_time = time.time()
       process = psutil.Process()
       
       # Get initial resource usage
       initial_memory = process.memory_info().rss / 1024 / 1024  # MB
       initial_cpu = process.cpu_percent()
       
       print(f"Starting batch processing at {time.strftime('%Y-%m-%d %H:%M:%S')}")
       print(f"Initial memory usage: {initial_memory:.1f} MB")
       print(f"Processing {len(bam_files)} files...")
       
       # Process files
       Fragment.bams_to_frags(bam_files, vcf_files, output_dir, is_nanopore=False)
       
       # Get final resource usage
       end_time = time.time()
       final_memory = process.memory_info().rss / 1024 / 1024  # MB
       final_cpu = process.cpu_percent()
       
       processing_time = end_time - start_time
       memory_increase = final_memory - initial_memory
       
       print(f"\nProcessing completed at {time.strftime('%Y-%m-%d %H:%M:%S')}")
       print(f"Total processing time: {processing_time:.1f} seconds")
       print(f"Average time per file: {processing_time / len(bam_files):.1f} seconds")
       print(f"Final memory usage: {final_memory:.1f} MB")
       print(f"Memory increase: {memory_increase:.1f} MB")
       print(f"Peak CPU usage: {final_cpu:.1f}%")
       
       return {
           'total_time': processing_time,
           'time_per_file': processing_time / len(bam_files),
           'memory_usage': final_memory,
           'memory_increase': memory_increase,
           'cpu_usage': final_cpu
       }

See Also
--------

* :doc:`basic_workflow` - Basic fragmentomics workflow
* :doc:`mutation_analysis` - Mutation analysis techniques
* :doc:`fragmentomics_scores` - Computing fragmentomics scores
* :doc:`simulation_examples` - Simulation examples
* :func:`pyfraglib.Fragment.bams_to_frags` - Parallel BAM processing
* :func:`pyfraglib.Fragment.from_bams` - Multi-sample fragment extraction