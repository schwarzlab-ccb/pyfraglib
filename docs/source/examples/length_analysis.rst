Fragment Length Analysis
========================

Analyze fragment length distributions using Gaussian mixture models and statistical methods:

.. code-block:: python

   from pyfraglib import Fragment, FragFile
   from pyfraglib.math import fit_gmm, plot_gmm, goodness_of_fit_stats
   from pyfraglib.lengths import fragment_length_plot, fragment_length_gmm
   import numpy as np
   import matplotlib.pyplot as plt
   
   # Load fragments
   fragments = Fragment.from_bam("sample.bam", "variants.vcf")
   
   # Extract fragment lengths (exclude bogus fragments)
   lengths = np.array([f.length for f in fragments if not f.is_bogus])
   
   # Basic statistics
   print(f"Total fragments: {len(lengths)}")
   print(f"Mean length: {np.mean(lengths):.1f} bp")
   print(f"Median length: {np.median(lengths):.1f} bp")
   print(f"Standard deviation: {np.std(lengths):.1f} bp")
   print(f"Length range: {np.min(lengths)} - {np.max(lengths)} bp")

Gaussian Mixture Model Analysis
-------------------------------

Fit Gaussian mixture models to identify fragment populations:

.. code-block:: python

   from pyfraglib.math import fit_gmm, plot_gmm, goodness_of_fit_stats
   import numpy as np
   
   # Load fragments
   fragments = Fragment.from_bam("sample.bam", "variants.vcf")
   lengths = np.array([f.length for f in fragments if not f.is_bogus], dtype=np.float64)
   
   # Create GMM configuration
   gmm_config = {
       "number_of_gaussians": 3,
       "subsample_percentage": 0.1,
       "means_lower_bounds": [50, 150, 250],
       "means_upper_bounds": [120, 200, 350],
       "std_lower_bounds": [10, 15, 20],
       "std_upper_bounds": [30, 40, 60],
       "initial_means": [80, 167, 300],
       "initial_stds": [20, 25, 40],
       "initial_weights": [0.4, 0.4, 0.2]
   }
   
   # Save configuration
   import json
   with open("gmm_config.json", "w") as f:
       json.dump(gmm_config, f, indent=2)
   
   # Fit GMM
   result, n_components, params, data_subset = fit_gmm(lengths, "gmm_config.json")
   
   # Extract results
   means = params[:n_components]
   stds = params[n_components:2*n_components]
   weights = params[2*n_components:]
   
   print("GMM Results:")
   for i in range(n_components):
       print(f"  Component {i+1}: mean={means[i]:.1f}, std={stds[i]:.1f}, weight={weights[i]:.3f}")
   
   # Calculate goodness of fit
   gof_stats = goodness_of_fit_stats(lengths, params, n_components)
   print(f"\nGoodness of fit:")
   print(f"  Kolmogorov-Smirnov statistic: {gof_stats['kolmogorov_smirnov_statistic']:.4f}")
   print(f"  Kolmogorov-Smirnov p-value: {gof_stats['kolmogorov_smirnov_p_value']:.4f}")
   print(f"  Wasserstein distance: {gof_stats['wasserstein_distance']:.4f}")
   print(f"  Jensen-Shannon divergence: {gof_stats['jensen_shannon_divergence']:.4f}")
   
   # Plot GMM
   plot_gmm(lengths, n_components, params, "output/", "sample_gmm")

High-Level Length Analysis
--------------------------

Use the built-in length analysis functions:

.. code-block:: python

   from pyfraglib import Fragment
   from pyfraglib.lengths import fragment_length_plot, fragment_length_gmm
   
   # Load fragments
   fragments = Fragment.from_bam("sample.bam", "variants.vcf")
   
   # Create length histogram plots
   fragment_length_plot(fragments, "output/", "sample")
   
   # Fit GMM with configuration file
   fragment_length_gmm(fragments, "gmm_config.json", "output/", "sample")

Comparative Length Analysis
---------------------------

Compare fragment lengths between different conditions:

.. code-block:: python

   from pyfraglib import Fragment
   import numpy as np
   import matplotlib.pyplot as plt
   import seaborn as sns
   from scipy import stats
   
   # Load different samples
   samples = [
       ("healthy", "healthy.bam", "healthy.vcf"),
       ("cancer", "cancer.bam", "cancer.vcf"),
       ("treated", "treated.bam", "treated.vcf")
   ]
   
   # Extract lengths for each sample
   sample_lengths = {}
   for sample_name, bam_file, vcf_file in samples:
       fragments = Fragment.from_bam(bam_file, vcf_file)
       lengths = [f.length for f in fragments if not f.is_bogus]
       sample_lengths[sample_name] = lengths
   
   # Statistical comparison
   print("Length statistics by sample:")
   for sample_name, lengths in sample_lengths.items():
       print(f"  {sample_name}: mean={np.mean(lengths):.1f}, std={np.std(lengths):.1f}")
   
   # Perform pairwise t-tests
   samples_list = list(sample_lengths.keys())
   for i in range(len(samples_list)):
       for j in range(i+1, len(samples_list)):
           sample1, sample2 = samples_list[i], samples_list[j]
           lengths1, lengths2 = sample_lengths[sample1], sample_lengths[sample2]
           
           # Perform t-test
           t_stat, p_value = stats.ttest_ind(lengths1, lengths2)
           print(f"  {sample1} vs {sample2}: t={t_stat:.3f}, p={p_value:.2e}")
   
   # Visualization
   fig, axes = plt.subplots(2, 2, figsize=(12, 10))
   
   # Box plot
   axes[0, 0].boxplot([sample_lengths[name] for name in samples_list], labels=samples_list)
   axes[0, 0].set_title('Fragment Length Distribution')
   axes[0, 0].set_ylabel('Length (bp)')
   
   # Violin plot
   data_for_violin = []
   labels_for_violin = []
   for sample_name, lengths in sample_lengths.items():
       data_for_violin.extend(lengths)
       labels_for_violin.extend([sample_name] * len(lengths))
   
   import pandas as pd
   df = pd.DataFrame({'length': data_for_violin, 'sample': labels_for_violin})
   sns.violinplot(data=df, x='sample', y='length', ax=axes[0, 1])
   axes[0, 1].set_title('Fragment Length Distribution (Violin)')
   
   # Histogram overlay
   for sample_name, lengths in sample_lengths.items():
       axes[1, 0].hist(lengths, bins=50, alpha=0.7, label=sample_name, density=True)
   axes[1, 0].set_xlabel('Fragment Length (bp)')
   axes[1, 0].set_ylabel('Density')
   axes[1, 0].set_title('Fragment Length Histograms')
   axes[1, 0].legend()
   
   # Cumulative distribution
   for sample_name, lengths in sample_lengths.items():
       sorted_lengths = np.sort(lengths)
       y = np.arange(1, len(sorted_lengths) + 1) / len(sorted_lengths)
       axes[1, 1].plot(sorted_lengths, y, label=sample_name)
   axes[1, 1].set_xlabel('Fragment Length (bp)')
   axes[1, 1].set_ylabel('Cumulative Probability')
   axes[1, 1].set_title('Cumulative Distribution')
   axes[1, 1].legend()
   
   plt.tight_layout()
   plt.savefig('length_comparison.png', dpi=300, bbox_inches='tight')
   plt.show()

Fragment Length vs Mutation Status
-----------------------------------

Analyze how fragment lengths differ between mutated and wildtype fragments:

.. code-block:: python

   from pyfraglib import Fragment
   import numpy as np
   import matplotlib.pyplot as plt
   from scipy import stats
   
   # Load fragments with mutation annotation
   fragments = Fragment.from_bam("tumor.bam", "somatic_variants.vcf")
   
   # Separate lengths by mutation status
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
   print(f"Mean length (mutated): {np.mean(mutated_lengths):.1f} bp")
   print(f"Mean length (wildtype): {np.mean(wildtype_lengths):.1f} bp")
   
   # Perform t-test
   t_stat, p_value = stats.ttest_ind(mutated_lengths, wildtype_lengths)
   print(f"T-test: t={t_stat:.3f}, p={p_value:.2e}")
   
   # Effect size (Cohen's d)
   pooled_std = np.sqrt(((len(mutated_lengths) - 1) * np.var(mutated_lengths) + 
                        (len(wildtype_lengths) - 1) * np.var(wildtype_lengths)) /
                       (len(mutated_lengths) + len(wildtype_lengths) - 2))
   cohens_d = (np.mean(mutated_lengths) - np.mean(wildtype_lengths)) / pooled_std
   print(f"Cohen's d: {cohens_d:.3f}")
   
   # Visualization
   fig, axes = plt.subplots(1, 3, figsize=(15, 5))
   
   # Histogram comparison
   axes[0].hist(wildtype_lengths, bins=50, alpha=0.7, label='Wildtype', density=True)
   axes[0].hist(mutated_lengths, bins=50, alpha=0.7, label='Mutated', density=True)
   axes[0].set_xlabel('Fragment Length (bp)')
   axes[0].set_ylabel('Density')
   axes[0].set_title('Length Distribution: Mutated vs Wildtype')
   axes[0].legend()
   
   # Box plot
   axes[1].boxplot([wildtype_lengths, mutated_lengths], labels=['Wildtype', 'Mutated'])
   axes[1].set_ylabel('Fragment Length (bp)')
   axes[1].set_title('Length Distribution Box Plot')
   
   # Q-Q plot
   from scipy.stats import probplot
   probplot(mutated_lengths, dist="norm", plot=axes[2])
   axes[2].set_title('Q-Q Plot: Mutated Fragment Lengths')
   
   plt.tight_layout()
   plt.savefig('mutation_length_analysis.png', dpi=300, bbox_inches='tight')
   plt.show()

Advanced GMM Analysis
---------------------

Explore different numbers of Gaussian components:

.. code-block:: python

   from pyfraglib import Fragment
   from pyfraglib.math import fit_gmm, goodness_of_fit_stats
   import numpy as np
   import matplotlib.pyplot as plt
   
   # Load fragments
   fragments = Fragment.from_bam("sample.bam", "variants.vcf")
   lengths = np.array([f.length for f in fragments if not f.is_bogus], dtype=np.float64)
   
   # Test different numbers of components
   n_components_range = range(1, 6)
   aic_scores = []
   bic_scores = []
   ks_stats = []
   
   for n in n_components_range:
       # Create configuration for n components
       config = {
           "number_of_gaussians": n,
           "subsample_percentage": 0.1,
           "means_lower_bounds": [50 + i*50 for i in range(n)],
           "means_upper_bounds": [150 + i*50 for i in range(n)],
           "std_lower_bounds": [10] * n,
           "std_upper_bounds": [50] * n,
           "initial_means": [100 + i*50 for i in range(n)],
           "initial_stds": [25] * n,
           "initial_weights": [1.0/n] * n
       }
       
       # Save configuration
       with open(f"gmm_config_{n}.json", "w") as f:
           json.dump(config, f, indent=2)
       
       # Fit GMM
       result, n_comp, params, data_subset = fit_gmm(lengths, f"gmm_config_{n}.json")
       
       # Calculate information criteria
       log_likelihood = -result.fun
       aic = 2 * (3 * n) - 2 * log_likelihood
       bic = np.log(len(lengths)) * (3 * n) - 2 * log_likelihood
       
       aic_scores.append(aic)
       bic_scores.append(bic)
       
       # Calculate goodness of fit
       gof_stats = goodness_of_fit_stats(lengths, params, n)
       ks_stats.append(gof_stats['kolmogorov_smirnov_statistic'])
   
   # Plot model selection criteria
   fig, axes = plt.subplots(1, 3, figsize=(15, 5))
   
   axes[0].plot(n_components_range, aic_scores, 'o-', label='AIC')
   axes[0].set_xlabel('Number of Components')
   axes[0].set_ylabel('AIC Score')
   axes[0].set_title('Akaike Information Criterion')
   axes[0].grid(True)
   
   axes[1].plot(n_components_range, bic_scores, 'o-', label='BIC', color='orange')
   axes[1].set_xlabel('Number of Components')
   axes[1].set_ylabel('BIC Score')
   axes[1].set_title('Bayesian Information Criterion')
   axes[1].grid(True)
   
   axes[2].plot(n_components_range, ks_stats, 'o-', label='KS Statistic', color='green')
   axes[2].set_xlabel('Number of Components')
   axes[2].set_ylabel('KS Statistic')
   axes[2].set_title('Kolmogorov-Smirnov Statistic')
   axes[2].grid(True)
   
   plt.tight_layout()
   plt.savefig('model_selection.png', dpi=300, bbox_inches='tight')
   plt.show()
   
   # Find optimal number of components
   optimal_aic = n_components_range[np.argmin(aic_scores)]
   optimal_bic = n_components_range[np.argmin(bic_scores)]
   optimal_ks = n_components_range[np.argmin(ks_stats)]
   
   print(f"Optimal number of components:")
   print(f"  AIC: {optimal_aic}")
   print(f"  BIC: {optimal_bic}")
   print(f"  KS: {optimal_ks}")

See Also
--------

* :doc:`basic_workflow` - Basic fragmentomics workflow
* :doc:`mutation_analysis` - Mutation analysis techniques
* :doc:`fragmentomics_scores` - Computing fragmentomics scores
* :func:`pyfraglib.math.fit_gmm` - Gaussian mixture model fitting
* :func:`pyfraglib.math.goodness_of_fit_stats` - Goodness of fit statistics