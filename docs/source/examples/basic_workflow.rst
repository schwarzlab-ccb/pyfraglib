Basic Workflow
==============

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

This workflow demonstrates the essential steps for fragmentomics analysis:

1. **Fragment Extraction**: Process BAM files with optional VCF for mutation annotation
2. **Data Persistence**: Save fragments to compressed .frag files for efficient storage
3. **Data Loading**: Load fragments from .frag files for analysis
4. **Quality Assessment**: Check fragment quality and mutation rates
5. **Statistical Analysis**: Fit Gaussian mixture models to fragment lengths

Best Practices
--------------

* Always check for bogus fragments before analysis
* Use .frag files for efficient storage and faster loading
* Close FragFile objects to free memory
* Consider using multiprocessing for large datasets with :func:`pyfraglib.Fragment.bams_to_frags`

See Also
--------

* :doc:`mutation_analysis` - Detailed mutation analysis workflows
* :doc:`fragmentomics_scores` - Computing WPS and motif diversity
* :doc:`batch_processing` - Processing multiple samples efficiently