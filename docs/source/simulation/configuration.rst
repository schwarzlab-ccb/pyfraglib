Configuration Files
===================

Configuration files provide a flexible way to specify simulation parameters and enable reproducible simulation workflows. The pyfraglib simulation module supports JSON configuration files with comprehensive parameter specification.

Basic Configuration Structure
-----------------------------

All configuration files follow this basic structure:

.. code-block:: json

   {
       "fasta_path": "/path/to/reference.fasta",
       "output_name": "simulation_name",
       "regions": [
           {"chr": "chr1", "start": 1000000, "end": 2000000},
           {"chr": "chr2", "start": 3000000, "end": 4000000}
       ],
       "n_fragments": 50000,
       "simulation_mode": "basic"
   }

Required Parameters
-------------------

All configuration files must specify:

* **fasta_path** (string): Path to indexed reference FASTA file
* **output_name** (string): Name for output files
* **regions** (array): Genomic regions to simulate from
* **n_fragments** (integer): Number of fragments to generate
* **simulation_mode** (string): Type of simulation ("basic", "tissue_mixture", "cancer_progression", "fetal_fraction")

Basic Simulation Configuration
------------------------------

Simple single-tissue simulation:

.. code-block:: json

   {
       "fasta_path": "/data/reference/hg38.fasta",
       "output_name": "synthetic_basic",
       "regions": [
           {"chr": "chr1", "start": 1000000, "end": 2000000},
           {"chr": "chr2", "start": 3000000, "end": 4000000},
           {"chr": "chr3", "start": 5000000, "end": 6000000}
       ],
       "n_fragments": 25000,
       "simulation_mode": "basic",
       "length_distribution": {
           "components": [
               {"mean": 145, "std": 20, "weight": 0.3},
               {"mean": 167, "std": 15, "weight": 0.5},
               {"mean": 190, "std": 30, "weight": 0.2}
           ],
           "bounds": {
               "min_length": 50,
               "max_length": 500
           }
       }
   }

Tissue Mixture Configuration
----------------------------

Multi-tissue simulation with specified fractions:

.. code-block:: json

   {
       "fasta_path": "/data/reference/hg38.fasta",
       "output_name": "tissue_mixture_simulation",
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
           "hematopoietic": {
               "mean_length": 167,
               "std_length": 25,
               "cleavage_preferences": {
                   "CCCA": 1.2,
                   "TTTC": 1.1,
                   "AAAG": 0.9
               }
           },
           "liver": {
               "mean_length": 175,
               "std_length": 30,
               "cleavage_preferences": {
                   "CCCA": 1.3,
                   "GGGA": 1.2,
                   "TTTT": 0.8
               }
           },
           "tumor": {
               "mean_length": 155,
               "std_length": 35,
               "cleavage_preferences": {
                   "CCCA": 0.8,
                   "TTTC": 1.4,
                   "GGGA": 1.1
               }
           }
       }
   }

Cancer Progression Configuration
---------------------------------

Longitudinal cancer simulation with increasing tumor fractions:

.. code-block:: json

   {
       "fasta_path": "/data/reference/hg38.fasta",
       "output_name": "cancer_progression",
       "regions": [
           {"chr": "chr1", "start": 1000000, "end": 2000000},
           {"chr": "chr2", "start": 3000000, "end": 4000000}
       ],
       "n_fragments": 50000,
       "simulation_mode": "cancer_progression",
       "timepoints": [
           {
               "name": "baseline",
               "tumor_fraction": 0.01,
               "time_months": 0
           },
           {
               "name": "month_3",
               "tumor_fraction": 0.05,
               "time_months": 3
           },
           {
               "name": "month_6",
               "tumor_fraction": 0.15,
               "time_months": 6
           },
           {
               "name": "month_12",
               "tumor_fraction": 0.30,
               "time_months": 12
           }
       ],
       "background_tissues": {
           "hematopoietic": 0.85,
           "liver": 0.15
       }
   }

Fetal Fraction Configuration
----------------------------

NIPT simulation with varying fetal contributions:

.. code-block:: json

   {
       "fasta_path": "/data/reference/hg38.fasta",
       "output_name": "fetal_fraction_simulation",
       "regions": [
           {"chr": "chr13", "start": 1000000, "end": 2000000},
           {"chr": "chr18", "start": 1000000, "end": 2000000},
           {"chr": "chr21", "start": 1000000, "end": 2000000}
       ],
       "n_fragments": 200000,
       "simulation_mode": "fetal_fraction",
       "fetal_fractions": [0.05, 0.10, 0.15, 0.20, 0.25],
       "maternal_tissues": {
           "hematopoietic": 0.9,
           "placenta": 0.1
       },
       "gestational_age": {
           "weeks": 20,
           "length_modulation": {
               "fetal_shortening": 0.95,
               "maternal_stability": 1.0
           }
       }
   }

Advanced Configuration Options
------------------------------

Comprehensive configuration with all available parameters:

.. code-block:: json

   {
       "fasta_path": "/data/reference/hg38.fasta",
       "output_name": "comprehensive_simulation",
       "regions": [
           {"chr": "chr1", "start": 1000000, "end": 10000000},
           {"chr": "chr2", "start": 1000000, "end": 10000000},
           {"chr": "chr3", "start": 1000000, "end": 10000000}
       ],
       "n_fragments": 500000,
       "simulation_mode": "tissue_mixture",
       "tissue_fractions": {
           "hematopoietic": 0.70,
           "liver": 0.15,
           "tumor": 0.10,
           "placenta": 0.05
       },
       "length_distribution": {
           "global_params": {
               "min_length": 50,
               "max_length": 600,
               "outlier_fraction": 0.01
           },
           "tissue_specific": {
               "hematopoietic": {
                   "components": [
                       {"mean": 145, "std": 20, "weight": 0.3},
                       {"mean": 167, "std": 15, "weight": 0.5},
                       {"mean": 190, "std": 30, "weight": 0.2}
                   ]
               },
               "liver": {
                   "components": [
                       {"mean": 150, "std": 25, "weight": 0.2},
                       {"mean": 175, "std": 20, "weight": 0.6},
                       {"mean": 200, "std": 35, "weight": 0.2}
                   ]
               },
               "tumor": {
                   "components": [
                       {"mean": 140, "std": 30, "weight": 0.4},
                       {"mean": 165, "std": 25, "weight": 0.4},
                       {"mean": 200, "std": 40, "weight": 0.2}
                   ]
               }
           }
       },
       "nuclease_profiles": {
           "global_preferences": {
               "AT_rich_regions": 1.1,
               "GC_rich_regions": 0.9,
               "repetitive_elements": 0.8
           },
           "tissue_specific": {
               "hematopoietic": {
                   "cleavage_bias": {"A": 0.25, "T": 0.25, "G": 0.25, "C": 0.25},
                   "motif_preferences": {
                       "CCCA": 1.2, "TTTC": 1.1, "AAAG": 0.9, "GGGT": 1.0
                   }
               },
               "liver": {
                   "cleavage_bias": {"A": 0.3, "T": 0.3, "G": 0.2, "C": 0.2},
                   "motif_preferences": {
                       "CCCA": 1.3, "GGGA": 1.2, "TTTT": 0.8, "AAAC": 1.1
                   }
               },
               "tumor": {
                   "cleavage_bias": {"A": 0.2, "T": 0.4, "G": 0.2, "C": 0.2},
                   "motif_preferences": {
                       "CCCA": 0.8, "TTTC": 1.4, "GGGA": 1.1, "AAAT": 1.2
                   }
               }
           }
       },
       "quality_control": {
           "min_fragment_quality": 20,
           "max_n_content": 0.1,
           "bogus_fragment_rate": 0.05
       },
       "output_options": {
           "save_intermediate": true,
           "compression_level": 6,
           "include_metadata": true
       }
   }

Region Specification
--------------------

Genomic regions can be specified in multiple formats:

BED-style regions:

.. code-block:: json

   {
       "regions": [
           {"chr": "chr1", "start": 1000000, "end": 2000000},
           {"chr": "chr2", "start": 3000000, "end": 4000000}
       ]
   }

Chromosome-wide regions:

.. code-block:: json

   {
       "regions": [
           {"chr": "chr1", "start": 0, "end": -1},
           {"chr": "chr2", "start": 0, "end": -1}
       ]
   }

Gene-based regions:

.. code-block:: json

   {
       "regions": [
           {"gene": "BRCA1", "padding": 10000},
           {"gene": "TP53", "padding": 5000}
       ]
   }

Length Distribution Parameters
------------------------------

Detailed length distribution specification:

.. code-block:: json

   {
       "length_distribution": {
           "components": [
               {
                   "mean": 145,
                   "std": 20,
                   "weight": 0.3,
                   "bounds": {"min": 100, "max": 200}
               },
               {
                   "mean": 167,
                   "std": 15,
                   "weight": 0.5,
                   "bounds": {"min": 130, "max": 220}
               },
               {
                   "mean": 190,
                   "std": 30,
                   "weight": 0.2,
                   "bounds": {"min": 150, "max": 300}
               }
           ],
           "global_bounds": {
               "min_length": 50,
               "max_length": 500
           },
           "outlier_handling": {
               "clip_outliers": true,
               "outlier_threshold": 3.0
           }
       }
   }

Nuclease Profile Configuration
------------------------------

Specify nuclease cleavage preferences:

.. code-block:: json

   {
       "nuclease_profiles": {
           "default": {
               "name": "standard_dnase",
               "cleavage_bias": {
                   "A": 0.25,
                   "T": 0.25,
                   "G": 0.25,
                   "C": 0.25
               },
               "motif_preferences": {
                   "CCCA": 1.2,
                   "TTTC": 1.1,
                   "AAAG": 0.9,
                   "GGGT": 1.0
               },
               "sequence_context": {
                   "purine_rich": 1.1,
                   "pyrimidine_rich": 0.9,
                   "mixed": 1.0
               }
           },
           "tissue_specific": {
               "tumor": {
                   "name": "tumor_dnase",
                   "cleavage_bias": {
                       "A": 0.2,
                       "T": 0.4,
                       "G": 0.2,
                       "C": 0.2
                   },
                   "motif_preferences": {
                       "CCCA": 0.8,
                       "TTTC": 1.4,
                       "GGGA": 1.1,
                       "AAAT": 1.2
                   }
               }
           }
       }
   }

Quality Control Parameters
--------------------------

Configure quality control settings:

.. code-block:: json

   {
       "quality_control": {
           "fragment_quality": {
               "min_mapq": 20,
               "max_n_content": 0.1,
               "min_complexity": 0.3
           },
           "bogus_fragments": {
               "rate": 0.05,
               "types": ["n_content", "low_complexity", "invalid_coords"]
           },
           "validation": {
               "length_distribution_check": true,
               "motif_diversity_check": true,
               "tissue_fraction_tolerance": 0.05
           }
       }
   }

Output Configuration
--------------------

Control output format and options:

.. code-block:: json

   {
       "output_options": {
           "file_format": "frag",
           "compression": {
               "enabled": true,
               "level": 6,
               "algorithm": "gzip"
           },
           "metadata": {
               "include_config": true,
               "include_timestamp": true,
               "include_version": true
           },
           "intermediate_files": {
               "save_by_tissue": true,
               "save_by_region": false,
               "cleanup_after": true
           }
       }
   }

Configuration Validation
------------------------

Validate configuration files before simulation:

.. code-block:: python

   from pyfraglib.simulation import validate_config
   import json
   
   # Load configuration
   with open("simulation_config.json", "r") as f:
       config = json.load(f)
   
   # Validate configuration
   try:
       validation_result = validate_config(config)
       if validation_result["valid"]:
           print("Configuration is valid")
       else:
           print("Configuration errors:")
           for error in validation_result["errors"]:
               print(f"  - {error}")
   except Exception as e:
       print(f"Configuration validation failed: {e}")

Configuration Templates
-----------------------

Pre-defined templates for common scenarios:

.. code-block:: python

   from pyfraglib.simulation import ConfigTemplate
   
   # Create template for liquid biopsy
   liquid_biopsy_template = ConfigTemplate.liquid_biopsy(
       fasta_path="reference.fasta",
       output_name="liquid_biopsy",
       tumor_fraction=0.1,
       n_fragments=100000
   )
   
   # Create template for NIPT
   nipt_template = ConfigTemplate.nipt(
       fasta_path="reference.fasta",
       output_name="nipt_simulation",
       fetal_fraction=0.15,
       n_fragments=200000
   )
   
   # Create template for method validation
   validation_template = ConfigTemplate.method_validation(
       fasta_path="reference.fasta",
       output_name="validation_set",
       sample_types=["healthy", "cancer", "treated"],
       n_fragments=50000
   )
   
   # Save templates
   with open("liquid_biopsy_config.json", "w") as f:
       json.dump(liquid_biopsy_template, f, indent=2)

Command Line Usage
------------------

Use configuration files with the command line interface:

.. code-block:: bash

   # Basic simulation
   pyfrag.py simulate --config basic_config.json --out-dir output/
   
   # Tissue mixture simulation
   pyfrag.py simulate --config tissue_mixture_config.json --out-dir output/
   
   # Cancer progression simulation
   pyfrag.py simulate --config cancer_progression_config.json --out-dir output/
   
   # Validate configuration before running
   pyfrag.py simulate --config config.json --validate-only
   
   # Override specific parameters
   pyfrag.py simulate --config config.json --out-dir output/ --n-fragments 100000

Best Practices
--------------

1. **Configuration Management**:
   - Use version control for configuration files
   - Include metadata about simulation purpose
   - Use descriptive output names

2. **Parameter Selection**:
   - Start with default tissue profiles
   - Adjust parameters based on validation
   - Use realistic fragment counts

3. **Validation**:
   - Always validate configurations before large simulations
   - Compare simulated data with real data
   - Check tissue fraction accuracy

4. **Reproducibility**:
   - Set random seeds for reproducible results
   - Document configuration rationale
   - Save configuration with results

See Also
--------

* :doc:`overview` - Simulation overview
* :doc:`fragment_simulator` - Basic fragment simulation
* :doc:`tissue_mixture` - Tissue mixture simulation
* :doc:`../examples/simulation_examples` - Practical examples
* :class:`pyfraglib.simulation.ConfigTemplate` - Configuration templates