"""
Tissue Mixture cfDNA Simulation
===============================

Advanced cfDNA simulation module for modeling complex biological scenarios
involving multiple tissue sources, disease progression, and clinical
conditions. This module extends the base FragmentSimulator to handle realistic
multi-tissue mixtures commonly encountered in liquid biopsy applications.

Key Applications
----------------
- **Cancer Detection**: Tumor-derived cfDNA mixed with normal tissue background
- **Organ Transplant Monitoring**: Donor cfDNA in recipient plasma
- **Disease Progression**: Longitudinal cfDNA changes over time
- **Therapy Response**: cfDNA dynamics during treatment

Biological Foundation
---------------------
The simulator incorporates tissue-specific characteristics based on:

1. **Fragment Size Distributions**: Tissue-specific size patterns
   - Tumor cfDNA: Variable size distributions
   - Normal tissue: Standard nucleosomal patterns (~167bp)

2. **Nucleosome Positioning**: Tissue-specific chromatin organization
   - Different nucleosome spacing patterns
   - Tissue-specific chromatin accessibility
   - Methylation-dependent fragmentation

3. **End Motif Preferences**: Tissue-specific cleavage patterns
   - Methylation-dependent DNASE1L3 activity
   - Tissue-specific nuclease expression levels
   - Disease-associated fragmentation changes

4. **Clinical Signatures**: Disease-specific alterations
   - Cancer progression signatures
   - Fetal fraction dynamics
   - Therapy response patterns

Example Usage
-------------
.. code-block:: python

    from pyfraglib.simulator.tissue_mixture_simulator import (
        TissueMixtureSimulator, TissueProfile
    )

    # Initialize simulator
    simulator = TissueMixtureSimulator("/path/to/reference.fasta")

    # Simulate cancer cfDNA mixture
    cancer_fragments = simulator.simulate_tissue_mixture(
        tissue_types=["hematopoietic", "tumor"],
        tissue_fractions=[0.95, 0.05],  # 5% tumor fraction
        total_fragments=10000,
        genomic_regions=[("chr1", 1000000, 1100000)]
    )

    # Cancer progression study
    progression_data = simulator.simulate_cancer_progression(
        normal_profile="hematopoietic",
        tumor_fractions=[0.01, 0.05, 0.15, 0.30],
        time_points=["baseline", "3months", "6months", "12months"],
        fragments_per_timepoint=20000,
        genomic_regions=[("chr1", 1000000, 2000000)]
    )

Tissue Profiles
---------------
Predefined tissue profiles based on literature:

- **hematopoietic**: Blood cell-derived cfDNA (normal background)
- **liver**: Hepatocyte-derived cfDNA (organ damage monitoring)
- **tumor**: Generic tumor-derived cfDNA (cancer detection)

Custom tissue profiles can be defined using the TissueProfile dataclass.

References
----------
- Wan et al. (2017) Nature Reviews: cfDNA fragmentation patterns
- Cristiano et al. (2019) Nature: Fragmentomics for cancer detection
- Sun et al. (2019) PNAS: Tissue-specific cfDNA fragmentation

License
-------
Copyright (C) 2025 Daniel Schütte, daniel.schuette@iccb-cologne.org

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version. This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
more details. You should have received a copy of the GNU General Public
License along with this program. If not, see <https://www.gnu.org/licenses/>.
"""
import logging

import numpy as np

from dataclasses import dataclass
from typing import Final
from pyfraglib.core import fail
from pyfraglib.fragment import FragmentList
from pyfraglib.simulator.fragment_simulator import FragmentSimulator, \
                                                   NucleaseProfile

LOGGER_NAME: Final[str] = "pyfraglib.simulator.mixture"


@dataclass
class TissueProfile:
    """
    Comprehensive tissue-specific fragmentation parameter profile.

    This class encapsulates all tissue-specific parameters that influence
    cfDNA fragmentation patterns. Each tissue type has unique characteristics
    that affect fragment size distributions, nucleosome positioning, chromatin
    accessibility, and nuclease activity patterns.

    Notes
    -----
    Tissue profiles are based on experimental observations from:
    - MNase-seq data for nucleosome positioning
    - ATAC-seq/DNase-seq for chromatin accessibility
    - Bisulfite sequencing for methylation patterns
    - Clinical cfDNA studies for fragment size distributions

    The parameters interact to produce realistic fragmentation:
    - Higher chromatin openness increases base cleavage probability
    - Tighter nucleosome spacing affects fragment size periodicity
    - Methylation levels influence DNASE1L3 activity and CC motif preference
    - End motif preferences modulate cleavage site selection

    Examples
    --------
    >>> # Create custom tissue profile
    >>> custom_tissue = TissueProfile(
    ...     name="custom_cancer",
    ...     fragment_size_distribution={
    ...         "mean": 155, "std": 15, "short_fraction": 0.3
    ...     },
    ...     nucleosome_spacing=180,
    ...     chromatin_openness=0.8,
    ...     end_motif_preferences={"CC": 1.4, "AT": 1.1},
    ...     methylation_level=0.4
    ... )
    """
    #: Tissue identifier (e.g., "hematopoietic", "liver", "tumor").
    #: Used for logging and identification purposes.
    name: str

    #: Parameters defining fragment size distribution:
    #:
    #: - "mean": Mean fragment size (typically 140-170bp)
    #: - "std": Standard deviation of size distribution
    #: - "short_fraction": Proportion of short fragments (<150bp)
    fragment_size_distribution: dict[str, float]

    #: Average nucleosome spacing in base pairs. Typical values:
    #:
    #: - Normal tissues: 185-190bp
    #: - Tumor tissues: 175-185bp (more compact)
    #: - Fetal tissues: 175-180bp
    nucleosome_spacing: float

    #: Relative chromatin accessibility (0.0-1.0 scale):
    #:
    #: - 0.0: Completely closed chromatin
    #: - 0.5: Moderately accessible
    #: - 1.0: Highly accessible (e.g., hematopoietic)
    chromatin_openness: float

    #: Relative preferences for specific end motifs:
    #:
    #: - Keys: Dinucleotide motifs ("CC", "CG", "AT", etc.)
    #: - Values: Multiplicative factors (1.0 = neutral, >1.0 = preferred)
    end_motif_preferences: dict[str, float]

    #: Average DNA methylation level (0.0-1.0):
    #:
    #: - Affects DNASE1L3 activity (higher methylation = more activity)
    #: - Typical values: 0.3-0.8 depending on tissue and disease state
    methylation_level: float


TISSUE_PROFILES = {
    "hematopoietic": TissueProfile(
        name="hematopoietic",
        fragment_size_distribution={
            "mean": 167, "std": 10, "short_fraction": 0.1
        },
        nucleosome_spacing=185,
        chromatin_openness=0.7,
        end_motif_preferences={"CC": 1.2, "CG": 1.1},
        methylation_level=0.7
    ),
    "liver": TissueProfile(
        name="liver",
        fragment_size_distribution={
            "mean": 165, "std": 12, "short_fraction": 0.15
        },
        nucleosome_spacing=190,
        chromatin_openness=0.5,
        end_motif_preferences={"CC": 1.3, "CT": 1.1},
        methylation_level=0.75
    ),
    "tumor": TissueProfile(
        name="tumor",
        fragment_size_distribution={
            "mean": 150, "std": 18, "short_fraction": 0.25
        },
        nucleosome_spacing=175,
        chromatin_openness=0.8,
        end_motif_preferences={"AA": 1.2, "TT": 1.2},
        methylation_level=0.5  # variable methylation
    )
}


@dataclass
class DiseaseSignature:
    """
    Disease-specific fragmentation pattern alterations.

    This class defines systematic changes in cfDNA fragmentation patterns
    associated with specific diseases or pathological conditions. These
    signatures can be applied to modify baseline tissue profiles to model
    disease states.

    Notes
    -----
    Disease signatures are typically applied stochastically:

    - Not all fragments from diseased tissue show the signature
    - Application probability often correlates with disease severity
    - Multiple signatures may be combined for complex conditions

    Common disease patterns:

    - **Cancer**: Shorter fragments, reduced periodicity, altered motifs
    - **Inflammation**: Increased apoptotic fragmentation patterns
    - **Organ damage**: Tissue-specific fragmentation changes

    Examples
    --------
    >>> cancer_sig = DiseaseSignature(
    ...     name="aggressive_cancer",
    ...     size_shift=-20.0,
    ...     periodicity_change=0.7,
    ...     preferred_cut_sites=["CATG", "GATC"],
    ...     aberrant_ends={"CC": 1.5, "AA": 0.8}
    ... )
    """
    #: Disease identifier (e.g., "breast_cancer", "lung_cancer", "sepsis").
    name: str

    #: Systematic shift in mean fragment size (in base pairs):
    #:
    #: - Negative values: Shorter fragments (common in cancer)
    #: - Positive values: Longer fragments (rare)
    #: - Typical range: -30 to +10 bp
    size_shift: float

    #: Alteration in 10bp nucleosomal periodicity strength:
    #:
    #: - 1.0: No change in periodicity
    #: - <1.0: Reduced periodicity (disrupted nucleosome positioning)
    #: - >1.0: Enhanced periodicity
    periodicity_change: float

    #: Disease-specific preferred cleavage motifs.
    #: List of DNA motifs that show increased cleavage in disease state.
    preferred_cut_sites: list[str]

    #: Abnormal end motif frequencies in disease:
    #:
    #: - Keys: Motif sequences
    #: - Values: Relative frequency changes (1.0 = no change)
    aberrant_ends: dict[str, float]


class TissueMixtureSimulator(FragmentSimulator):
    """
    cfDNA simulator for multi-tissue mixtures and clinical scenarios.

    This class extends the base FragmentSimulator to model complex biological
    scenarios where cfDNA originates from multiple tissue sources with
    different characteristics.

    The simulator handles tissue mixing, disease progression modeling, and
    clinical condition simulation with biologically accurate parameters
    derived from literature and experimental observations.

    Key Features
    ------------
    Multi-tissue Mixing:

    - Proportional mixing of different tissue sources
    - Tissue-specific fragmentation characteristics
    - Noise modeling

    Cancer Progression:

    - Longitudinal tumor fraction changes
    - Cancer-specific fragmentation signatures

    Clinical Realism:

    - Biological and technical noise addition
    - Fragment loss modeling (extraction efficiency)
    - End repair artifacts simulation

    Examples
    --------
    Basic tissue mixture simulation:

    >>> simulator = TissueMixtureSimulator(\"/path/to/reference.fasta\")
    >>> fragments = simulator.simulate_tissue_mixture(
    ...     tissue_types=[\"hematopoietic\", \"tumor\"],
    ...     tissue_fractions=[0.90, 0.10],
    ...     total_fragments=10000,
    ...     genomic_regions=[(\"chr1\", 1000000, 1100000)]
    ... )

    Cancer progression study:

    >>> progression = simulator.simulate_cancer_progression(
    ...     normal_profile=\"hematopoietic\",
    ...     tumor_fractions=[0.01, 0.05, 0.15],
    ...     time_points=[\"diagnosis\", \"3months\", \"6months\"],
    ...     fragments_per_timepoint=20000,
    ...     genomic_regions=[(\"chr1\", 1000000, 2000000)],
    ...     size_shift=-20.0, short_enrichment=0.25
    ... )

    References
    ----------
    - Wan et al. (2017) Nature Reviews Genetics: cfDNA fragmentation patterns
    - Sun et al. (2019) PNAS: Tissue-specific fragmentation signatures
    - Cristiano et al. (2019) Nature: Cancer detection via fragmentomics
    """
    def __init__(self, fasta_path: str, **kwargs) -> None:  # type: ignore
        super().__init__(fasta_path, **kwargs)  # type: ignore
        self.tissue_profiles = TISSUE_PROFILES.copy()
        self.logger = logging.getLogger(LOGGER_NAME)

    def add_custom_tissue(self, tissue_profile: TissueProfile) -> None:
        """
        Add a custom tissue profile to the simulator's available profiles.

        This method allows users to define custom tissue types with specific
        fragmentation characteristics beyond the predefined profiles. Custom
        profiles can model rare tissues, disease-specific alterations, or
        experimental conditions.

        Parameters
        ----------
        tissue_profile : TissueProfile
            Complete tissue profile definition including name, fragment size
            distribution, nucleosome spacing, chromatin openness, end motif
            preferences, and methylation level.

        Examples
        --------
        >>> simulator = TissueMixtureSimulator("/path/to/reference.fasta")
        >>> kidney_profile = TissueProfile(
        ...     name="kidney",
        ...     fragment_size_distribution={
        ...         "mean": 162, "std": 11, "short_fraction": 0.12
        ...     },
        ...     nucleosome_spacing=188,
        ...     chromatin_openness=0.6,
        ...     end_motif_preferences={"CG": 1.1, "GC": 1.1},
        ...     methylation_level=0.72
        ... )
        >>> simulator.add_custom_tissue(kidney_profile)
        >>> # Now "kidney" can be used in tissue_types lists
        """
        self.tissue_profiles[tissue_profile.name] = tissue_profile
        self.logger.info(f"Added custom tissue profile: {tissue_profile.name}")

    def simulate_tissue_mixture(
        self, tissue_types: list[str], tissue_fractions: list[float],
        total_fragments: int, genomic_regions: list[tuple[str, int, int]],
        add_noise: bool = True
    ) -> FragmentList:
        """
        Generate realistic cfDNA mixture from multiple tissue sources.

        This method creates a biologically realistic mixture of cfDNA fragments
        originating from different tissue types, each contributing a specified
        fraction of the total fragments. The method applies tissue-specific
        fragmentation characteristics and can optionally add biological noise.

        Parameters
        ----------
        tissue_types : list[str]
            List of tissue type names to include in the mixture. Must be
            present in tissue_profiles (either predefined or custom).
        tissue_fractions : list[float]
            Fraction contribution from each tissue type. Must sum to 1.0
            (within 0.001 tolerance). Order corresponds to tissue_types.
        total_fragments : int
            Total number of fragments to generate across all tissues.
            Individual tissue contributions are calculated proportionally.
        genomic_regions : list[tuple[str, int, int]]
            List of genomic regions (chromosome, start, end) to sample
            fragments from. Fragments are distributed across regions.
        add_noise : bool, default=True
            Whether to add realistic biological and technical noise:
            - Random fragment loss (~5%)
            - Size variation (±2bp standard deviation)
            - End repair artifacts (~2%)

        Returns
        -------
        FragmentList
            Combined fragment list containing fragments from all tissues
            with tissue-specific characteristics. Fragments are shuffled
            to simulate realistic mixing.

        Raises
        ------
        SystemExit
            If tissue fractions don't sum to 1.0 or if unknown tissue
            types are requested.

        Notes
        -----
        Each tissue contributes fragments according to its profile:
        - Fragment size distribution specific to tissue type
        - Nuclease activity profile based on tissue characteristics
        - End motif preferences reflecting tissue methylation state
        - Chromatin accessibility affecting cleavage patterns

        The simulation process:
        1. Calculate fragment count per tissue (proportional allocation)
        2. Generate tissue-specific nuclease profiles
        3. Simulate fragments for each tissue using base FragmentSimulator
        4. Combine fragments from all tissues
        5. Apply biological noise if requested

        Examples
        --------
        Cancer detection mixture (5% tumor fraction):

        >>> fragments = simulator.simulate_tissue_mixture(
        ...     tissue_types=["hematopoietic", "tumor"],
        ...     tissue_fractions=[0.95, 0.05],
        ...     total_fragments=50000,
        ...     genomic_regions=[("chr1", 1000000, 1500000)],
        ...     add_noise=True
        ... )
        """
        if abs(sum(tissue_fractions) - 1.0) > 0.001:
            fail("Tissue fractions must sum to 1.0")

        if len(tissue_types) != len(tissue_fractions):
            fail("Tissue number must match number of fractions")

        sim_info: dict[str, object] = dict(zip(tissue_types, tissue_fractions))
        self.logger.info(f"Simulating mixture: {sim_info}")

        all_fragments = FragmentList()
        for tissue, fraction in zip(tissue_types, tissue_fractions):
            if tissue not in self.tissue_profiles:
                fail(f"Unknown tissue type: {tissue}")

            profile = self.tissue_profiles[tissue]
            num_frags: int = int(total_fragments * fraction)
            if num_frags == 0:
                continue

            frags_per_region = num_frags // len(genomic_regions)
            for chrom, start, end in genomic_regions:
                nuclease_profile = self._get_tissue_nuclease_profile(profile)
                tissue_fragments = self._simulate_tissue_fragments(
                    chrom, start, end,
                    frags_per_region,
                    profile,
                    nuclease_profile
                )
                for fragment in tissue_fragments:
                    all_fragments.append(fragment)

        if add_noise:
            all_fragments = self._add_biological_noise(all_fragments)

        return all_fragments

    def _get_tissue_nuclease_profile(
        self, tissue: TissueProfile
    ) -> NucleaseProfile:
        """
        Generate tissue-specific nuclease activity profile.

        This method converts tissue characteristics into nuclease activity
        parameters that affect fragmentation patterns. The conversion is
        based on known relationships between tissue properties and nuclease
        expression/activity levels.

        Parameters
        ----------
        tissue : TissueProfile
            Tissue profile containing chromatin and methylation
            characteristics.

        Returns
        -------
        NucleaseProfile
            Nuclease activity profile customized for the tissue type.

        Notes
        -----
        Tissue-to-nuclease parameter mapping:

        **DNase I Activity**:
        - Scaled by chromatin openness
        - Higher in accessible tissues (e.g., hematopoietic)
        - Lower in compact tissues (e.g., liver)

        **DNase1L3 Activity**:
        - Influenced by methylation levels
        - Higher methylation enhances activity
        - Primary cfDNA fragmentation nuclease

        **DFFB Activity**:
        - Constant baseline level (1.0)
        - May be elevated in apoptotic conditions
        - Affects linker region cleavage

        **Motif Preferences**:
        - Incorporates tissue-specific end motif preferences
        - Methylation-dependent CC preferences
        - Tissue-specific sequence biases
        """
        dnase1_activity = 1.0 * tissue.chromatin_openness
        dnase1l3_activity = 1.0 + 0.5 * (tissue.methylation_level - 0.5)
        dffb_activity = 1.0
        dnase1l3_prefs = {"CC": 1.5 * tissue.methylation_level}
        dnase1l3_prefs.update(tissue.end_motif_preferences)

        return NucleaseProfile(
            dnase1_activity=dnase1_activity,
            dnase1l3_activity=dnase1l3_activity,
            dffb_activity=dffb_activity,
            dnase1l3_motif_preference=dnase1l3_prefs
        )

    def _simulate_tissue_fragments(
        self, chrom: str, start: int, end: int, num_fragments: int,
        tissue: TissueProfile, nuclease_profile: NucleaseProfile
    ) -> FragmentList:
        """
        Generate fragments with tissue-specific characteristics.

        This method creates fragments using the base FragmentSimulator
        but with parameters customized for a specific tissue type. It
        serves as a bridge between tissue profiles and the core simulation
        engine.

        Parameters
        ----------
        chrom : str
            Chromosome name for fragment generation.
        start : int
            Start position of genomic region.
        end : int
            End position of genomic region.
        num_fragments : int
            Number of fragments to generate for this tissue.
        tissue : TissueProfile
            Complete tissue profile with all characteristics.
        nuclease_profile : NucleaseProfile
            Tissue-specific nuclease activity profile.

        Returns
        -------
        FragmentList
            Fragment collection with tissue-specific properties.

        Notes
        -----
        The method translates tissue profile parameters into simulation
        parameters:
        - Fragment size distribution -> fragment_size_params
        - Tissue name -> tissue_type for simulation context
        - Nuclease profile -> cleavage preferences
        - All other tissue characteristics handled by base simulator
        """
        fragment_size_params = {
            "mean": tissue.fragment_size_distribution["mean"],
            "std": tissue.fragment_size_distribution["std"],
            "short_fraction": tissue.fragment_size_distribution[
                "short_fraction"
            ]
        }
        fragments = self.simulate_fragments(
            chrom, start, end,
            num_fragments,
            tissue_type=tissue.name,
            nuclease_profile=nuclease_profile,
            fragment_size_params=fragment_size_params
        )

        return fragments

    def _add_biological_noise(self, fragments: FragmentList) -> FragmentList:
        """
        Apply biological and technical noise to fragment collection.

        This method simulates various sources of noise and artifacts that
        occur during cfDNA extraction, library preparation, and sequencing.
        These effects are important for creating realistic training data
        and testing the robustness of analytical methods.

        Parameters
        ----------
        fragments : FragmentList
            Input fragment list to which noise will be applied.

        Returns
        -------
        FragmentList
            Modified fragment list with applied noise effects.

        Notes
        -----
        Applied noise sources:

        **Fragment Loss (5% random loss)**:
        - Simulates inefficient DNA extraction
        - PCR amplification bias
        - Size-dependent purification losses

        **Size Variation (±2bp standard deviation)**:
        - Reflects measurement uncertainty
        - Size calling accuracy limitations
        - Biological size variability

        **End Repair Artifacts (2% of fragments)**:
        - Terminal nucleotide damage/modification
        - Library preparation artifacts
        - Sequencing error effects
        - Represented as 'N' in terminal position

        These noise levels are based on typical experimental observations
        and can be adjusted for specific applications or platforms.
        """
        noisy_fragments = FragmentList()

        for fragment in fragments:
            if np.random.random() < 0.05:
                continue

            size_noise = np.random.normal(0, 2)
            fragment.length = int(max(50, fragment.length + size_noise))
            fragment.end_pos = fragment.start_pos + fragment.length

            if np.random.random() < 0.02:
                fragment.end5p = fragment.end5p[:-1] + "N"

            noisy_fragments.append(fragment)

        return noisy_fragments

    def simulate_cancer_progression(
        self, normal_profile: str, tumor_fractions: list[float],
        time_points: list[str], fragments_per_timepoint: int,
        genomic_regions: list[tuple[str, int, int]],
        size_shift: float = -15.0, short_enrichment: float = 0.2
    ) -> dict[str, FragmentList]:
        """
        Model longitudinal cancer progression through cfDNA changes.

        This method simulates the evolution of cfDNA characteristics over
        the course of cancer progression, modeling increasing tumor fractions
        and the development of cancer-specific fragmentation signatures.
        It is useful for studying early detection, monitoring disease
        progression, and therapy response.

        Parameters
        ----------
        normal_profile : str
            Name of normal tissue profile to use as background (typically
            "hematopoietic" for blood-based cfDNA).
        tumor_fractions : list[float]
            Tumor-derived cfDNA fraction at each time point. Values should
            be in [0, 1] range and typically increase over time for
            progression studies.
        time_points : list[str]
            Labels for each time point (e.g., ["baseline", "3months"]).
            Must have same length as tumor_fractions.
        fragments_per_timepoint : int
            Number of fragments to generate at each time point. Consistent
            across time points for comparative analysis.
        genomic_regions : list[tuple[str, int, int]]
            Genomic regions (chromosome, start, end) to sample fragments
            from. Same regions used for all time points.
        size_shift : float, default=-15.0
            Systematic shift in fragment size (in base pairs).
            Negative values create shorter fragments (typical for cancer).
            Typical range: -30 to +10 bp.
        short_enrichment : float, default=0.2
            Enrichment factor for short fragments (<150bp).
            Higher values increase selection for very short fragments.
            Typical range: 0.1-0.3.

        Returns
        -------
        dict[str, FragmentList]
            Dictionary mapping time point labels to FragmentList objects.
            Each FragmentList contains the cfDNA profile for that time point
            with appropriate tumor fraction and cancer signatures.

        Raises
        ------
        SystemExit
            If tumor_fractions and time_points have different lengths.

        Notes
        -----
        Cancer progression modeling incorporates:

        **Tumor Fraction Changes**:
        - Proportional mixing of normal and tumor tissues
        - Realistic tumor fraction ranges (0.01-0.50 typical)
        - Time-dependent tumor evolution

        **Cancer Signatures**:
        - Fragment size shifts (typically shorter fragments)
        - Altered end motif preferences
        - Modified nucleosomal periodicity
        - Cancer-type specific patterns

        **Biological Realism**:
        - Progressive signature strength with tumor fraction
        - Stochastic application of cancer signatures
        - Background noise maintenance

        Examples
        --------
        Early stage cancer progression:

        >>> progression = simulator.simulate_cancer_progression(
        ...     normal_profile="hematopoietic",
        ...     tumor_fractions=[0.01, 0.03, 0.08, 0.15],
        ...     time_points=["diagnosis", "1month", "3months", "6months"],
        ...     fragments_per_timepoint=25000,
        ...     genomic_regions=[("chr1", 1000000, 2000000)],
        ...     size_shift=-20.0, short_enrichment=0.25
        ... )

        Therapy response monitoring:

        >>> response = simulator.simulate_cancer_progression(
        ...     normal_profile="hematopoietic",
        ...     tumor_fractions=[0.20, 0.15, 0.08, 0.03],  # Decreasing
        ...     time_points=["pre_treatment", "1cycle", "3cycles", "6cycles"],
        ...     fragments_per_timepoint=30000,
        ...     genomic_regions=[("chr1", 1000000, 2000000)],
        ...     size_shift=-15.0, short_enrichment=0.2
        ... )
        """
        if len(tumor_fractions) != len(time_points):
            fail("Tumor fractions must match time points")

        results = {}
        for tp, tf in zip(time_points, tumor_fractions):
            self.logger.info(f"Simulating {tp} with tumor fraction {tf}")

            tissue_types = [normal_profile, "tumor"]
            tissue_fractions = [1 - tf, tf]
            fragments = self.simulate_tissue_mixture(
                tissue_types,
                tissue_fractions,
                fragments_per_timepoint,
                genomic_regions
            )

            if tf > 0:
                fragments = self._add_cancer_signatures(
                    fragments, tf, size_shift, short_enrichment
                )

            results[tp] = fragments

        return results

    def _add_cancer_signatures(
        self, fragments: FragmentList, tumor_fraction: float,
        size_shift: float = -15.0, short_enrichment: float = 0.2
    ) -> FragmentList:
        """
        Apply cancer-specific fragmentation signatures to fragment collection.

        This method modifies fragments to incorporate known cancer-associated
        fragmentation patterns. The signatures are applied stochastically
        based on the tumor fraction, reflecting the mixed nature of cfDNA
        in cancer patients.

        Parameters
        ----------
        fragments : FragmentList
            Input fragment collection to modify.
        tumor_fraction : float
            Proportion of tumor-derived fragments (0.0-1.0).
            Determines probability of signature application.
        size_shift : float, default=-15.0
            Systematic shift in fragment size (in base pairs).
            Negative values create shorter fragments (typical for cancer).
            Typical range: -30 to +10 bp.
        short_enrichment : float, default=0.2
            Enrichment factor for short fragments (<150bp).
            Higher values increase selection for very short fragments.
            Typical range: 0.1-0.3.

        Returns
        -------
        FragmentList
            Modified fragment list with cancer signatures applied.

        Notes
        -----
        Cancer-specific signatures applied:

        **Fragment Size Changes**:
        - Systematic size shifts using size_shift parameter
        - Applied stochastically based on tumor fraction
        - Creates shorter fragments typical of cancer cfDNA

        **Short Fragment Enrichment**:
        - Increased proportion of <150bp fragments
        - Enhanced selection controlled by short_enrichment parameter
        - Reflects cancer-associated fragmentation patterns

        **Application Strategy**:
        - Signatures applied probabilistically to fragments
        - Application probability = tumor_fraction
        - Some fragments may receive multiple modifications
        - Maintains realistic heterogeneity

        The signatures are based on clinical observations from:
        - Liquid biopsy studies
        - Fragmentomics research
        - Cancer-specific cfDNA patterns
        """

        modified_fragments = FragmentList()
        for fragment in fragments:
            # Apply cancer signatures to simulate tumor-derived cfDNA characteristics
            if np.random.random() < tumor_fraction:
                new_length = int(fragment.length + size_shift
                                 * np.random.random())
                new_length = max(50, new_length)
                fragment.length = new_length
                fragment.end_pos = fragment.start_pos + new_length
                if (fragment.length < 150 and
                        np.random.random() < short_enrichment):
                    modified_fragments.append(fragment)

            modified_fragments.append(fragment)

        return modified_fragments
