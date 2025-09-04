"""
Tissue Mixture cfDNA Simulation
===============================

Simulation module for modeling more complex biological scenarios involving
multiple tissue sources and cancer progression. This module extends the base
``FragmentSimulator`` to handle multi-tissue mixtures.

Key Applications
----------------
- **Cancer Detection**: Tumor-derived cfDNA mixed with normal tissue background
- **Disease Progression**: Longitudinal cfDNA changes over time

Example Usage
-------------
.. code-block:: python

    from pyfraglib.simulator.tissue_mixture_simulator import (
        TissueMixtureSimulator, TissueProfile
    )
    from pyfraglib.simulator.fragment_simulator import SequenceContextGenerator

    # Initialize sequence generator and simulator
    seq_gen = SequenceContextGenerator("/path/to/reference.fasta")
    simulator = TissueMixtureSimulator(seq_gen)

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

    # Clean up
    seq_gen.close()

Tissue Profiles
---------------
Predefined tissue profiles:

- **hematopoietic**: Blood cell-derived cfDNA (normal background)
- **liver**: Hepatocyte-derived cfDNA (organ damage monitoring)
- **tumor**: Generic tumor-derived cfDNA (cancer detection)

Custom tissue profiles can be defined using the ``TissueProfile`` dataclass.

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
                                                   NucleaseProfile, \
                                                   SequenceContextGenerator

LOGGER_NAME: Final[str] = "pyfraglib.simulator.mixture"


@dataclass
class TissueProfile:
    """
    Tissue-specific fragmentation parameter profile.

    This class encapsulates tissue-specific parameters that influence cfDNA
    fragmentation patterns like fragment size distributions, chromatin
    accessibility, and nuclease activity patterns.

    Note
    ----
    Currently, chromatin accessibility can only be set through the alpha and
    beta parameters of beta distribution that generates the default nucleosome
    map. In future releases, data-based nucleosome maps will be properly
    handled

    Examples
    --------
    >>> custom_tissue = TissueProfile(
    ...     name="custom_cancer",
    ...     fragment_size_distribution={
    ...         "mean": 155, "std": 15, "mono_fraction": 0.8
    ...     },
    ...     nuclease_activities={
    ...         "dnase1_activity": 1.2,
    ...         "dnase1l3_activity": 1.0,
    ...         "dffb_activity": 1.5
    ...     },
    ...     chromatin_beta_params=(2.5, 4.5)  # Open, disrupted chromatin
    ... )
    """
    #: Tissue identifier (e.g., "hematopoietic", "liver", "tumor").
    name: str

    #: Parameters defining fragment size distribution. See
    #: `_generate_fragment_sizes` for possible parameters.
    fragment_size_distribution: dict[str, float]

    #: Tissue-specific nuclease activity levels. Available nucleases:
    #:
    #: - dnase1_activity: DNase I activity level (typical range: 0.5-1.5)
    #: - dnase1l3_activity: DNase1L3 activity level (typical range: 0.8-1.5)
    #: - dffb_activity: DFFB activity level (typical range: 0.1-2.0)
    nuclease_activities: dict[str, float]

    #: Beta distribution parameters for nucleosome occupancy:
    #:
    #: - alpha: Shape parameter α (higher = more occupied nucleosomes)
    #: - beta: Shape parameter β (higher = more open chromatin)
    #: - Typical combinations: open (2,5), normal (4,3), compact (6,2)
    chromatin_beta_params: tuple[float, float]


TISSUE_PROFILES = {
    "hematopoietic": TissueProfile(
        name="hematopoietic",
        fragment_size_distribution={
            "mean": 167, "std": 10, "mono_fraction": 0.85,
            "di_mean": 2*167, "di_std": 30
        },
        nuclease_activities={
            "dnase1_activity": 1.2,  # Higher due to open chromatin
            "dnase1l3_activity": 1.3,  # High baseline cfDNA nuclease
            "dffb_activity": 0.2      # Low apoptotic activity
        },
        chromatin_beta_params=(3.5, 3.5)  # Balanced, accessible chromatin
    ),
    "liver": TissueProfile(
        name="liver",
        fragment_size_distribution={
            "mean": 165, "std": 12, "mono_fraction": 0.80,
            "di_mean": 2*165, "di_std": 35
        },
        nuclease_activities={
            "dnase1_activity": 0.8,   # Lower due to compact chromatin
            "dnase1l3_activity": 1.4,  # High methylation-enhanced activity
            "dffb_activity": 0.5      # Moderate tissue damage
        },
        chromatin_beta_params=(5.0, 2.5)  # Compact, well-organized chromatin
    ),
    "tumor": TissueProfile(
        name="tumor",
        fragment_size_distribution={
            "mean": 150, "std": 18, "mono_fraction": 0.92,
            "di_mean": 2*150, "di_std": 30
        },
        nuclease_activities={
            "dnase1_activity": 1.4,   # Very high due to open chromatin
            "dnase1l3_activity": 0.9,  # Reduced due to lower methylation
            "dffb_activity": 2.0      # High apoptotic/necrotic activity
        },
        chromatin_beta_params=(2.0, 5.0)  # Very open, disrupted chromatin
    )
}


class TissueMixtureSimulator(FragmentSimulator):
    """
    cfDNA simulator for multi-tissue mixtures. This class extends the base
    FragmentSimulator.

    Examples
    --------
    Basic tissue mixture simulation:

    >>> seq_gen = SequenceContextGenerator(\"/path/to/reference.fasta\")
    >>> simulator = TissueMixtureSimulator(seq_gen)
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
    """
    def __init__(
        self, sequence_generator: SequenceContextGenerator,
        **kwargs: dict[str, object]
    ) -> None:
        super().__init__(sequence_generator, **kwargs)  # type: ignore
        self.tissue_profiles = TISSUE_PROFILES.copy()
        self.logger = logging.getLogger(LOGGER_NAME)

    def add_custom_tissue(self, tissue_profile: TissueProfile) -> None:
        """
        Add a custom tissue profile to the simulator's available profiles.

        This method allows users to define custom tissue types with specific
        fragmentation characteristics beyond the predefined profiles.

        Parameters
        ----------
        tissue_profile : TissueProfile
            Complete tissue profile definition including name, fragment size
            distribution, nucleosome spacing, chromatin openness, end motif
            preferences, and methylation level.
        """
        self.tissue_profiles[tissue_profile.name] = tissue_profile
        self.logger.info(f"Added custom tissue profile: {tissue_profile.name}")

    def simulate_tissue_mixture(
        self, tissue_types: list[str], tissue_fractions: list[float],
        total_fragments: int, genomic_regions: list[tuple[str, int, int]],
        add_noise: bool = True
    ) -> FragmentList:
        """
        Generate cfDNA mixture from multiple tissue sources.

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

            # @NOTE(ds): Fragments are _note_ sampled from regions proportional
            # to region length! This might be unexpected for certain types of
            # analysis (if e.g. many fragments come from short regions).
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

        Parameters
        ----------
        tissue : TissueProfile
            Tissue profile containing nuclease activity specifications.

        Returns
        -------
        NucleaseProfile
            Nuclease activity profile with tissue-specific activity levels
            and default sequence preferences.

        Note
        ----
        Future releases might change this API to be more customizable, i.e.
        to set different preferences on the nucleases, too.
        """
        return NucleaseProfile(
            dnase1_activity=tissue.nuclease_activities["dnase1_activity"],
            dnase1l3_activity=tissue.nuclease_activities["dnase1l3_activity"],
            dffb_activity=tissue.nuclease_activities["dffb_activity"]
        )

    def _simulate_tissue_fragments(
        self, chrom: str, start: int, end: int, num_fragments: int,
        tissue: TissueProfile, nuclease_profile: NucleaseProfile
    ) -> FragmentList:
        """
        Generate fragments with tissue-specific characteristics.

        This method creates fragments using a tissue-specific FragmentSimulator
        with custom nucleosome map based on the tissue's chromatin state.

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
        Creates tissue-specific nucleosome map using beta distribution
        parameters, then uses a dedicated simulator instance for this tissue.
        """
        alpha, beta = tissue.chromatin_beta_params
        tissue_nucleosome_map = FragmentSimulator.generate_nucleosome_map(
            self.sequence_generator, alpha, beta
        )
        tissue_simulator = FragmentSimulator(
            sequence_generator=self.sequence_generator,
            nucleosome_map=tissue_nucleosome_map
        )
        fragment_size_params = {
            "mean": tissue.fragment_size_distribution["mean"],
            "std": tissue.fragment_size_distribution["std"],
            "di_mean": tissue.fragment_size_distribution["di_mean"],
            "di_std": tissue.fragment_size_distribution["di_std"],
            "mono_fraction": tissue.fragment_size_distribution["mono_fraction"]
        }

        fragments = tissue_simulator.simulate_fragments(
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
        - Fragment Loss (5% random loss)
        - Size Variation (±2bp standard deviation)**:
        - End Repair Artifacts (2% of fragments)**:
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
        genomic_regions: list[tuple[str, int, int]]
    ) -> dict[str, FragmentList]:
        """
        Model longitudinal cancer progression.

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

        Examples
        --------

        >>> progression = simulator.simulate_cancer_progression(
        ...     normal_profile="hematopoietic",
        ...     tumor_fractions=[0.01, 0.03, 0.08, 0.15],
        ...     time_points=["diagnosis", "1month", "3months", "6months"],
        ...     fragments_per_timepoint=25000,
        ...     genomic_regions=[("chr1", 1000000, 2000000)]
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

            results[tp] = fragments

        return results
