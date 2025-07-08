# This file is part of `pyfraglib`, a software suite to calculate fragmentomics
# features from cfDNA and perform downstream analyses.
#
# Copyright (C) 2024 Daniel Schütte, daniel.schuette@iccb-cologne.org
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version. This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
# more details. You should have received a copy of the GNU General Public
# License along with this program. If not, see <https://www.gnu.org/licenses/>.
import logging

import numpy as np

from dataclasses import dataclass
from typing import Final
from pyfraglib.fragment import FragmentList
from pyfraglib.simulator import FragmentSimulator
from pyfraglib.simulator.fragment_simulator import NucleaseProfile

LOGGER_NAME: Final[str] = "pyfraglib.simulator.mixture"


@dataclass
class TissueProfile:
    """Defines tissue-specific fragmentation characteristics"""
    name: str
    fragment_size_distribution: dict[str, float]
    nucleosome_spacing: float
    chromatin_openness: float
    end_motif_preferences: dict[str, float]
    methylation_level: float  # average methylation, affects DNASE1L3


# @NOTE(ds): The following tissue profiles are literature-based and not
# biologically validated!
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
    "placenta": TissueProfile(
        name="placenta",
        fragment_size_distribution={
            "mean": 143, "std": 15, "short_fraction": 0.3
        },
        nucleosome_spacing=180,
        chromatin_openness=0.6,
        end_motif_preferences={"CG": 1.4, "GC": 1.2},
        methylation_level=0.4  # hypomethylated
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
    """Disease-specific fragmentation alterations"""
    name: str
    size_shift: float
    periodicity_change: float
    preferred_cut_sites: list[str]
    aberrant_ends: dict[str, float]


class TissueMixtureSimulator(FragmentSimulator):
    """
    Extends `FragmentSimulator` to handle multiple tissue sources and disease
    states.
    """

    def __init__(self, *args, **kwargs):  # type: ignore
        super().__init__(*args, **kwargs)  # type: ignore
        self.tissue_profiles = TISSUE_PROFILES.copy()
        self.logger = logging.getLogger(LOGGER_NAME)

    def add_custom_tissue(self, tissue_profile: TissueProfile):  # type: ignore
        """Add a custom tissue profile to the simulator."""
        self.tissue_profiles[tissue_profile.name] = tissue_profile
        self.logger.info(f"Added custom tissue profile: {tissue_profile.name}")

    def simulate_tissue_mixture(
        self, tissue_types: list[str], tissue_fractions: list[float],
        total_fragments: int, genomic_regions: list[tuple[str, int, int]],
        add_noise: bool = True
    ) -> FragmentList:
        """
        Simulate cfDNA from multiple tissue sources.

        Args:
            tissue_types: List of tissue names
            tissue_fractions: Fraction from each tissue (must sum to 1)
            total_fragments: Total number of fragments to generate
            genomic_regions: List of (chrom, start, end) tuples to sample from
            add_noise: Add biological/technical noise

        Returns:
            Mixed FragmentList with tissue annotations
        """
        if abs(sum(tissue_fractions) - 1.0) > 0.001:
            raise ValueError("Tissue fractions must sum to 1.0")

        if len(tissue_types) != len(tissue_fractions):
            raise ValueError("Tissue number must match number of fractions")

        sim_info: dict[str, object] = dict(zip(tissue_types, tissue_fractions))
        self.logger.info(f"Simulating mixture: {sim_info}")

        all_fragments = FragmentList()
        for tissue, fraction in zip(tissue_types, tissue_fractions):
            if tissue not in self.tissue_profiles:
                raise ValueError(f"Unknown tissue type: {tissue}")

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
                    fragment._tissue_origin = tissue  # type: ignore
                all_fragments.extend(tissue_fragments)  # type: ignore

        if add_noise:
            all_fragments = self._add_biological_noise(all_fragments)

        return all_fragments

    def _get_tissue_nuclease_profile(
        self, tissue: TissueProfile
    ) -> NucleaseProfile:
        """Generate tissue-specific nuclease activity profile."""
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
        """Generate fragments with tissue-specific characteristics."""
        old_mono_mean = self.mono_nuc_mean  # type: ignore
        old_mono_std = self.mono_nuc_std  # type: ignore

        self.mono_nuc_mean = tissue.fragment_size_distribution["mean"]
        self.mono_nuc_std = tissue.fragment_size_distribution["std"]

        fragments = self.simulate_fragments(
            chrom, start, end,
            num_fragments,
            tissue_type=tissue.name,
            nuclease_profile=nuclease_profile
        )

        self.mono_nuc_mean = old_mono_mean
        self.mono_nuc_std = old_mono_std

        return fragments

    def _add_biological_noise(self, fragments: FragmentList) -> FragmentList:
        """Add realistic biological and technical noise to fragments."""
        noisy_fragments = FragmentList()

        for fragment in fragments:
            # @NOTE(ds): Random 5% fragment loss (inefficient extraction,
            # degradation, etc.).
            if np.random.random() < 0.05:
                continue

            size_noise = np.random.normal(0, 2)
            fragment.length = int(max(50, fragment.length + size_noise))
            fragment.end_pos = fragment.start_pos + fragment.length

            # @NOTE(ds): 2% of fragments suffer from end repair artifacts.
            if np.random.random() < 0.02:
                fragment.end5p = fragment.end5p[:-1] + "N"

            noisy_fragments.append(fragment)

        return noisy_fragments

    def simulate_cancer_progression(
        self, normal_profile: str, tumor_fractions: list[float],
        time_points: list[str], fragments_per_timepoint: int,
        genomic_regions: list[tuple[str, int, int]],
        cancer_type: str = "generic"
    ) -> dict[str, FragmentList]:
        """
        Simulate cfDNA across cancer progression time points.

        Args:
            normal_profile: Normal tissue type (e.g., "hematopoietic")
            tumor_fractions: Tumor fraction at each time point
            time_points: Time point labels
            fragments_per_timepoint: Fragments to generate per time point
            genomic_regions: Regions to sample from
            cancer_type: Specific cancer type for signatures

        Returns:
            Dictionary mapping time points to FragmentLists
        """
        if len(tumor_fractions) != len(time_points):
            raise ValueError("Tumor fractions must match time points")

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
                    fragments, cancer_type, tf
                )

            results[tp] = fragments

        return results

    def _add_cancer_signatures(
        self, fragments: FragmentList, cancer_type: str, tumor_fraction: float
    ) -> FragmentList:
        """Add cancer-specific fragmentation signatures."""
        cancer_signatures = {
            "breast": {"size_shift": -15, "short_enrichment": 0.2},
            "lung": {"size_shift": -20, "short_enrichment": 0.25},
            "colorectal": {"size_shift": -10, "short_enrichment": 0.15},
            "generic": {"size_shift": -15, "short_enrichment": 0.2}
        }

        sig = cancer_signatures.get(cancer_type, cancer_signatures["generic"])

        modified_fragments = FragmentList()
        for fragment in fragments:
            if (hasattr(fragment, "_tissue_origin") and  # type: ignore
                    fragment._tissue_origin == "tumor"):  # type: ignore
                new_length = int(fragment.length + sig["size_shift"]
                                 * np.random.random())
                new_length = max(50, new_length)
                fragment.length = new_length
                fragment.end_pos = fragment.start_pos + new_length
                if (fragment.length < 150 and
                        np.random.random() < sig["short_enrichment"]):
                    modified_fragments.append(fragment)

            modified_fragments.append(fragment)

        return modified_fragments

    def simulate_fetal_fraction(
        self, fetal_fraction: float, total_fragments: int,
        genomic_regions: list[tuple[str, int, int]], gestational_age: int = 20
    ) -> FragmentList:
        """
        Simulate maternal plasma cfDNA with fetal fraction (for NIPT).

        Args:
            fetal_fraction: Fraction of fetal DNA (0-1)
            total_fragments: Total fragments to generate
            genomic_regions: Regions to sample
            gestational_age: Weeks of gestation (affects fragment size)

        Returns:
            Mixed maternal and fetal fragments
        """
        fetal_profile = TissueProfile(
            name="fetal",
            fragment_size_distribution={
                "mean": 143 - (20 - gestational_age) * 0.5,  # size up with age
                "std": 12,
                "short_fraction": 0.35 - (gestational_age - 10) * 0.01
            },
            nucleosome_spacing=175,
            chromatin_openness=0.7,
            end_motif_preferences={"CG": 1.3, "GC": 1.3},
            methylation_level=0.3  # hypomethylated
        )

        self.add_custom_tissue(fetal_profile)
        return self.simulate_tissue_mixture(
            ["hematopoietic", "fetal"],
            [1 - fetal_fraction, fetal_fraction],
            total_fragments,
            genomic_regions
        )
