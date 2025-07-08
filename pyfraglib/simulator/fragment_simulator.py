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
import numpy.typing as npt

from dataclasses import dataclass
from intervaltree import IntervalTree
from typing import Final

from pyfraglib.fragment import Fragment, FragmentList
from pyfraglib.core import get_chromosome_length, fail

LOGGER_NAME: Final[str] = "pyfraglib.simulator"


@dataclass
class NucleosomeMap:
    """Container for nucleosome positioning data."""
    positions: dict[str, npt.NDArray[np.int64]]  # chr -> nucleosome centers
    occupancy: dict[str, npt.NDArray[np.float64]]  # chr -> occupancy scores


@dataclass  # type: ignore
class ChromatinState:  # type: ignore
    """Represents chromatin accessibility state."""
    open_regions: IntervalTree  # type: ignore
    tf_binding_sites: IntervalTree  # type: ignore
    ctcf_sites: IntervalTree  # type: ignore


@dataclass
class NucleaseProfile:
    """Nuclease activity parameters."""
    dnase1_activity: float = 1.0
    dnase1l3_activity: float = 1.0
    dffb_activity: float = 1.0
    dnase1_motif_preference: dict[str, float] | None = None
    dnase1l3_motif_preference: dict[str, float] | None = None

    def __post_init__(self) -> None:
        if self.dnase1_motif_preference is None:
            self.dnase1_motif_preference = {
                "CG": 0.8, "GC": 0.8,  # Prefers naked DNA
                "AT": 1.0, "TA": 1.0,
                "AA": 1.1, "TT": 1.1
            }
        if self.dnase1l3_motif_preference is None:
            self.dnase1l3_motif_preference = {
                "CC": 1.5, "GG": 1.5,  # Strong CC preference
                "CG": 1.2, "GC": 1.2,
                "CT": 1.1, "TC": 1.1
            }


class FragmentSimulator:
    """
    Simulates cfDNA fragments based on biological parameters.

    This simulator models the stepwise fragmentation process:
    1. Intracellular digestion (DFFB-mediated, nucleosome-guided)
    2. Extracellular digestion (DNASE1/DNASE1L3-mediated)
    """
    MONO_NUC_MEAN: Final[int] = 167
    MONO_NUC_STD: Final[int] = 10
    DI_NUC_MEAN: Final[int] = 167*2
    DI_NUC_STD: Final[int] = 15
    LINKER_MEAN: Final[int] = 20
    LINKER_STD: Final[int] = 5
    PERIODICITY_AMPLITUDE: Final[float] = 0.1
    PERIODICITY_PHASE: Final[float] = 0.0
    MIN_FRAGMENT_SIZE: Final[int] = 40
    MAX_FRAGMENT_SIZE: Final[int] = 900

    def __init__(
        self, genome_ref: str = "hg19",
        nucleosome_map: NucleosomeMap | None = None,
        chromatin_state: ChromatinState | None = None
    ) -> None:
        """
        Initialize the fragment simulator.

        Args:
            genome_ref: Reference genome version ("hg19" or "hg38")
            nucleosome_map: Custom nucleosome positioning data
            chromatin_state: Chromatin accessibility information
        """
        self.logger: logging.Logger = logging.getLogger(LOGGER_NAME)
        self.genome_ref: str = genome_ref
        self.nucleosome_map: NucleosomeMap = (
            nucleosome_map or self._default_nucleosome_model()
        )
        self.chromatin_state: ChromatinState = (
            chromatin_state or self._default_chromatin_state()
        )

    def _default_nucleosome_model(self) -> NucleosomeMap:
        """
        Generate a default nucleosome map based on sequence features.
        Uses a simplified model with ~147bp nucleosomes and variable linkers.
        """
        self.logger.info("Generating default nucleosome model.")

        positions: dict[str, npt.NDArray[np.int64]] = {}
        occupancy: dict[str, npt.NDArray[np.float64]] = {}
        chromosomes: list[str] = [
            "1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
            "11", "12", "13", "14", "15", "16", "17", "18", "19",
            "20", "21", "22", "X", "Y"
        ]

        for chrom in chromosomes:
            chrom_len: int = get_chromosome_length(chrom, self.genome_ref)
            base_positions: npt.NDArray[np.int64] = np.arange(
                100, chrom_len - 100, 185  # 147bp nucleosome + ~38bp linker
            )
            noise: npt.NDArray[np.float64] = np.random.normal(
                0, 10, len(base_positions)
            )
            positions[chrom] = \
                base_positions + noise.astype(np.int64)  # type: ignore
            occupancy[chrom] = \
                np.random.beta(5, 2, len(base_positions))  # type: ignore

        return NucleosomeMap(positions=positions, occupancy=occupancy)

    def _default_chromatin_state(self) -> ChromatinState:
        """
        Generate default chromatin accessibility regions. More sophisticated
        simulations could load these values from ATAC-seq or another, similar
        data modality.
        """
        self.logger.info("Generating default chromatin state.")

        open_regions: IntervalTree = IntervalTree()  # type: ignore
        tf_sites: IntervalTree = IntervalTree()  # type: ignore
        ctcf_sites: IntervalTree = IntervalTree()  # type: ignore

        for _ in range(100):
            start: int = np.random.randint(1000, 10000000)
            end: int = start + np.random.randint(200, 1000)
            open_regions.addi(start, end)  # type: ignore

        return ChromatinState(
            open_regions=open_regions,  # type: ignore
            tf_binding_sites=tf_sites,  # type: ignore
            ctcf_sites=ctcf_sites  # type: ignore
        )

    def _get_cleavage_probability(
        self, position: int, chrom: str, nuclease_profile: NucleaseProfile
    ) -> float:
        """
        Calculate probability of cleavage at a given position.

        Considers:
        - Nucleosome occupancy
        - Chromatin accessibility
        - Nuclease preferences
        """
        nuc_positions: npt.NDArray[np.int64] = (
            self.nucleosome_map.positions.get(
                chrom, np.array([])  # type: ignore
            )
        )
        nuc_occupancy: npt.NDArray[np.float64] = (
            self.nucleosome_map.occupancy.get(
                chrom, np.array([])  # type: ignore
            )
        )

        base_prob: float
        if len(nuc_positions) > 0:
            distances: npt.NDArray[np.float64] = np.abs(
                nuc_positions - position  # type: ignore
            )
            nearest_idx: int = int(np.argmin(distances))
            distance_to_nuc: float = \
                float(distances[nearest_idx])  # type: ignore
            occupancy: float = \
                float(nuc_occupancy[nearest_idx])  # type: ignore

            # @NOTE(ds): Probability decreases near nucleosome centers.
            if distance_to_nuc < 73:  # within nucleosome core
                base_prob = 0.1 * (1 - occupancy)
            elif distance_to_nuc < 100:  # nucleosome edge
                base_prob = 0.3
            else:  # linker region
                base_prob = 0.8
        else:
            base_prob = 0.5

        if self.chromatin_state.open_regions.overlaps(  # type: ignore
            position
        ):
            base_prob *= 1.5

        total_activity: float = (
            nuclease_profile.dnase1_activity +
            nuclease_profile.dnase1l3_activity +
            nuclease_profile.dffb_activity
        )

        return min(base_prob * total_activity / 3.0, 1.0)

    def _generate_fragment_sizes(
        self,
        num_fragments: int,
        size_distribution: str = "normal",
        fragment_size_params: dict[str, float] | None = None
    ) -> npt.NDArray[np.int64]:
        """
        Generate fragment sizes based on biological distributions.

        Args:
            num_fragments: Number of fragments to generate
            size_distribution: Type of dist ("normal", "cancer", "fetal")
            fragment_size_params: Optional dict of means and stds to use

        Returns:
            Array of fragment sizes
        """
        if fragment_size_params:
            mono_mean = fragment_size_params["mean"]
            mono_std = fragment_size_params["std"]
            mono_fraction = 1.0 - fragment_size_params.get(
                "short_fraction", 0.25
            )
        else:
            mono_mean = self.MONO_NUC_MEAN
            mono_std = self.MONO_NUC_STD
            mono_fraction = 0.75

        sizes: npt.NDArray[np.float64]
        if size_distribution == "normal":
            num_mono: int = int(num_fragments * mono_fraction)
            num_di: int = num_fragments - num_mono

            mono_sizes: npt.NDArray[np.float64] = np.random.normal(
                mono_mean, mono_std, num_mono
            )
            di_sizes: npt.NDArray[np.float64] = np.random.normal(
                mono_mean * 2, mono_std * 1.5, num_di
            )
            sizes = np.concatenate([mono_sizes, di_sizes])  # type: ignore

        elif size_distribution == "cancer":
            mean_shift: int = -20
            mono_sizes = np.random.normal(
                mono_mean + mean_shift,
                mono_std * 1.2,
                int(num_fragments * 0.85)
            )
            di_sizes = np.random.normal(
                mono_mean * 2 + mean_shift,
                mono_std * 1.3,
                int(num_fragments * 0.15)
            )
            sizes = np.concatenate([mono_sizes, di_sizes])  # type: ignore

        elif size_distribution == "fetal":
            mean_shift = -25
            mono_sizes = np.random.normal(
                mono_mean + mean_shift,
                mono_std * 0.8,
                int(num_fragments * 0.9)
            )
            di_sizes = np.random.normal(
                mono_mean * 2 + mean_shift,
                mono_std,
                int(num_fragments * 0.1)
            )
            sizes = np.concatenate([mono_sizes, di_sizes])  # type: ignore

        else:
            fail(f"Unknown size distribution: {size_distribution}")

        sizes = self._add_periodicity(sizes)
        sizes_int: npt.NDArray[np.int64] = np.clip(
            sizes, self.MIN_FRAGMENT_SIZE, self.MAX_FRAGMENT_SIZE
        ).astype(np.int64)

        return sizes_int

    def _add_periodicity(
        self, sizes: npt.NDArray[np.float64]
    ) -> npt.NDArray[np.float64]:
        """Add 10bp periodicity to fragment sizes."""
        phase_factor: npt.NDArray[np.float64] = np.sin(
                2 * np.pi * sizes / 10 + self.PERIODICITY_PHASE  # type: ignore
        )
        modulation: npt.NDArray[np.float64] = (
            1 + self.PERIODICITY_AMPLITUDE * phase_factor  # type: ignore
        )
        return sizes * modulation  # type: ignore

    def _generate_end_motif(
        self, sequence_context: str, nuclease_profile: NucleaseProfile
    ) -> str:
        """
        Generate end motif based on nuclease preferences.

        Args:
            sequence_context: Surrounding sequence (at least 4bp)
            nuclease_profile: Active nucleases and their preferences

        Returns:
            4-mer end motif
        """
        if len(sequence_context) < 4:
            return "NNNN"

        motif: str = sequence_context[:4]
        if (nuclease_profile.dnase1l3_activity >
                nuclease_profile.dnase1_activity):
            if np.random.random() < 0.3:
                motif = "CC" + motif[2:]  # DNASE1L3 prefers CC motifs
        return motif

    def simulate_fragments(
        self, chrom: str, start: int, end: int, num_fragments: int,
        tissue_type: str = "healthy",
        nuclease_profile: NucleaseProfile | None = None,
        fragment_size_params: dict[str, float] | None = None
    ) -> FragmentList:
        """
        Simulate cfDNA fragments from a genomic region.

        Args:
            chrom: Chromosome name
            start: Start position
            end: End position
            num_fragments: Number of fragments to generate
            tissue_type: Source tissue type
            nuclease_profile: Nuclease activity parameters
            fragment_size_params: Custom size distribution parameters

        Returns:
            FragmentList containing simulated fragments
        """
        self.logger.info(
            f"Simulating {num_fragments} fragments from {chrom}:{start}-{end}."
        )

        if nuclease_profile is None:
            nuclease_profile = NucleaseProfile()

        tissue_to_size_dist: dict[str, str] = {
            "healthy": "normal",
            "hematopoietic": "normal",
            "liver": "normal",
            "placenta": "fetal",
            "tumor": "cancer",
            "fetal": "fetal"
        }
        size_dist: str = tissue_to_size_dist.get(tissue_type, "normal")
        sizes: npt.NDArray[np.int64] = self._generate_fragment_sizes(
            num_fragments, size_dist, fragment_size_params
        )

        max_size: int = int(sizes.max())  # type: ignore
        positions: npt.NDArray[np.int64] = np.random.randint(
            start, end - max_size, num_fragments
        )

        fragment_list: FragmentList = FragmentList()
        pos: int
        size: int
        for i, (pos, size) in enumerate(zip(positions, sizes)):  # type: ignore
            cleave_prob: float = self._get_cleavage_probability(
                pos, chrom, nuclease_profile
            )

            if np.random.random() > cleave_prob:
                continue  # fragment not released

            end5p: str = self._generate_end_motif(
                "ATCG" * 10, nuclease_profile
            )
            end3p: str = self._generate_end_motif(
                "CGAT" * 10, nuclease_profile
            )

            fragment: Fragment = Fragment.create_simulated(
                start_pos=int(pos),
                end_pos=int(pos + size),
                chrom=chrom,
                length=int(size),
                end5p=end5p,
                end3p=end3p,
                is_forward=bool(
                    np.random.choice([True, False])  # type: ignore
                ),
                is_mutated=False
            )

            fragment_list.append(fragment)

        self.logger.info(
            f"Generated {fragment_list.length()} fragments."
        )
        return fragment_list
