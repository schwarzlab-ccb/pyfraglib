"""
Biologically Realistic cfDNA Fragment Simulation
================================================

This module provides simulation capabilities for cell-free DNA (cfDNA)
fragments based on experimentally observed fragmentation biology. The
simulator tries to achieve realism by integrating multiple layers of biological
complexity. It might be used e.g. for method validation or training data
generation.

Key Features
------------
- **Real genomic sequences**: Uses reference FASTA files to provide sequence
  context
- **Nucleosome positioning**: Models DNA protection by nucleosome
  through incorporating nucleosome maps
- **Chromatin accessibility**: Tissue-specific chromatin states and
  accessibility, DNA shielding by transcription factors
- **Nuclease specificity**: Experimentally-determined nuclease cleavage
  preferences
- **Fragment size distributions**: Adjustable size patterns with periodicity
- **End motif generation**: Sequence-context aware end motif simulation

Biological Foundation
---------------------
Core parts of the simulator are based on the following assumptions:

1. **Nuclease Activities**:

    - DNASE1L3: plasma nuclease, CC dinucleotide preference
      (Serpas et al. 2019 PNAS)
    - DNASE1: endonuclease, AT-rich accessibility preference, predominantly
      produces T ends
      (Lazarovici et al. 2013 PNAS, Han et al. 2020 Am J Hum Genet.)
    - DFFB: Nucleosome linker preference, produces A ends
      (Han et al. 2020 Am J Hum Genet.)

2. **Fragment Size Patterns**:

   - Mono-nucleosomal peak (~167 bp) with 10 bp periodicity
   - Di-nucleosomal components (~334 bp)
     (reviewed in Lo et al. 2021 Science)

3. **Chromatin Context**:

   - Nucleosome positioning affects cleavage accessibility
   - Transcription factor binding site protection
     (reviewed in Lo et al. 2021 Science)

Example Usage
-------------
.. code-block:: python

    from pyfraglib.simulator.fragment_simulator import (
        FragmentSimulator, NucleaseProfile, SequenceContextGenerator
    )

    # Initialize sequence generator with reference genome
    seq_gen = SequenceContextGenerator("/path/to/reference.fasta")

    # Initialize simulator with sequence generator
    simulator = FragmentSimulator(seq_gen)

    # Define nuclease activity profile
    nuclease_profile = NucleaseProfile(
        dnase1_activity=1.0,
        dnase1l3_activity=1.2,
        dffb_activity=0.5
    )

    # Simulate fragments from a genomic region
    fragments = simulator.simulate_fragments(
        chrom="1",
        start=1000000,
        end=1100000,
        num_fragments=10000,
        tissue_type="hematopoietic",
        nuclease_profile=nuclease_profile
    )

    # Process generated fragments
    print(f"Generated {fragments.length()} cfDNA fragments")

    # Clean up
    seq_gen.close()

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
import pysam

import numpy as np
import numpy.typing as npt

from dataclasses import dataclass
from functools import lru_cache
from intervaltree import IntervalTree
from typing import Final
from pyfraglib.fragment import Fragment, FragmentList
from pyfraglib.core import fail

LOGGER_NAME: Final[str] = "pyfraglib.simulator"


@dataclass
class NucleosomeMap:
    """
    Container for nucleosome positioning data across chromosomes.

    This class stores nucleosome positioning information, e.g. derived from
    experimental data. The positioning is used to model realistic nucleosome
    protection effects during cfDNA fragmentation.

    Notes
    -----
    Nucleosome occupancy scores affect the probability of DNA cleavage:

    - High occupancy (>0.8): Strong protection, low cleavage probability
    - Medium occupancy (0.3-0.8): Moderate protection
    - Low occupancy (<0.3): Weak protection, higher cleavage probability

    The positions array should be sorted in ascending order for efficient
    binary search operations during simulation.

    Examples
    --------
    >>> import numpy as np
    >>> nuc_map = NucleosomeMap(
    ...     positions={"chr1": np.array([1000, 1185, 1370], dtype=np.int64)},
    ...     occupancy={"chr1": np.array([0.8, 0.6, 0.9], dtype=np.float64)}
    ... )
    >>> print(f"Nucleosome at position {nuc_map.positions['chr1'][0]}")
    Nucleosome at position 1000
    """
    positions: dict[str, npt.NDArray[np.int64]]  # chr -> nucleosome centers
    occupancy: dict[str, npt.NDArray[np.float64]]  # chr -> occupancy scores


@dataclass  # type: ignore
class ChromatinState:  # type: ignore
    """
    Represents tissue-specific chromatin accessibility and regulatory states.

    This class encapsulates chromatin accessibility information that affects
    cfDNA fragmentation patterns. Different chromatin states lead to varying
    nuclease accessibility and cleavage probability.

    Notes
    -----
    Each IntervalTree should contain intervals with genomic coordinates:

    - Chromosome-specific trees for efficient spatial queries
    - 0-based coordinate system following standard conventions
    - Additional metadata can be stored as interval data

    The chromatin state significantly influences:

    - Base cleavage probability (open regions: ~0.6, closed regions: ~0.1)
    - Transcription factor protection effects (30% of normal cleavage)
    - Tissue-specific fragmentation patterns

    Examples
    --------
    >>> from intervaltree import IntervalTree
    >>> chr1_open = IntervalTree()
    >>> chr1_open.addi(1000, 2000)  # Open chromatin region
    >>> chr1_tf = IntervalTree()
    >>> chr1_tf.addi(1500, 1520)  # Protected TF binding site
    >>> chromatin = ChromatinState(
    ...     open_regions={"chr1": chr1_open},
    ...     tf_binding_sites={"chr1": chr1_tf},
    ...     ctcf_sites={"chr1": IntervalTree()}
    ... )
    """
    #: Genomic intervals representing accessible chromatin regions (fragments
    #: originating from these regions have higher cleavage probability).
    open_regions: dict[str, IntervalTree]  # type: ignore

    #: Transcription factor binding sites that provide protection from
    #: nuclease cleavage. These regions typically show reduced fragmentation.
    tf_binding_sites: dict[str, IntervalTree]  # type: ignore

    #: CTCF binding sites that create specific chromatin architecture and
    #: affect local nuclease accessibility patterns.
    ctcf_sites: dict[str, IntervalTree]  # type: ignore


@dataclass
class NucleaseProfile:
    """
    Nuclease activity and sequence preference parameters.

    This class defines the activity levels and sequence preferences of 3
    different nucleases involved in cfDNA fragmentation. The parameters are
    based on experimental literature and can be customized for different
    biological conditions or disease states.

    Notes
    -----
    Activity levels are relative and multiplicative:

    - 1.0 represents normal physiological activity
    - Values >1.0 increase nuclease activity
    - Values <1.0 decrease nuclease activity
    - 0.0 completely disables the nuclease

    Motif preferences are multiplicative factors applied to base cleavage
    probability based on local sequence context. Higher values increase
    cleavage probability for sequences containing the motif.

    Examples
    --------
    >>> # Healthy individual profile
    >>> healthy_profile = NucleaseProfile(
    ...     dnase1_activity=1.0,
    ...     dnase1l3_activity=1.2,
    ...     dffb_activity=0.1
    ... )
    >>>
    >>> # Cancer patient with increased apoptosis
    >>> cancer_profile = NucleaseProfile(
    ...     dnase1_activity=1.1,
    ...     dnase1l3_activity=1.0,
    ...     dffb_activity=2.0
    ... )
    """
    #: Relative activity level of DNase I. DNase I shows preference for AT-rich
    #: accessible regions.
    dnase1_activity: float = 1.0

    #: Relative activity level of DNase I Like 3 (DNase1L3). This is the
    #: primary nuclease responsible for cfDNA fragmentation in healthy
    #: individuals. Shows strong CC dinucleotide preference.
    dnase1l3_activity: float = 1.0

    #: Relative activity level of the apoptotic nuclease DNA Fragmentation
    #: Factor B (DFFB). More active during cell death and shows preference
    #: for nucleosome linker regions.
    dffb_activity: float = 1.0

    #: Sequence motif preferences for DNase I. If None, uses default
    #: preferences based on literature (AT-rich sequences preferred).
    dnase1_motif_preference: dict[str, float] | None = None

    #: Sequence motif preferences for DNase1L3. If None, uses default
    #: preferences (strong CC dinucleotide preference).
    dnase1l3_motif_preference: dict[str, float] | None = None

    #: Sequence motif preferences for DFFB. If None, uses default
    #: preferences (slight A preference).
    dffb_motif_preference: dict[str, float] | None = None

    def __post_init__(self) -> None:
        if not (0.0 <= self.dnase1_activity):
            fail(f"DNase I activity < 0 ({self.dnase1_activity})")
        if not (0.0 <= self.dnase1l3_activity):
            fail(f"DNase1L3 activity < 0 ({self.dnase1l3_activity})")
        if not (0.0 <= self.dffb_activity):
            fail(f"DFFB activity < 0 ({self.dffb_activity})")

        # Validate that at least one nuclease is active
        total_activity = (self.dnase1_activity +
                          self.dnase1l3_activity +
                          self.dffb_activity)
        if total_activity <= 0.01:
            fail("At least one nuclease must have activity > 0.01")

        if self.dnase1_motif_preference is None:
            self.dnase1_motif_preference = {
                "CG": 0.8, "GC": 0.8,  # Lower preference for CpG sites
                "AT": 1.2, "TA": 1.2,  # Clear AT preference
                "AA": 1.1, "TT": 1.1   # Slight preference for AT-rich
            }
        if self.dnase1l3_motif_preference is None:
            self.dnase1l3_motif_preference = {
                "CC": 1.5,             # Strong CC preference
                "GG": 1.2,             # Moderate GG preference
                "CG": 1.2, "GC": 1.2,  # Moderate CpG preference
                "CT": 1.1, "TC": 1.1   # Slight pyrimidine preference
            }
        if self.dffb_motif_preference is None:
            self.dffb_motif_preference = {
                "A": 1.1,              # Slight A preference
                "T": 1.05,             # Slight T preference
                "C": 0.95, "G": 0.95   # Slightly reduced C/G preference
            }


class SequenceContextGenerator:
    """
    Genomic sequence context provider for cfDNA fragment simulation.

    This class interfaces with reference FASTA files to provide authentic
    genomic sequence context for fragment end motif generation. The class
    implements sequence caching and provides methods for extracting sequence
    context around genomic positions, enabling simulation of cfDNA fragments
    with somewhat realistic end motifs.

    Notes
    -----
    The FASTA file must be indexed with samtools faidx before use.

    Examples
    --------
    >>> seq_gen = SequenceContextGenerator(\"/path/to/reference.fasta\")
    >>> context = seq_gen.generate_sequence_context(\"chr1\", 1000000, 20)
    >>> print(f\"Sequence context: {context}\")
    >>> seq_gen.close()

    Using as context manager:

    >>> with SequenceContextGenerator(\"/path/to/reference.fasta\") as seq_gen:
    ...     context = seq_gen.generate_sequence_context(\"chr1\", 1000000, 20)
    """
    def __init__(self, fasta_path: str) -> None:
        """
        Initialize sequence generator with reference FASTA file.

        Args:
            fasta_path: Path to reference FASTA file (must be indexed)
        """
        self.logger = logging.getLogger(LOGGER_NAME + ".sequence")
        self.fasta_path = fasta_path

        try:
            self.fasta = pysam.FastaFile(fasta_path)
        except Exception as e:
            fail(f"Could not open FASTA file '{fasta_path}': {e}")

        self.logger.info(f"Loaded reference FASTA: {fasta_path}")

    def generate_sequence_context(
        self, chrom: str, position: int, length: int = 20
    ) -> str:
        """
        Extract genomic sequence context around a specified position.

        This method retrieves genomic sequence from the reference FASTA file,
        centered around the specified position. It handles edge cases at
        chromosome boundaries by padding with 'N' nucleotides when necessary.

        Parameters
        ----------
        chrom : str
            Chromosome name (e.g., "chr1", "1", "X"). Must match chromosome
            naming convention in the FASTA file.
        position : int
            Genomic position (0-based coordinate system) around which to
            extract sequence context.
        length : int, default=20
            Total length of sequence to extract, centered on the position.
            Must be positive.

        Returns
        -------
        str
            DNA sequence string of specified length, uppercase nucleotides.
            Contains 'N' characters for undefined or unavailable regions.

        Raises
        ------
        Exception
            If chromosome is not found in FASTA file or sequence extraction
            fails. In such cases, returns a string of 'N' characters as
            fallback.

        Notes
        -----
        The method extracts sequence symmetrically around the position:
        - Start position: position - length//2
        - End position: start + length

        Boundary handling:
        - If start < 0: sequence is extracted from position 0
        - If end > chromosome length: sequence is padded with 'N'

        Examples
        --------
        >>> seq_gen = SequenceContextGenerator("reference.fasta")
        >>> context = seq_gen.generate_sequence_context("chr1", 1000000, 20)
        >>> print(f"Context: {context}")  # e.g., "ATCGATCGATCGATCGATCG"
        >>> len(context)
        """
        start = max(0, position - length // 2)
        end = start + length

        try:
            sequence = self.fasta.fetch(chrom, start, end)
            if len(sequence) < length:
                padding_needed = length - len(sequence)
                sequence = sequence + "N" * padding_needed
            return sequence.upper()

        except Exception as e:
            self.logger.warning(
                f"Could not extract sequence from {chrom}:{start}-{end}: {e}"
            )
            return "N" * length

    def get_chromosome_length(self, chrom: str) -> int:
        """
        Retrieve chromosome length from the reference FASTA file.

        Parameters
        ----------
        chrom : str
            Chromosome name to query. Must match naming convention in FASTA.

        Returns
        -------
        int
            Length of the chromosome in base pairs. Returns 0 if chromosome
            is not found or an error occurs.
        """
        try:
            return self.fasta.get_reference_length(chrom)
        except Exception:
            return 0

    def get_available_chromosomes(self) -> list[str]:
        """
        List all available chromosome references in the FASTA file.

        Returns
        -------
        list[str]
            List of chromosome names available in the FASTA file, in the
            order they appear in the file header.
        """
        return list(self.fasta.references)

    def close(self) -> None:
        """
        Close the FASTA file handle and release resources.
        """
        if hasattr(self, "fasta"):
            self.fasta.close()

    def __del__(self) -> None:
        """
        Destructor to ensure FASTA file is properly closed.
        """
        self.close()


class FragmentSimulator:
    """
    cfDNA fragment simulator.

    This simulator integrates experimental observations to model cfDNA
    fragmentation:

    1. **Nucleosome positioning**: Nucleosome maps with occupancy scores
    2. **Chromatin accessibility**: Tissue-specific chromatin states
    3. **Nuclease specificity**: Experimentally-determined preferences
    4. **Sequence context**: Genomic sequences from reference FASTA
    5. **TF binding protection**: Transcription factor binding site shielding

    Main references are Lo et al. 2021 Science, Cristiano et al. 2019 Nature.
    """
    MONO_NUC_MEAN: Final[int] = 167
    MONO_NUC_STD: Final[int] = 10
    MONO_FRACTION: Final[float] = 0.85
    DI_NUC_MEAN: Final[int] = 167*2
    DI_NUC_STD: Final[int] = 15
    PERIODICITY_AMPLITUDE: Final[float] = 0.1
    MIN_FRAGMENT_SIZE: Final[int] = 40
    MAX_FRAGMENT_SIZE: Final[int] = 900

    _prob_cache: dict[tuple[str, int], float]

    def __init__(
        self, sequence_generator: SequenceContextGenerator,
        nucleosome_map: NucleosomeMap | None = None,
        chromatin_state: ChromatinState | None = None
    ) -> None:
        """
        Initialize the fragment simulator.

        Args:
            sequence_generator: Sequence context generator (must be constructed
                with indexed FASTA file)
            nucleosome_map: Custom nucleosome positioning data. If None,
                generates default nucleosome map.
            chromatin_state: Chromatin accessibility information. If None,
                generates default chromatin state.
        """
        self.logger: logging.Logger = logging.getLogger(LOGGER_NAME)
        self.sequence_generator: SequenceContextGenerator = sequence_generator
        self.nucleosome_map: NucleosomeMap = (
            nucleosome_map or self._default_nucleosome_model()
        )
        self.chromatin_state: ChromatinState = (
            chromatin_state or self._default_chromatin_state()
        )

        self._prob_cache: dict[tuple[str, int], float] = {}

    @lru_cache(maxsize=8192)  # type: ignore
    def _get_cached_sequence_context(self, chrom: str, position: int) -> str:
        """
        Get sequence context with LRU caching for performance.

        Fetches 1200bp sequence context centered around genomic position.

        Args:
            chrom: Chromosome name
            position: Genomic position (rounded to nearest 1kb for caching)

        Returns:
            1200bp sequence context string
        """
        rounded_pos = (position // 1000) * 1000
        try:
            return self.sequence_generator.generate_sequence_context(
                chrom, rounded_pos + 600, length=1200
            )
        except Exception:
            return "N" * 1200

    @staticmethod
    def generate_nucleosome_map(
        sequence_generator: SequenceContextGenerator,
        beta_alpha: float = 4.0,
        beta_beta: float = 3.0
    ) -> NucleosomeMap:
        """
        Generate a nucleosome map with customizable chromatin state.

        This static method creates a nucleosome map using beta distribution
        parameters to control chromatin accessibility. Can be used by both
        the base simulator and tissue-specific simulators.

        Parameters
        ----------
        sequence_generator : SequenceContextGenerator
            Sequence generator to get chromosome information.
        beta_alpha : float, default=4.0
            Alpha parameter for beta distribution. Higher values increase
            nucleosome occupancy (more compact chromatin).
        beta_beta : float, default=3.0
            Beta parameter for beta distribution. Higher values decrease
            nucleosome occupancy (more open chromatin).

        Returns
        -------
        NucleosomeMap
            Nucleosome map with specified chromatin characteristics.

        Notes
        -----
        Common beta parameter combinations:
        - Open chromatin: (2, 5) - low occupancy, high accessibility
        - Normal chromatin: (4, 3) - balanced occupancy
        - Compact chromatin: (6, 2) - high occupancy, low accessibility
        """
        positions: dict[str, npt.NDArray[np.int64]] = {}
        occupancy: dict[str, npt.NDArray[np.float64]] = {}

        available_chroms = sequence_generator.get_available_chromosomes()
        for chrom in available_chroms:
            chrom_len: int = sequence_generator.get_chromosome_length(chrom)
            if chrom_len < 10_000:  # Skip very short contigs
                continue

            base_positions: npt.NDArray[np.int64] = np.arange(
                100, chrom_len - 100, 185  # 147 bp nucleosome + 38 bp linker
            )
            noise: npt.NDArray[np.float64] = np.random.normal(
                0, 10, len(base_positions)
            )
            positions[chrom] = \
                base_positions + noise.astype(np.int64)  # type: ignore
            occupancy[chrom] = np.random.beta(  # type: ignore
                beta_alpha, beta_beta, len(base_positions)
            )

        return NucleosomeMap(positions=positions, occupancy=occupancy)

    def _default_nucleosome_model(self) -> NucleosomeMap:
        """
        Generate a default nucleosome map using standard parameters.
        """
        self.logger.info("Generating default nucleosome model.")
        return self.generate_nucleosome_map(self.sequence_generator)

    def _default_chromatin_state(self) -> ChromatinState:
        """
        Generate default chromatin accessibility regions. More sophisticated
        simulations could load these values from ATAC-seq or another, similar
        data modality.
        """
        self.logger.info("Generating default chromatin state.")

        open_regions_dict: dict[str, IntervalTree] = {}  # type: ignore
        tf_sites_dict: dict[str, IntervalTree] = {}  # type: ignore
        ctcf_sites_dict: dict[str, IntervalTree] = {}  # type: ignore

        available_chroms = self.sequence_generator.get_available_chromosomes()
        for chrom in available_chroms:
            chrom_length = self.sequence_generator.get_chromosome_length(chrom)
            if chrom_length == 0:
                continue

            chrom_open = IntervalTree()  # type: ignore
            for _ in range(1000):
                start: int = np.random.randint(
                    1000, min(chrom_length - 1000, 10000000)
                )
                end: int = start + np.random.randint(200, 1000)
                if end < chrom_length:
                    chrom_open.addi(start, end)  # type: ignore

            open_regions_dict[chrom] = chrom_open  # type: ignore
            tf_sites_dict[chrom] = IntervalTree()  # type: ignore
            ctcf_sites_dict[chrom] = IntervalTree()  # type: ignore

        return ChromatinState(
            open_regions=open_regions_dict,  # type: ignore
            tf_binding_sites=tf_sites_dict,  # type: ignore
            ctcf_sites=ctcf_sites_dict  # type: ignore
        )

    def _get_cleavage_probability(
        self, position: int, chrom: str, nuclease_profile: NucleaseProfile,
        tissue_factor: float
    ) -> float:
        """
        Calculate realistic cleavage probability at a genomic position.

        Uses a multiplicative model integrating multiple biological factors.

        Args:
            position: Genomic position (0-based)
            chrom: Chromosome name
            nuclease_profile: Active nucleases and their preferences
            tissue_factor: Tissue-specific cleavage factor

        Returns:
            Cleavage probability between 0.001 and 1.0
        """
        base_prob = 0.1  # default is closed chromatin
        if chrom in self.chromatin_state.open_regions:  # type: ignore
            overlaps: list[object] = \
                self.chromatin_state.open_regions[  # type: ignore
                    chrom
                ].at(position)
            if overlaps:
                base_prob = 0.6

        nucleosome_factor = self._get_nucleosome_protection(position, chrom)
        nuclease_factor = self._get_nuclease_activity_factor(
            position, chrom, nuclease_profile
        )
        tf_protection = self._get_tf_protection_factor(position, chrom)

        total_prob = (base_prob * nucleosome_factor * tissue_factor *
                      nuclease_factor * tf_protection)

        return min(max(total_prob, 0.001), 1.0)

    def _get_nucleosome_protection(self, position: int, chrom: str) -> float:
        """
        Calculate nucleosome protection factor at a genomic position.

        Optimized with caching.
        """
        if chrom not in self.nucleosome_map.positions:
            return 1.0  # no protection

        actual_distance, cached_occupancy = \
            self._get_nearest_nucleosome_info(chrom, position)

        # @NOTE(ds): Nucleosome protection: LOWER values = HIGHER protection
        # (these are multiplicative factors, so <1.0 reduces cleavage
        # probability).
        if actual_distance <= 73:  # nucleosome core
            return 0.05 + (1.0 - cached_occupancy) * 0.15
        elif actual_distance <= 120:
            return 0.20 + (1.0 - cached_occupancy) * 0.30
        else:  # linker region
            return 0.70 + (1.0 - cached_occupancy) * 0.30

    @lru_cache(maxsize=16384)  # type: ignore
    def _get_nearest_nucleosome_info(
        self, chrom: str, position: int
    ) -> tuple[float, float]:
        """
        Get distance and occupancy for nearest nucleosome to genomic position.

        Uses LRU caching with 100bp resolution for performance.

        Returns:
            Tuple of (distance_to_nearest_nucleosome, occupancy_score)
        """
        cache_pos = (position // 100) * 100

        nucleosome_positions = self.nucleosome_map.positions[chrom]
        nucleosome_occupancy = self.nucleosome_map.occupancy[chrom]

        nearest_idx = np.searchsorted(nucleosome_positions, cache_pos)
        candidates = []
        if nearest_idx > 0:
            candidates.append(nearest_idx - 1)
        if nearest_idx < len(nucleosome_positions):
            candidates.append(nearest_idx)

        if not candidates:
            return 100.0, 0.5

        distances = np.abs(  # type: ignore
            nucleosome_positions[candidates] - cache_pos  # type: ignore
        )
        best_candidate_idx = np.argmin(distances)  # type: ignore
        best_idx = candidates[best_candidate_idx]

        cached_distance: float = \
            float(distances[best_candidate_idx])  # type: ignore
        cached_occupancy: float = \
            float(nucleosome_occupancy[best_idx])  # type: ignore
        actual_distance = abs(cached_distance + (position - cache_pos))

        return actual_distance, cached_occupancy

    def _get_tissue_accessibility_factor(self, tissue_type: str) -> float:
        """
        Get tissue-specific chromatin accessibility factor.
        """
        tissue_factors = {
            "healthy": 1.0,
            "hematopoietic": 1.2,  # More open chromatin
            "liver": 0.9,          # Slightly more compact
            "tumor": 1.4           # Highly accessible
        }
        return tissue_factors.get(tissue_type, 1.0)

    def _get_nuclease_activity_factor(
        self, position: int, chrom: str, nuclease_profile: NucleaseProfile
    ) -> float:
        """
        Calculate nuclease-specific activity factor based on sequence
        preferences.

        Uses configured motif preferences from nuclease profile to determine
        cleavage probability at genomic position based on local sequence
        context.
        """
        context = self._get_cached_sequence_context(chrom, position)
        cache_key = (chrom, position // 1000 * 1000)
        offset = position - cache_key[1]
        if offset >= 10 and offset < len(context) - 10:
            local_context = context[offset-10:offset+10]
        else:
            if offset < 10:
                local_context = context[:20] if len(context) >= 20 else context
            else:
                local_context = \
                    context[-20:] if len(context) >= 20 else context

        if len(local_context) == 0:
            return 1.0  # no sequence context available

        total_nuclease_activity = (nuclease_profile.dnase1_activity +
                                   nuclease_profile.dnase1l3_activity +
                                   nuclease_profile.dffb_activity)
        if total_nuclease_activity == 0:
            return 1.0

        nuclease_factors = []
        if (nuclease_profile.dnase1_activity > 0 and
                nuclease_profile.dnase1_motif_preference):
            dnase1_factor = self._calculate_sequence_preference_factor(
                local_context, nuclease_profile.dnase1_motif_preference
            )
            nuclease_factors.append(
                (dnase1_factor, nuclease_profile.dnase1_activity)
            )

        if (nuclease_profile.dnase1l3_activity > 0 and
                nuclease_profile.dnase1l3_motif_preference):
            dnase1l3_factor = self._calculate_sequence_preference_factor(
                local_context, nuclease_profile.dnase1l3_motif_preference
            )
            nuclease_factors.append(
                (dnase1l3_factor, nuclease_profile.dnase1l3_activity)
            )

        if nuclease_profile.dffb_activity > 0:
            if nuclease_profile.dffb_motif_preference:
                dffb_sequence_factor = \
                    self._calculate_sequence_preference_factor(
                        local_context, nuclease_profile.dffb_motif_preference
                    )
            else:
                dffb_sequence_factor = 1.0

            nucleosome_protection = self._get_nucleosome_protection(
                position, chrom
            )
            dffb_linker_preference = 0.5 + 1.5 * nucleosome_protection
            dffb_factor = dffb_sequence_factor * dffb_linker_preference
            nuclease_factors.append(
                (dffb_factor, nuclease_profile.dffb_activity)
            )

        if not nuclease_factors:
            return 1.0

        weighted_activity = sum(
            factor * activity for factor, activity in nuclease_factors
        )
        return weighted_activity / total_nuclease_activity

    def _calculate_sequence_preference_factor(
        self, sequence: str, motif_preferences: dict[str, float]
    ) -> float:
        """
        Calculate cleavage preference factor based on motifs present in
        sequence.

        Parameters
        ----------
        sequence : str
            Local DNA sequence context
        motif_preferences : dict[str, float]
            Motif preferences where >1.0 = favored, <1.0 = disfavored

        Returns
        -------
        float
            Multiplicative factor for cleavage probability (1.0 = neutral)
        """
        if not motif_preferences or len(sequence) == 0:
            return 1.0

        preference_factor = 1.0
        sequence_upper = sequence.upper()

        for motif, preference in motif_preferences.items():
            motif_upper = motif.upper()
            if len(motif_upper) == 1:
                motif_count = sequence_upper.count(motif_upper)
                motif_frequency = motif_count / len(sequence_upper)
                preference_factor *= (
                    1.0 + (preference - 1.0) * motif_frequency
                )
            elif len(motif_upper) > 1:
                motif_count = 0
                for i in range(len(sequence_upper) - len(motif_upper) + 1):
                    if sequence_upper[i:i + len(motif_upper)] == motif_upper:
                        motif_count += 1

                if motif_count > 0:
                    # We saturate motif effects at 3 occurrences.
                    effect_strength = min(1.0, motif_count / 3.0)
                    preference_factor *= (
                        1.0 + (preference - 1.0) * effect_strength
                    )

        return max(0.1, preference_factor)

    def _get_tf_protection_factor(self, position: int, chrom: str) -> float:
        """
        Calculate transcription factor binding site protection.

        TF binding sites are protected from nuclease cleavage.
        """
        if chrom not in self.chromatin_state.tf_binding_sites:  # type: ignore
            return 1.0

        overlaps: list[object] = (
            self.chromatin_state.tf_binding_sites[  # type: ignore
                chrom
            ].at(position)
        )
        if overlaps:
            return 0.3
        else:
            return 1.0

    def _generate_fragment_sizes(
        self, num_fragments: int,
        fragment_size_params: dict[str, float] | None = None
    ) -> npt.NDArray[np.int64]:
        """
        Generate fragment sizes based on biological distributions.

        Args:
            num_fragments: Number of fragments to generate
            fragment_size_params: Optional dict of means and stds to use

        Returns:
            Array of fragment sizes

        Note:
            Users can provide the following fragment size parameters:
            - "mean": mononucleosomal mean
            - "std": mononucleosomal standard deviation
            - "mono_fraction": mononucleosomal fraction
            - "di_mean": dinucleosomal mean
            - "di_std": dinucleosomal standard deviation
            - "size_shift": size shift (applied to both peaks)
        """
        if fragment_size_params:
            mono_mean = fragment_size_params["mean"]
            mono_std = fragment_size_params["std"]
            mono_fraction = fragment_size_params.get("mono_fraction", 0.25)
            di_mean = fragment_size_params.get("di_mean", mono_mean*2)
            di_std = fragment_size_params.get("di_std", mono_std)
            size_shift = fragment_size_params.get("size_shift", 0)
        else:
            mono_mean = self.MONO_NUC_MEAN
            mono_std = self.MONO_NUC_STD
            di_mean = self.DI_NUC_MEAN
            di_std = self.DI_NUC_STD
            mono_fraction = self.MONO_FRACTION
            size_shift = 0

        num_mono: int = int(num_fragments * mono_fraction)
        num_di: int = num_fragments - num_mono

        sizes: npt.NDArray[np.float64]
        mono_sizes: npt.NDArray[np.float64] = np.random.normal(
            mono_mean+size_shift, mono_std, num_mono
        )
        di_sizes: npt.NDArray[np.float64] = np.random.normal(
            di_mean+size_shift, di_std, num_di
        )
        sizes = np.concatenate([mono_sizes, di_sizes])  # type: ignore

        sizes = self._add_periodicity(sizes, mono_mean)
        sizes_int: npt.NDArray[np.int64] = np.clip(
            sizes, self.MIN_FRAGMENT_SIZE, self.MAX_FRAGMENT_SIZE
        ).astype(np.int64)

        return sizes_int

    def _add_periodicity(
        self, sizes: npt.NDArray[np.float64], mono_mean: float
    ) -> npt.NDArray[np.float64]:
        """
        Add 10 bp periodicity centered around the mono-nucleosomal peak.

        Applies biologically realistic 10 bp periodicity primarily to fragments
        in the mono-nucleosomal range. The periodicity strength decreases with
        distance from the configured mono-nucleosomal peak.

        Args:
            sizes: Fragment sizes to modulate
            mono_mean: Mean of mono-nucleosomal peak (from configuration)

        Returns:
            Fragment sizes with added periodicity
        """
        # Calculate optimal phase to enhance the configured mono-nucleosomal
        # peak. We want mono_mean to hit a sine wave maximum (pi/2).
        optimal_phase = np.pi/2 - (2 * np.pi * mono_mean / 10) % (2 * np.pi)
        phase_factor: npt.NDArray[np.float64] = np.sin(
            2 * np.pi * sizes / 10 + optimal_phase  # type: ignore
        )

        # Weight periodicity by proximity to mono-nucleosomal range.
        distance_from_mono = np.abs(sizes - mono_mean)  # type: ignore
        periodicity_weight = np.exp(
            -(distance_from_mono ** 2) / (2 * (50 ** 2))  # type: ignore
        )

        weighted_amplitude = \
            self.PERIODICITY_AMPLITUDE * periodicity_weight  # type: ignore
        modulation: npt.NDArray[np.float64] = (
            1 + weighted_amplitude * phase_factor  # type: ignore
        )

        return sizes * modulation  # type: ignore

    def _generate_end_motif(
        self, chrom: str, position: int, motif_length: int = 4
    ) -> str:
        """
        Generate realistic end motif based on genomic sequence context and
        nuclease preferences. Optimized with sequence caching.
        """
        context = self._get_cached_sequence_context(chrom, position)
        cache_key = (chrom, position // 1000 * 1000)
        offset = position - cache_key[1]
        if offset >= motif_length and offset < len(context) - motif_length:
            base_motif = context[offset:offset + motif_length]
        else:
            base_motif = "ATCGATCGATCG"[:motif_length]

        return base_motif

    def simulate_fragments(
        self, chrom: str, start: int, end: int, num_fragments: int,
        tissue_type: str = "healthy",
        nuclease_profile: NucleaseProfile | None = None,
        fragment_size_params: dict[str, float] | None = None
    ) -> FragmentList:
        """
        Simulate cfDNA fragments from a genomic region with full biological
        realism.

        Uses realistic nucleosome positioning, chromatin accessibility, and
        nuclease-specific cleavage preferences based on experimental data.

        Args:
            chrom: Chromosome name
            start: Start position
            end: End position
            num_fragments: Number of fragments to generate
            tissue_type: Source tissue type affecting chromatin state
            nuclease_profile: Nuclease activity parameters
            fragment_size_params: Custom size distribution parameters

        Returns:
            FragmentList containing biologically realistic simulated fragments
        """

        self.logger.info(
            f"Simulating {num_fragments} fragments from {chrom}:{start}-{end}."
        )

        if nuclease_profile is None:
            nuclease_profile = NucleaseProfile()

        sizes: npt.NDArray[np.int64] = self._generate_fragment_sizes(
            num_fragments, fragment_size_params
        )
        max_size: int = int(sizes.max())  # type: ignore
        max_start_pos = max(start, end - max_size)
        if max_start_pos <= start:
            self.logger.warning(
                f"Region {chrom}:{start}-{end} too small for largest fragment "
                f"size {max_size} bp. Using minimal range."
            )
            max_start_pos = min(start + 1, end - 1)

        if max_start_pos <= start or end - start < 10:
            self.logger.error(f"Region {chrom}:{start}-{end} is too small")
            return FragmentList()

        positions: npt.NDArray[np.int64] = np.random.randint(
            start, max_start_pos, num_fragments
        )

        fragment_list: FragmentList = FragmentList()
        tissue_factor = self._get_tissue_accessibility_factor(tissue_type)
        batch_size = 1000
        for batch_start in range(0, num_fragments, batch_size):
            batch_end = min(batch_start + batch_size, num_fragments)
            batch_positions = positions[batch_start:batch_end]  # type: ignore
            batch_sizes = sizes[batch_start:batch_end]  # type: ignore

            cleave_probs = np.array([  # type: ignore
                self._get_cleavage_probability(
                    int(pos), chrom, nuclease_profile,  # type: ignore
                    tissue_factor
                ) for pos in batch_positions  # type: ignore
            ])

            random_vals = \
                np.random.random(len(batch_positions))  # type: ignore
            accepted_mask = random_vals <= cleave_probs  # type: ignore

            accepted_positions = batch_positions[accepted_mask]  # type: ignore
            accepted_sizes = batch_sizes[accepted_mask]  # type: ignore
            for pos, size in zip(  # type: ignore
                accepted_positions, accepted_sizes  # type: ignore
            ):
                end5p: str = self._generate_end_motif(
                    chrom=chrom,
                    position=int(pos)  # type: ignore
                )
                end3p: str = self._generate_end_motif(
                    chrom=chrom,
                    position=int(pos + size)  # type: ignore
                )

                fragment: Fragment = Fragment.create_simulated(
                    start_pos=int(pos),  # type: ignore
                    end_pos=int(pos + size),  # type: ignore
                    chrom=chrom,
                    length=int(size),  # type: ignore
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
