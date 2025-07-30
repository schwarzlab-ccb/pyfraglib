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

Additional references can be found in the respective classes and methods.

Example Usage
-------------
.. code-block:: python

    from pyfraglib.simulator.fragment_simulator import (
        FragmentSimulator, NucleaseProfile
    )

    # Initialize simulator with reference genome
    simulator = FragmentSimulator("/path/to/reference.fasta")

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
    experimental data such as MNase-seq or NOMe-seq experiments. The
    positioning data is used to model realistic nucleosome protection effects
    during cfDNA fragmentation.

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
    >>> open_regions = IntervalTree()
    >>> open_regions.addi(1000, 2000)  # Open chromatin region
    >>> tf_sites = IntervalTree()
    >>> tf_sites.addi(1500, 1520)  # Protected TF binding site
    >>> chromatin = ChromatinState(
    ...     open_regions=open_regions,
    ...     tf_binding_sites=tf_sites,
    ...     ctcf_sites=IntervalTree()
    ... )
    """
    #: Genomic intervals representing accessible chromatin regions (e.g.,
    #: from ATAC-seq or DNase-seq data). Fragments originating from these
    #: regions have higher cleavage probability.
    open_regions: IntervalTree  # type: ignore

    #: Transcription factor binding sites that provide protection from
    #: nuclease cleavage. These regions typically show reduced fragmentation.
    tf_binding_sites: IntervalTree  # type: ignore

    #: CTCF binding sites that create specific chromatin architecture and
    #: affect local nuclease accessibility patterns.
    ctcf_sites: IntervalTree  # type: ignore


@dataclass
class NucleaseProfile:
    """
    Nuclease activity and sequence preference parameters.

    This class defines the activity levels and sequence preferences of
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
    #: Relative activity level of DNase I. Higher values increase overall
    #: fragmentation in accessible chromatin regions. DNase I shows preference
    #: for AT-rich accessible regions.
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

    def __post_init__(self) -> None:
        if self.dnase1_motif_preference is None:
            self.dnase1_motif_preference = {
                "CG": 0.8, "GC": 0.8,  # Lower preference for CpG sites
                "AT": 1.0, "TA": 1.0,  # Neutral AT preference
                "AA": 1.1, "TT": 1.1   # Slight preference for AT-rich
            }
        if self.dnase1l3_motif_preference is None:
            self.dnase1l3_motif_preference = {
                "CC": 1.5, "GG": 1.5,  # Strong CC/GG preference
                "CG": 1.2, "GC": 1.2,  # Moderate CpG preference
                "CT": 1.1, "TC": 1.1   # Slight pyrimidine preference
            }


class SequenceContextGenerator:
    """
    Genomic sequence context provider for realistic cfDNA fragment simulation.

    This class interfaces with reference FASTA files to provide authentic
    genomic sequence context for fragment end motif generation. It eliminates
    the need for synthetic sequence generation while preserving biological
    knowledge about nuclease cleavage preferences and sequence-dependent
    fragmentation patterns.

    The class implements efficient sequence caching and provides methods for
    extracting sequence context around genomic positions, enabling realistic
    simulation of cfDNA fragments with authentic end motifs derived from actual
    genomic sequences.

    Notes
    -----
    The FASTA file must be indexed with samtools faidx before use. The class
    automatically handles edge cases such as chromosome boundaries and missing
    sequences by padding with 'N' nucleotides.

    Sequence caching is implemented to improve performance during large-scale
    simulations, with automatic cache size management to prevent memory
    overflow.

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
        self, chrom: str, position: int, length: int = 20,
        tissue_type: str = "normal"
    ) -> str:
        """
        Extract genomic sequence context around a specified position.

        This method retrieves authentic genomic sequence from the reference
        FASTA file, centered around the specified position. It handles edge
        cases at chromosome boundaries by padding with 'N' nucleotides when
        necessary.

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
        tissue_type : str, default="normal"
            Tissue type for context-specific adjustments. Currently ignored
            but kept for API compatibility.

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

        Examples
        --------
        >>> seq_gen = SequenceContextGenerator("reference.fasta")
        >>> chr1_len = seq_gen.get_chromosome_length("chr1")
        >>> print(f"Chromosome 1 length: {chr1_len:,} bp")
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

        Examples
        --------
        >>> seq_gen = SequenceContextGenerator("reference.fasta")
        >>> chroms = seq_gen.get_available_chromosomes()
        >>> print(f"Available chromosomes: {chroms[:5]}...")  # First 5
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

    def get_cleavage_motif(
        self, sequence_context: str, position: int,
        nuclease_profile: NucleaseProfile, motif_length: int = 4
    ) -> str:
        """
        Extract biologically realistic cleavage motif from sequence context.

        This method generates cleavage motifs that reflect the sequence
        preferences of different nucleases involved in cfDNA fragmentation.
        The motifs are derived from actual genomic sequence when possible,
        with nuclease-specific modifications applied based on experimental
        literature.

        Parameters
        ----------
        sequence_context : str
            Genomic sequence context from which to extract the motif.
            Should be longer than motif_length to provide sufficient context.
        position : int
            Position within the sequence_context to extract motif from.
            Must be valid within the sequence bounds.
        nuclease_profile : NucleaseProfile
            Active nuclease profile defining cleavage preferences and
            activity levels for different nucleases.
        motif_length : int, default=4
            Length of the cleavage motif to extract. Typically 4-6 nucleotides
            for cfDNA end motif analysis.

        Returns
        -------
        str
            Cleavage motif string of specified length, potentially modified
            based on nuclease preferences. Contains uppercase DNA nucleotides.

        Notes
        -----
        If the position is invalid or sequence is insufficient, a synthetic
        motif is generated based on nuclease preferences and genomic
        nucleotide frequencies.

        Examples
        --------
        >>> seq_gen = SequenceContextGenerator(\"reference.fasta\")
        >>> profile = NucleaseProfile(dnase1l3_activity=1.5)
        >>> motif = seq_gen.get_cleavage_motif(\"ATCGATCGATCG\", 4, profile, 4)
        >>> print(f\"Cleavage motif: {motif}\")
        """
        if (position < 0 or
                position + motif_length > len(sequence_context) or
                len(sequence_context) < motif_length):
            return self._generate_nuclease_preferred_motif(
                motif_length, nuclease_profile
            )

        base_motif = sequence_context[position:position + motif_length]
        return self._apply_nuclease_preferences(base_motif, nuclease_profile)

    def _generate_nuclease_preferred_motif(
        self, length: int, nuclease_profile: NucleaseProfile
    ) -> str:
        """
        Generate motif based on nuclease cleavage preferences.
        """
        motif: list[str] = []

        # @NOTE(ds): Nucleotide frequencies in human genome
        # (Lander et al. 2001 Nature).
        base_weights = {"A": 0.296, "C": 0.204, "G": 0.204, "T": 0.296}
        if (nuclease_profile.dnase1l3_activity >
                nuclease_profile.dnase1_activity):
            base_weights["C"] *= 1.5
            base_weights["G"] *= 1.3
        elif (nuclease_profile.dnase1_activity >
                nuclease_profile.dnase1l3_activity):
            base_weights["A"] *= 1.2
            base_weights["T"] *= 1.2

        if nuclease_profile.dffb_activity > 1.0:
            base_weights = {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}

        for _ in range(length):
            bases = list(base_weights.keys())
            weights = list(base_weights.values())
            total_weight = sum(weights)
            probs = [w / total_weight for w in weights]
            base = str(np.random.choice(bases, p=probs))  # type: ignore
            motif.append(base)

        return "".join(motif)

    def _apply_nuclease_preferences(
        self, motif: str, nuclease_profile: NucleaseProfile
    ) -> str:
        """
        Apply nuclease-specific sequence preferences to motif.

        Based on experimental cleavage site analysis:
        - Serpas et al. (2019) PNAS: DNASE1L3 CC preference
        - Han et al. (2020) AJHG: Comprehensive nuclease characterization
        """
        if len(motif) < 2:
            return motif

        # @NOTE(ds): DNASE1L3 strongly prefers CC dinucleotides.
        if (nuclease_profile.dnase1l3_activity >
                nuclease_profile.dnase1_activity):
            if np.random.random() < 0.3:
                if len(motif) >= 2:
                    pos = np.random.randint(0, len(motif) - 1)
                    motif_list = list(motif)
                    motif_list[pos] = "C"
                    motif_list[pos + 1] = "C"
                    return "".join(motif_list)

        if nuclease_profile.dffb_activity > 1.0:
            if np.random.random() < 0.2:
                boundary_motifs = ["CATG", "GATC", "ATAT", "GCGC"]
                if len(motif) >= 4:
                    return str(
                        np.random.choice(boundary_motifs)  # type: ignore
                    )

        return motif


class FragmentSimulator:
    """
    Biologically realistic cfDNA fragment simulator for scientific research.

    This simulator integrates multiple experimental observations to model
    cfDNA fragmentation:

    1. **Nucleosome positioning**: Nucleosome maps with occupancy scores
    2. **Chromatin accessibility**: Tissue-specific chromatin states
    3. **Nuclease specificity**: Experimentally-determined preferences
    4. **Sequence context**: Real genomic sequences from reference FASTA
    5. **TF binding protection**: Transcription factor binding site shielding

    Suitable e.g. for generating training data for machine learning models.
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

    _prob_cache: dict[tuple[str, int], float]

    def __init__(
        self, fasta_path: str,
        nucleosome_map: NucleosomeMap | None = None,
        chromatin_state: ChromatinState | None = None,
        sequence_generator: SequenceContextGenerator | None = None
    ) -> None:
        """
        Initialize the fragment simulator.

        Args:
            fasta_path: Path to reference FASTA file (must be indexed)
            nucleosome_map: Custom nucleosome positioning data
            chromatin_state: Chromatin accessibility information
            sequence_generator: Custom sequence context generator
        """
        self.logger: logging.Logger = logging.getLogger(LOGGER_NAME)
        self.fasta_path: str = fasta_path
        self.sequence_generator: SequenceContextGenerator = (
            sequence_generator or SequenceContextGenerator(fasta_path)
        )
        self.nucleosome_map: NucleosomeMap = (
            nucleosome_map or self._default_nucleosome_model()
        )
        self.chromatin_state: ChromatinState = (
            chromatin_state or self._default_chromatin_state()
        )

        self._prob_cache = {}
        self._sequence_cache: dict[tuple[str, int], str] = {}
        self._nucleosome_cache: dict[tuple[str, int], tuple[float, float]] = {}

    def _default_nucleosome_model(self) -> NucleosomeMap:
        """
        Generate a default nucleosome map based on chromosome data from FASTA.
        Uses a simplified model with approximately 147 bp nucleosomes and
        variable linkers.
        """
        self.logger.info("Generating default nucleosome model.")

        positions: dict[str, npt.NDArray[np.int64]] = {}
        occupancy: dict[str, npt.NDArray[np.float64]] = {}

        available_chroms = self.sequence_generator.get_available_chromosomes()
        for chrom in available_chroms:
            chrom_len: int = (
                self.sequence_generator.get_chromosome_length(chrom)
            )
            if chrom_len < 10_000:  # Skip very short contigs
                continue

            base_positions: npt.NDArray[np.int64] = np.arange(
                100, chrom_len - 100, 185  # 147 bp nucleosome + ~38 bp linker
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
        self, position: int, chrom: str, nuclease_profile: NucleaseProfile,
        tissue_type: str = "healthy", tissue_factor: float | None = None
    ) -> float:
        """
        Calculate biologically realistic cleavage probability at a genomic
        position.

        Integrates multiple biological factors:
        - Nucleosome positioning and occupancy
        - Chromatin accessibility state
        - Tissue-specific chromatin properties
        - Nuclease-specific preferences (DNASE1, DNASE1L3, DFFB)

        Args:
            position: Genomic position (0-based)
            chrom: Chromosome name
            nuclease_profile: Active nucleases and their preferences
            tissue_type: Tissue type affecting chromatin accessibility

        Returns:
            Cleavage probability between 0.0 and 1.0
        """
        # Start with base accessibility from chromatin state
        base_prob = 0.1  # Background cleavage probability

        # Check if position is in accessible chromatin regions
        if chrom in self.chromatin_state.open_regions:  # type: ignore
            overlaps = (
                self.chromatin_state.open_regions[  # type: ignore
                    chrom
                ][position]
            )
            if overlaps:  # type: ignore
                base_prob = 0.6  # Higher probability in open chromatin

        nucleosome_factor = self._get_nucleosome_protection(position, chrom)
        if tissue_factor is None:
            tissue_factor = self._get_tissue_accessibility_factor(tissue_type)

        nuclease_factor = self._get_nuclease_activity_factor(
            position, chrom, nuclease_profile
        )
        total_prob = (base_prob * nucleosome_factor * tissue_factor *
                      nuclease_factor)

        tf_protection = self._get_tf_protection_factor(position, chrom)
        total_prob *= tf_protection
        return min(max(total_prob, 0.001), 1.0)

    def _get_nucleosome_protection(self, position: int, chrom: str) -> float:
        """
        Calculate nucleosome protection factor at a genomic position.

        Optimized with caching and efficient spatial indexing.
        """
        if chrom not in self.nucleosome_map.positions:
            return 0.3  # Default for unmapped chromosomes

        # @NOTE(ds): Cache nucleosome calculations in 100 bp windows.
        cache_key = (chrom, position // 100)
        if cache_key in self._nucleosome_cache:
            cached_distance, cached_occupancy = (
                self._nucleosome_cache[cache_key])
        else:
            nucleosome_positions = self.nucleosome_map.positions[chrom]
            nucleosome_occupancy = self.nucleosome_map.occupancy[chrom]

            nearest_idx = np.searchsorted(
                nucleosome_positions, cache_key[1] * 100
            )
            candidates = []
            if nearest_idx > 0:
                candidates.append(nearest_idx - 1)
            if nearest_idx < len(nucleosome_positions):
                candidates.append(nearest_idx)

            if not candidates:
                cached_distance, cached_occupancy = 200, 0.5
            else:
                center_pos = cache_key[1] * 100
                distances = np.abs(  # type: ignore
                    nucleosome_positions[candidates]  # type: ignore
                    - center_pos
                )
                best_candidate_idx = np.argmin(distances)  # type: ignore
                best_idx = candidates[best_candidate_idx]
                cached_distance = distances[best_candidate_idx]  # type: ignore
                cached_occupancy = nucleosome_occupancy[best_idx]

            self._nucleosome_cache[cache_key] = \
                cached_distance, cached_occupancy

            if len(self._nucleosome_cache) > 50000:
                old_keys = list(self._nucleosome_cache.keys())[:12500]
                for old_key in old_keys:
                    del self._nucleosome_cache[old_key]

        actual_distance = abs(cached_distance +
                              (position - cache_key[1] * 100))
        if actual_distance <= 73:
            protection = 0.05 + (1.0 - cached_occupancy) * 0.2
        elif actual_distance <= 120:
            protection = 0.3 + (1.0 - cached_occupancy) * 0.4
        else:
            protection = 0.8 + (1.0 - cached_occupancy) * 0.2

        return protection

    def _get_tissue_accessibility_factor(self, tissue_type: str) -> float:
        """
        Get tissue-specific chromatin accessibility factor.

        Based on tissue-specific chromatin openness patterns from
        ATAC-seq and DNase-seq data (ENCODE Consortium).
        """
        tissue_factors = {
            "healthy": 1.0,
            "hematopoietic": 1.2,  # More open chromatin
            "liver": 0.9,          # Slightly more compact
            "placenta": 1.1,       # Moderately open
            "tumor": 1.4           # Highly accessible due to disrupted
        }
        return tissue_factors.get(tissue_type, 1.0)

    def _get_nuclease_activity_factor(
        self, position: int, chrom: str, nuclease_profile: NucleaseProfile
    ) -> float:
        """
        Calculate nuclease-specific activity factor based on local sequence
        context.
        """
        cache_key = (chrom, position // 1000 * 1000)

        if cache_key in self._sequence_cache:
            context = self._sequence_cache[cache_key]
        else:
            try:
                context = self.sequence_generator.fasta.fetch(
                    chrom, cache_key[1], cache_key[1] + 1200
                )
                self._sequence_cache[cache_key] = context
                if len(self._sequence_cache) > 5000:
                    oldest_keys = list(self._sequence_cache.keys())[:2000]
                    for old_key in oldest_keys:
                        del self._sequence_cache[old_key]
            except Exception:
                context = "N" * 1200

        offset = position - cache_key[1]
        if offset >= 10 and offset < len(context) - 30:
            local_context = context[offset-10:offset+10]
        else:
            if offset < 10:
                local_context = (context[:20] if len(context) >= 20
                                 else context)
            else:
                local_context = (context[-20:] if len(context) >= 20
                                 else context)

        if len(local_context) > 0:
            gc_content = ((local_context.count('G') +
                           local_context.count('C')) / len(local_context))
            cc_count = local_context.count('CC')
            at_content = ((local_context.count('A') +
                           local_context.count('T')) / len(local_context))
        else:
            gc_content = 0.42  # Human genome average
            cc_count = 1
            at_content = 0.58

        activity = 0.0
        if nuclease_profile.dnase1l3_activity > 0:
            dnase1l3_factor = 1.0 + (cc_count * 0.3) + (gc_content * 0.2)
            activity += nuclease_profile.dnase1l3_activity * dnase1l3_factor

        if nuclease_profile.dnase1_activity > 0:
            dnase1_factor = 1.0 + (at_content * 0.2)
            activity += nuclease_profile.dnase1_activity * dnase1_factor

        if nuclease_profile.dffb_activity > 0:
            nucleosome_factor = self._get_nucleosome_protection(position,
                                                                chrom)
            dffb_factor = 2.0 - nucleosome_factor
            activity += nuclease_profile.dffb_activity * dffb_factor

        total_nuclease = (nuclease_profile.dnase1_activity +
                          nuclease_profile.dnase1l3_activity +
                          nuclease_profile.dffb_activity)

        if total_nuclease > 0:
            return activity / total_nuclease
        else:
            return 1.0

    def _get_tf_protection_factor(self, position: int, chrom: str) -> float:
        """
        Calculate transcription factor binding site protection.

        TF binding sites are protected from nuclease cleavage.
        """
        if chrom not in self.chromatin_state.tf_binding_sites:  # type: ignore
            return 1.0

        overlaps = (
            self.chromatin_state.tf_binding_sites[  # type: ignore
                chrom
            ][position]
        )
        if overlaps:  # type: ignore
            return 0.3
        else:
            return 1.0

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
            size_distribution: Type of dist ("normal", "cancer")
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
        """Add 10 bp periodicity to fragment sizes."""
        phase_factor: npt.NDArray[np.float64] = np.sin(
                2 * np.pi * sizes / 10 + self.PERIODICITY_PHASE  # type: ignore
        )
        modulation: npt.NDArray[np.float64] = (
            1 + self.PERIODICITY_AMPLITUDE * phase_factor  # type: ignore
        )
        return sizes * modulation  # type: ignore

    def _generate_end_motif(
        self, chrom: str, position: int, nuclease_profile: NucleaseProfile,
        tissue_type: str = "normal", motif_length: int = 4
    ) -> str:
        """
        Generate realistic end motif based on genomic sequence context and
        nuclease preferences. Optimized with sequence caching.
        """
        cache_key = (chrom, position // 1000 * 1000)

        if cache_key in self._sequence_cache:
            context = self._sequence_cache[cache_key]
        else:
            try:
                context = self.sequence_generator.fasta.fetch(
                    chrom, cache_key[1], cache_key[1] + 1200
                )
                self._sequence_cache[cache_key] = context
                if len(self._sequence_cache) > 5000:
                    oldest_keys = list(self._sequence_cache.keys())[:2000]
                    for old_key in oldest_keys:
                        del self._sequence_cache[old_key]
            except Exception:
                context = "N" * 1200

        offset = position - cache_key[1]
        if offset >= motif_length and offset < len(context) - motif_length:
            base_motif = context[offset:offset + motif_length]
        else:
            base_motif = "ATCG"[:motif_length]

        return self.sequence_generator._apply_nuclease_preferences(
            base_motif, nuclease_profile
        )

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
        tissue_to_size_dist: dict[str, str] = {
            "healthy": "normal",
            "hematopoietic": "normal",
            "liver": "normal",
            "placenta": "normal",
            "tumor": "cancer"
        }
        size_dist: str = tissue_to_size_dist.get(tissue_type, "normal")
        sizes: npt.NDArray[np.int64] = self._generate_fragment_sizes(
            num_fragments, size_dist, fragment_size_params
        )

        max_size: int = int(sizes.max())  # type: ignore
        effective_end = max(start + max_size + 1, end - max_size)
        if effective_end <= start:
            self.logger.warning(
                f"Region {chrom}:{start}-{end} too small for fragment size "
                f"{max_size}"
            )
            effective_end = start + max_size + 1

        positions: npt.NDArray[np.int64] = np.random.randint(
            start, effective_end, num_fragments
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
                    tissue_type, tissue_factor
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
                    position=int(pos),  # type: ignore
                    nuclease_profile=nuclease_profile,
                    tissue_type=tissue_type
                )
                end3p: str = self._generate_end_motif(
                    chrom=chrom,
                    position=int(pos + size),  # type: ignore
                    nuclease_profile=nuclease_profile,
                    tissue_type=tissue_type
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
