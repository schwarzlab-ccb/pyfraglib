"""
Fragmentomics Scores
====================

This module provides algorithms for calculating fragmentomics scores from cfDNA
fragment data and visualizing them along (selected parts of) the genome.

Key Algorithms
--------------
- **Windowed Protection Score (WPS)**: Measures nucleosome positioning and
  chromatin accessibility
- **Motif Diversity**: Quantifies diversity of fragment end motifs using
  entropy measures

Windowed Protection Score (WPS)
-------------------------------
The WPS algorithm measures the degree to which genomic regions are protected
from nuclease digestion, providing insights into nucleosome positioning and
chromatin accessibility. Please refer to the API documentation for a detailed
explanation of how the WPS is calculated.

Motif Diversity
---------------
Fragment end motifs reflect nuclease cleavage preferences and can reveal
tissue-specific patterns. We provide functions to calculate diversity indices
for fragment end motifs (i.e. Shannon entropy and Simpson index).

Example Usage
-------------
Calculate WPS for genomic regions:

.. code-block:: python

    from pyfraglib.scores import windowed_protection_score
    from pyfraglib.fragment import Fragment

    # Load fragments
    fragments = Fragment.from_bam("sample.bam", "variants.vcf")

    # Calculate WPS for BED regions
    wps_scores = windowed_protection_score(fragments, "regions.bed")

    # Access results
    for region, score in wps_scores.items():
        print(f"Region {region}: WPS = {score:.3f}")

Calculate motif diversity:

.. code-block:: python

    from pyfraglib.scores import motif_diversity

    # Calculate Shannon entropy for 4-mer end motifs
    shannon_div = motif_diversity(fragments, kmer_len=4, index="shannon")
    print(f"Shannon entropy: {shannon_div:.4f}")

    # Calculate Simpson index
    simpson_div = motif_diversity(fragments, kmer_len=4, index="simpson")
    print(f"Simpson index: {simpson_div:.4f}")

Visualize scores:

.. code-block:: python

    from pyfraglib.scores import score_line_plot

    # Create line plot of WPS scores
    score_line_plot(wps_scores, "output/", "sample_wps",
                    xlabel="Genomic Position", ylabel="WPS Score",
                    title="Windowed Protection Score")

Performance Considerations
--------------------------
- Consider subsampling fragments for initial exploratory analysis
- WPS calculation is memory-intensive for large genomic regions

Functions
---------
- :func:`windowed_protection_score`: WPS implementation
- :func:`motif_diversity`: Calculate diversity indices for fragment end motifs
- :func:`score_line_plot`: Generate line plots for fragmentomics scores

License
-------
This file is part of ``pyfraglib``, a software suite to calculate fragmentomics
features from cfDNA and perform downstream analyses.

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
import typing
import os
import pysam
import matplotlib

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt

from collections import defaultdict
from pyfraglib.core import shannon_entropy, simpson_index, fail, \
                           get_chromosome_length, hg19_chromosomes, \
                           hg38_chromosomes, homogenize_contig_name
from pyfraglib.fragment import FragmentList
from pyfraglib import get_logger


def motif_diversity(
    fragments: FragmentList, name: str, index: str = "shannon"
) -> tuple[float, float]:
    """
    Calculate diversity indices for fragment end motifs.

    Args:
        fragments: List of fragments to analyze. Bogus fragments are
            automatically excluded from calculations.
        name: Sample or analysis name for logging purposes. Used in log
            messages to identify which sample is being processed.
        index: Diversity index to calculate. Options are:

            * ``"shannon"``: Shannon entropy (default)
            * ``"simpson"``: Simpson index

    Returns:
        tuple[float, float]: A tuple containing (5p_diversity, 3p_diversity).
        Returns (0.0, 0.0) if no valid fragments are found

    Raises:
        SystemExit: If an invalid diversity index is specified. Only "shannon"
            and "simpson" are supported.

    Note:
        * Uses 4-mer sequences from fragment ends for analysis
        * Automatically filters out motifs containing 'N' bases
        * Zero counts are excluded to avoid log(0) issues in calculations
        * Results are normalized and comparable across samples

    See Also:
        * :func:`pyfraglib.core.shannon_entropy` - Shannon entropy calculation
        * :func:`pyfraglib.core.simpson_index` - Simpson index calculation
        * :class:`FragmentList.count_endmotifs` - Low-level motif counting
    """
    logger: logging.Logger = get_logger()

    index_func: typing.Callable[[list[float]], float]
    if index == "shannon":
        index_func = shannon_entropy
    elif index == "simpson":
        index_func = simpson_index
    else:
        fail(f"motif_diversity(): unknown index function ``{index}``")

    logger.info(f"calculating motif diversity score ({index}) for {name}")

    ends_5p: defaultdict[str, int] = defaultdict(int)
    ends_3p: defaultdict[str, int] = defaultdict(int)

    ends_5p, ends_3p, num_frags = fragments.count_endmotifs(kmer_len=4)

    # @NOTE(ds): Filter out zero counts to avoid log(0) issues in diversity
    # score calculations.
    valid_5p_counts: list[int] = [val for val in ends_5p.values() if val > 0]
    valid_3p_counts: list[int] = [val for val in ends_3p.values() if val > 0]

    total_5p: int = sum(valid_5p_counts)
    proportions_5p: list[float] = \
        [val / total_5p for val in valid_5p_counts] if total_5p > 0 else []
    index_5p: float = index_func(proportions_5p) if proportions_5p else 0.0

    total_3p: int = sum(valid_3p_counts)
    proportions_3p: list[float] = \
        [val / total_3p for val in valid_3p_counts] if total_3p > 0 else []
    index_3p: float = index_func(proportions_3p) if proportions_3p else 0.0

    return (index_5p, index_3p)


def windowed_protection_score(
    fragments: FragmentList, regions: pysam.TabixFile, win_size: int = 120,
    genome: str = "hg19"
) -> pd.DataFrame:
    """
    Calculate Windowed Protection Score (WPS) for genomic regions.

    Computes the Windowed Protection Score, which measures the degree to which
    genomic regions are protected from nuclease digestion. The WPS algorithm
    calculates protection as the difference between spanning and ending
    fragments within genomic windows. Higher WPS values indicate regions with
    regular nucleosome positioning, while lower values suggest accessible
    chromatin or irregular positioning.

    Args:
        fragments: Collection of fragments to analyze. Bogus fragments are
            automatically excluded from WPS calculations.
        regions: Indexed BED file containing genomic regions of interest.
            Must be opened with pysam.TabixFile and properly formatted.
        win_size: Window size in base pairs for WPS calculation. Default is
            120bp, which corresponds to nucleosome-sized protection. Must be
            positive.
        genome: Reference genome version for chromosome length validation.
            Supported values are "hg19" and "hg38".

    Returns:
        pd.DataFrame: DataFrame containing WPS results with columns:

            * ``chrom``: Chromosome name
            * ``pos``: Genomic position
            * ``abs_pos``: Absolute position across the genome
            * ``wps``: Windowed Protection Score
            * ``depth``: Fragment coverage depth
            * ``spanning_frags``: Number of fragments spanning the window
            * ``ending_frags``: Number of fragments ending in the window
            * ``info``: Additional information from BED file

    Note:
        * Results include both raw WPS scores and supporting coverage metrics

    See Also:
        * :func:`windowed_protection_score_fast` - Optimized WPS implementation
        * :func:`score_line_plot` - Visualization of WPS results
    """
    return windowed_protection_score_fast(fragments, regions, win_size, genome)


def windowed_protection_score_fast(
    fragments: FragmentList, regions: pysam.TabixFile, win_size: int = 120,
    genome: str = "hg19"
) -> pd.DataFrame:
    """
    This fast WPS implementation replaces the previous ``pyfraglib`` WPS
    implementation. We removed the ``step_size`` argument because this
    algorithm does not need that. It iterates over fragments instead of genomic
    regions like so:

    > case 1 (fragment_size < window_size):
          [frag_start - win_size/2, frag_end + win_size/2] <- -1

    > case 2 (fragment_size >= window_size):
          [frag_start - win_size/2, frag_start + win_size/2) <- -1

          [frag_start + win_size/2, frag_end   - win_size/2] <- +1

          (frag_end   - win_size/2, frag_end   + win_size/2] <- -1

    The implementation below tries to be extremely careful about interval
    definitions to not introduce 1-off errors. The arguments are identical to
    :func:`windowed_protection_score`. See there for more information.
    """
    assert win_size > 0

    chrom_map_wps: dict[str, npt.NDArray[np.int64]] = \
        create_chromosome_map(genome)
    chrom_map_depth: dict[str, npt.NDArray[np.int64]] = \
        create_chromosome_map(genome)

    # @NOTE(ds): We also want to record the exact number of spanning and ending
    # fragments per window. As with the WPS, we record those numbers at the
    # window midpoint.
    chrom_map_end: dict[str, npt.NDArray[np.int64]] = \
        create_chromosome_map(genome)
    chrom_map_span: dict[str, npt.NDArray[np.int64]] = \
        create_chromosome_map(genome)

    for fragment in fragments:
        if fragment.is_bogus:
            continue

        frag_start: int = fragment.start_pos
        frag_end: int = fragment.end_pos  # 1 past end
        frag_len: int = fragment.length
        frag_chrom: str = homogenize_contig_name(fragment.chrom)

        win_half: int = win_size // 2
        istart: int = max(frag_start - win_half, 0)
        iend: int = min(frag_end + win_half,
                        get_chromosome_length(frag_chrom, genome))

        chrom_map_depth[frag_chrom][frag_start:frag_end] += 1
        if frag_len < win_size:
            chrom_map_wps[frag_chrom][istart:iend] -= 1
            chrom_map_end[frag_chrom][istart:iend] += 1
        else:
            ostart: int = frag_start + win_half
            oend: int = frag_end - win_half  # 1 past end

            chrom_map_wps[frag_chrom][istart:ostart] -= 1
            chrom_map_wps[frag_chrom][ostart:oend] += 1
            chrom_map_wps[frag_chrom][oend:iend] -= 1

            chrom_map_end[frag_chrom][istart:ostart] += 1
            chrom_map_span[frag_chrom][ostart:oend] += 1
            chrom_map_end[frag_chrom][oend:iend] += 1

    return chromosome_maps_to_df(chrom_map_wps, chrom_map_depth,
                                 chrom_map_span, chrom_map_end,
                                 regions, genome)


def precalc_size(regions: pysam.TabixFile) -> int:
    """
    Internal utility function.
    """
    it: int = 0
    for region in regions.fetch():
        istart: str
        iend: str
        _, istart, iend, _ = region.split()
        it += int(iend) - int(istart) + 1
    regions.reset()
    return it


def create_chromosome_map(
    genome: str = "hg19"
) -> dict[str, npt.NDArray[np.int64]]:
    """
    Internal utility function.
    """
    chrom_map: dict[str, npt.NDArray[np.int64]] = dict()
    name: str
    length: int
    chromosomes: list[tuple[str, int, str, str]]

    if genome == "hg19":
        chromosomes = hg19_chromosomes
    elif genome == "hg38":
        chromosomes = hg38_chromosomes
    else:
        fail(f"unknown genome ``{genome}`` requested")

    for name, length, _, _ in chromosomes:
        chrom_map[name] = np.zeros(length, dtype=np.int64)

    return chrom_map


def chromosome_maps_to_df(
    chrom_map_wps: dict[str, npt.NDArray[np.int64]],
    chrom_map_depth: dict[str, npt.NDArray[np.int64]],
    chrom_map_span: dict[str, npt.NDArray[np.int64]],
    chrom_map_end: dict[str, npt.NDArray[np.int64]],
    regions: pysam.TabixFile, genome: str = "hg19"
) -> pd.DataFrame:
    """
    Internal utility function.

    Note:
        The input BED file must be sorted with respect to chromosomes.
        E.g. listing chr1 regions, chr3 regions, and chr2 regions in this order
        will produce incorrect results with respect to the absolute genome
        position!
    """
    col_names: list[str] = ["chrom", "pos", "abs_pos", "wps", "depth",
                            "spanning_frags", "ending_frags", "info"]
    output_df: pd.DataFrame = pd.DataFrame(
        None, index=range(precalc_size(regions)),
        columns=col_names
    )

    region: str
    it: int = 0
    cur_chrom: str | None = None
    cum_pos: int = 0
    chromosome_wps: npt.NDArray[np.int64] | None = None
    chromosome_depth: npt.NDArray[np.int64] | None = None
    chromosome_span: npt.NDArray[np.int64] | None = None
    chromosome_end: npt.NDArray[np.int64] | None = None

    for region in regions.fetch():
        chrom: str
        istart: str
        iend: str
        info: str

        chrom, istart, iend, info = region.split()
        chrom = homogenize_contig_name(chrom)

        if not cur_chrom or \
                chromosome_wps is None or \
                chromosome_depth is None or \
                chromosome_span is None or \
                chromosome_end is None:
            cur_chrom = chrom
            chromosome_wps = chrom_map_wps[chrom]
            chromosome_depth = chrom_map_depth[chrom]
            chromosome_span = chrom_map_span[chrom]
            chromosome_end = chrom_map_end[chrom]
        elif chrom != cur_chrom:
            cum_pos += get_chromosome_length(cur_chrom, genome)
            cur_chrom = chrom
            chromosome_wps = chrom_map_wps[chrom]
            chromosome_depth = chrom_map_depth[chrom]
            chromosome_span = chrom_map_span[chrom]
            chromosome_end = chrom_map_end[chrom]

        rlen: int = int(iend) - int(istart)  # region length - 1
        rel_pos: npt.NDArray[np.int64] = np.arange(int(istart), int(iend) + 1)

        output_df.loc[it:(it+rlen), "chrom"] = chrom
        output_df.loc[it:(it+rlen), "pos"] = rel_pos
        output_df.loc[it:(it+rlen), "abs_pos"] = rel_pos + cum_pos
        output_df.loc[it:(it+rlen), "wps"] = chromosome_wps[rel_pos]
        output_df.loc[it:(it+rlen), "depth"] = chromosome_depth[rel_pos]
        output_df.loc[it:(it+rlen), "spanning_frags"] = \
            chromosome_span[rel_pos]
        output_df.loc[it:(it+rlen), "ending_frags"] = chromosome_end[rel_pos]
        output_df.loc[it:(it+rlen), "info"] = info

        it += rlen + 1

    regions.reset()

    return output_df


def score_line_plot(
    df: pd.DataFrame, name: str, out_dir: str, score: str = "wps",
    exclude_chroms: list[str] = ["Y", "M"], region: tuple[int, int] = (0, 0),
    genome: str = "hg19", log_transform: bool = False,
    plot_spanning_ending_frags: bool = False
) -> None:
    """
    Generate line plots for fragmentomics scores across the genome.

    Creates exploratory line plots showing fragmentomics scores across
    chromosomes or specific genomic regions. Supports multiple score types
    including WPS, fragment depth, and fragment ratios with customizable
    visualization options.

    The function generates genome-wide plots or focused plots for single
    chromosome regions.

    Args:
        df: DataFrame containing score data with required columns:
            ``chrom``, ``pos``, ``abs_pos``, and the specified score column.
            Additional columns ``spanning_frags`` and ``ending_frags`` required
            if ``plot_spanning_ending_frags=True``.
        name: Sample or analysis name used for the plot title and output
            filename.
        out_dir: Output directory path where the plot PNG file will be saved.
        score: Score type to plot. Supported options:

            * ``"wps"``: Windowed Protection Score (default)
            * ``"depth"``: Fragment coverage depth
            * ``"ratio_end_span"``: Ratio of ending to spanning fragments
            * ``"ratio_span_total"``: Ratio of spanning to total fragments
        exclude_chroms: List of chromosomes to exclude from plotting.
            Default excludes Y and mitochondrial chromosomes.
        region: Genomic region to focus on as (start, end) coordinates.
            Only applied when plotting a single chromosome. Default (0, 0)
            plots the entire chromosome.
        genome: Reference genome version for chromosome definitions.
            Supported values are "hg19" and "hg38".
        log_transform: Whether to apply logarithmic transformation to y-axis.
            Useful for scores with wide dynamic ranges.
        plot_spanning_ending_frags: Whether to overlay spanning and ending
            fragment counts on the plot. Only available for single-chromosome
            plots.

    Returns:
        None. The function saves a PNG plot to the specified output directory
        with filename ``{name}_{score}.png``.

    Raises:
        SystemExit: If an unsupported score type is specified.

    Note:
        * Automatically calculates ratio scores if not present in DataFrame

    Example:
        Create WPS line plot across the genome::

            import pandas as pd
            from pyfraglib.scores import score_line_plot

            # Assuming wps_results is a df from windowed_protection_score
            score_line_plot(
                df=wps_results,
                name="sample_001",
                out_dir="plots/",
                score="wps",
                exclude_chroms=["Y", "M"]
            )
            # Creates: plots/sample_001_wps.png

    See Also:
        * :func:`windowed_protection_score` - Generate WPS data for plotting
    """
    if score not in ["wps", "depth", "ratio_end_span", "ratio_span_total"]:
        fail(f"unknown score ``{score}`` requested")

    outpath: str = os.path.join(out_dir, f"{name}_{score}.png")
    get_logger().info(
        f"saving score plot for {name} to ``{outpath}``"
    )

    chromosomes: list[tuple[str, int, str, str]]
    if genome == "hg19":
        chromosomes = hg19_chromosomes
    elif genome == "hg38":
        chromosomes = hg38_chromosomes
    else:
        fail(f"unknown genome ``{genome}`` requested")

    if score == "ratio_end_span":
        # Avoid division by zero: set ratio to 0 when spanning_frags is 0.
        denominator = df["spanning_frags"].replace(0, 1)
        df["ratio_end_span"] = df["ending_frags"] / denominator
        df.loc[df["spanning_frags"] == 0, "ratio_end_span"] = (
            0.0
        )

    if score == "ratio_span_total":
        # Avoid division by zero: set ratio to 0 when total is 0.
        total_frags = df["ending_frags"] + df["spanning_frags"]
        denominator = total_frags.replace(0, 1)
        df["ratio_span_total"] = \
            df["spanning_frags"] / denominator
        df.loc[total_frags == 0, "ratio_span_total"] = 0.0

    fig: matplotlib.figure.Figure = plt.figure()
    plt.title(f"Sample {name}")

    cumsum: int = 0
    chrom_midpoints: list[int] = []
    chrom_names: list[str] = []
    color: str = "lightgrey"
    for name, length, _, _ in chromosomes:
        if name in exclude_chroms:
            continue

        chrom_midpoints.append(cumsum + length//2)
        chrom_names.append(name)

        plt.axvspan(cumsum, cumsum+length, color=color, alpha=0.3)
        if color == "lightgrey":
            color = "white"
        else:
            color = "lightgrey"

        condition = df["chrom"] == name
        plt.plot(
            df.loc[condition, "pos"] + cumsum,
            df.loc[condition, score],
            color="#1f77b4",
            label=score,
            linewidth=2.5
        )
        if plot_spanning_ending_frags:
            plt.plot(
                df.loc[condition, "pos"] + cumsum,
                df.loc[condition, "spanning_frags"],
                color="#beddf1",
                label="spanning fragments",
            )
            plt.plot(
                df.loc[condition, "pos"] + cumsum,
                df.loc[condition, "ending_frags"],
                color="#f8c57c",
                label="ending fragments",
            )

        cumsum += length

    if len(chrom_names) == 1:
        plt.xlim(region)
        plt.xticks([region[0], region[1], np.mean(region)],  # type: ignore
                   [region[0], region[1], chrom_names[0]])  # type: ignore
    else:
        plt.xlim(-cumsum*0.05, cumsum*1.05)
        staggered_chrom_names: list[str] = [
            f"{label}\n" if i % 2 == 0 else f"\n{label}"
            for i, label in enumerate(chrom_names)
        ]
        plt.xticks(chrom_midpoints, staggered_chrom_names, rotation=0)

    plt.xlabel("Chromosomes")
    if log_transform:
        plt.yscale("log")
        plt.ylabel("Logarithm of value")
    else:
        plt.ylabel(score.capitalize())

    if len(chrom_names) == 1:
        plt.legend()

    plt.tight_layout()
    fig.savefig(outpath, dpi=900)
