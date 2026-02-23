"""
Fragment Statistics and Visualization
=====================================

This module provides overview statistics and visualization for cfDNA fragment
lists.

Key Features
------------
- **Fragment Distribution Analysis**: Chromosome-wise fragment counts with
  mutation status visualization
- **End Motif Analysis**: k-mer frequency analysis for 5' and 3' fragment ends
- **Quality Metrics**: Statistical summaries including bogus and mutated
  fragment counts
- **CSV Export Functions**: Export fragment length distributions and end motif
  counts to CSV format for external analysis

License
-------
This file is part of ``pyfraglib``, a software suite to calculate
fragmentomics features from cfDNA and perform downstream analyses.

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
import os
import csv
import json
import logging
import itertools

import matplotlib
import matplotlib.pyplot as plt

from collections import defaultdict, Counter
from pyfraglib.core import fail
from pyfraglib.fragment import FragmentList
from pyfraglib import get_logger


def fragments_per_chromosome_barplot(
    fragments: FragmentList, out_dir: str, name: str
) -> None:
    """
    Generate chromosome-wise fragment distribution plots with mutation status.

    Creates a stacked bar chart showing fragment counts per chromosome, with
    separate visualization of mutated versus wildtype fragments. This provides
    an overview of fragment distribution across chromosomes and helps identify
    chromosomal bias in mutation patterns.

    Parameters
    ----------
    fragments : FragmentList
        Collection of fragments to analyze. Bogus fragments are automatically
        excluded from the analysis.
    out_dir : str
        Output directory path where the PNG file will be saved. Directory
        must exist and be writable.
    name : str
        Sample identifier used for plot title and output filename. Used to
        generate informative plot titles and unique output filenames.

    Returns
    -------
    A PNG file named ``{name}_mut_frags_per_chrom.png``.
    """
    logger: logging.Logger = get_logger()
    logger.info(
        "creating fragment-per-chromosome histograms in ``{}``".format(out_dir)
    )

    chrom_map_mut: defaultdict[str, int] = defaultdict(int)
    chrom_map_wt:  defaultdict[str, int] = defaultdict(int)
    num_frags = 0
    for frag in fragments:
        if frag.is_bogus:
            continue

        num_frags += 1
        if frag.is_mutated:
            chrom_map_mut[frag.chrom] += 1
            _ = chrom_map_wt[frag.chrom]
        else:
            chrom_map_wt[frag.chrom] += 1
            _ = chrom_map_mut[frag.chrom]

    fig = plt.figure()
    outpath = os.path.join(out_dir, "{}_mut_frags_per_chrom.png".format(name))

    assert chrom_map_wt.keys() == chrom_map_mut.keys()
    plt.bar(list(chrom_map_wt.keys()), list(chrom_map_wt.values()),
            label="wildtype fragments")
    if sum(chrom_map_mut.values()) != 0:
        plt.bar(list(chrom_map_mut.keys()), list(chrom_map_mut.values()),
                bottom=list(chrom_map_wt.values()), label="mutated fragments")

    plt.xticks(rotation=45)  # better rotate bc. we do not homogenize chr names
    plt.title("{}: Fragments per chromosome n={}".format(name, num_frags))
    plt.xlabel("Chromosomes")
    plt.ylabel("# of Fragments")
    plt.legend()

    fig.savefig(outpath, dpi=900)


def end_motifs_barplot(
    fragments: FragmentList, out_dir: str, name: str, kmer_len: int
) -> None:
    """
    Generate color-coded bar plots for fragment end motif frequencies.

    Creates separate bar plots for 5' and 3' end motif distributions, with
    base-specific color coding to visualize nuclease cleavage preferences.
    These plots reveal tissue-specific fragmentation patterns and cleavage
    site preferences.

    Parameters
    ----------
    fragments : FragmentList
        Fragment collection for motif analysis. Bogus fragments are
        automatically excluded.
    out_dir : str
        Output directory for PNG files. Must exist and be writable.
    name : str
        Sample identifier for plot titles and filenames.
    kmer_len : int
        Length of k-mer motifs to analyze. Typically 3 or 4 nucleotides.
        Determines the resolution of motif analysis.

    Returns
    -------
    ``{name}_k{kmer_len}_5p_frag_end_motifs.png`` - 5' end motifs
    ``{name}_k{kmer_len}_3p_frag_end_motifs.png`` - 3' end motifs
    """
    logger: logging.Logger = get_logger()
    logger.info(
        "creating end motif barplots in ``{}``".format(out_dir)
    )

    motifs_5p: defaultdict[str, int]
    motifs_3p: defaultdict[str, int]
    num_frags: int

    motifs_5p, motifs_3p, num_frags = fragments.count_endmotifs(kmer_len)
    assert motifs_5p.keys() == motifs_3p.keys()

    for label, motifs in [("5", motifs_5p), ("3", motifs_3p)]:
        mzip: zip[tuple[str, int]] = zip(motifs.keys(), motifs.values())
        szip: list[tuple[str, int]] = sorted(mzip, reverse=False)

        motifs_names: list[str] = []
        motifs_vals: list[int] = []
        m_name: str
        m_val: int
        for m_name, m_val in szip:
            motifs_names.append(m_name)
            motifs_vals.append(m_val)

        # @NOTE(ds): We assume sections of equal lengths!
        max_count: int = max(motifs_vals)
        sec_len: int = len(motifs_names) // 4

        fig: matplotlib.figure.Figure = plt.figure()
        outpath: str = os.path.join(
            out_dir, "{}_k{}_{}p_frag_end_motifs.png".format(
                name, kmer_len, label)
        )

        colors: list[str] = ["red", "green", "blue", "gold"]
        bases: list[str] = ["A", "C", "G", "T"]
        colors_bars: list[str] = list(itertools.chain(
            *[sec_len*[color] for color in colors]))

        plt.bar(motifs_names, motifs_vals, color=colors_bars)  # type: ignore
        for num, (color, base) in enumerate(zip(colors, bases)):
            x_start: float = sec_len*num-0.5
            x_end: float = sec_len*(num+1)-0.5
            x_mid: float = (x_start + x_end) / 2
            plt.hlines(y=max_count*1.15, color=color, linewidth=3,
                       xmin=x_start, xmax=x_end)
            plt.text(x_mid, max_count*1.12, base, ha="center", va="top",
                     color=color, fontsize=18)

        plt.title("{}: {}' End motif distribution n={}".format(
            name, label, num_frags))
        plt.xlabel("{}-mer End Motifs".format(kmer_len))
        plt.ylabel("# of Motifs")
        plt.xticks(rotation=80, fontsize=0.6)
        plt.tick_params(axis="x", length=0, pad=0.4)

        fig.savefig(outpath, dpi=1200, bbox_inches="tight")


def log_stats(
    fragments: FragmentList, logger: logging.Logger,
    out_dir: str, name: str
) -> None:
    """
    Generate overview fragment statistics with logging and JSON export.

    Calculates quality metrics for a fragment collection, logs them to the
    console and exports structured data to JSON format for downstream analysis
    and record keeping. The JSON file contains:

    - ``number_of_fragments``: Total fragment count
    - ``number_of_bogus_fragments``: Low-quality fragment count
    - ``number_of_mutated_fragments``: Mutation-carrying fragment count


    Parameters
    ----------
    fragments : FragmentList
        Fragment collection to analyze. All fragments are included in
        statistics, with separate counts for different quality categories.
    logger : logging.Logger
        Logger instance for console output. Used to report statistics.
    out_dir : str
        Output directory for JSON file. Must exist and be writable.
    name : str
        Sample identifier used for JSON filename and log messages.

    Returns
    -------
    A JSON file named ``{name}_frag_stats.json``.
    """
    num_bogus_frags: int = fragments.count_bogus_fragments()
    num_mut_frags: int = fragments.count_mutated_fragments()
    num_frags: int = fragments.length()

    if num_frags == 0:
        fail("no fragments found")

    logger.info(
        "identified {} fragments (bogus: {} ({}%), mutated: {} ({}%))".format(
            num_frags, num_bogus_frags,
            round(100*num_bogus_frags/num_frags, 3),
            num_mut_frags,
            round(100*num_mut_frags/num_frags, 3))
    )

    outpath: str = os.path.join(out_dir, "{}_frag_stats.json".format(name))
    with open(outpath, "w") as file:
        data: dict[str, object] = {
            "number_of_fragments": num_frags,
            "number_of_bogus_fragments": num_bogus_frags,
            "number_of_mutated_fragments": num_mut_frags
        }
        json.dump(data, file)


def export_length_distribution_csv(
    fragments: FragmentList, out_dir: str, name: str
) -> None:
    """
    Export fragment length distribution to CSV file.

    Creates a CSV file containing fragment length counts for statistical
    analysis and visualization in external tools. Each row represents a
    unique fragment length with its corresponding count. The column names are
    as follows:

    - ``fragment_length``: Fragment length in base pairs (sorted ascending)
    - ``count``: Number of fragments with that specific length

    Parameters
    ----------
    fragments : FragmentList
        Fragment collection to analyze. Bogus fragments are automatically
        excluded.
    out_dir : str
        Output directory path where the CSV file will be saved. Directory
        must exist and be writable.
    name : str
        Sample identifier used for the output filename. Used to generate
        unique CSV filenames for each sample.

    Returns
    -------
    A CSV file named ``{name}_length_distribution.csv``.
    """
    logger: logging.Logger = get_logger()
    lengths: list[int] = []
    for frag in fragments:
        if not frag.is_bogus:
            lengths.append(frag.length)

    if not lengths:
        fail("no valid fragments found for length distribution")

    length_counts: Counter[int] = Counter(lengths)
    sorted_lengths: list[tuple[int, int]] = sorted(length_counts.items())

    outpath: str = os.path.join(out_dir, f"{name}_length_distribution.csv")
    with open(outpath, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["fragment_length", "count"])
        for length, count in sorted_lengths:
            writer.writerow([length, count])

    logger.info(
        "exported {} unique fragment lengths to {}".format(
            len(sorted_lengths), outpath
        )
    )


def export_end_motifs_csv(
    fragments: FragmentList, out_dir: str, name: str, kmer_len: int
) -> None:
    """
    Export fragment end motif counts to CSV file.

    Creates a CSV file containing 5' and 3' end motif frequencies for
    downstream analysis and visualization in external tools. Each row
    represents a unique k-mer motif with counts for both fragment ends. The
    column names are as follows:

    - ``motif_5p``: 5' end k-mer sequence (e.g., "AAA", "AAC", etc.)
    - ``count_5p``: Number of fragments with this 5' end motif
    - ``motif_3p``: 3' end k-mer sequence (e.g., "AAA", "AAC", etc.)
    - ``count_3p``: Number of fragments with this 3' end motif

    Parameters
    ----------
    fragments : FragmentList
        Fragment collection for motif analysis. Bogus fragments are
        automatically excluded.
    out_dir : str
        Output directory path where the CSV file will be saved. Directory
        must exist and be writable.
    name : str
        Sample identifier used for the output filename. Used to generate
        unique CSV filenames for each sample.
    kmer_len : int
        Length of k-mer motifs to analyze. Typically 3 or 4 nucleotides.

    Returns
    -------
    A CSV file named ``{name}_k{kmer_len}_end_motifs.csv``.
    """
    logger: logging.Logger = get_logger()
    motifs_5p: defaultdict[str, int]
    motifs_3p: defaultdict[str, int]
    num_frags: int

    motifs_5p, motifs_3p, num_frags = fragments.count_endmotifs(kmer_len)
    if num_frags == 0:
        fail("no valid fragments found for end motif analysis")

    all_motifs: set[str] = set(motifs_5p.keys()) | set(motifs_3p.keys())
    sorted_motifs: list[str] = sorted(all_motifs)

    outpath: str = os.path.join(out_dir, f"{name}_k{kmer_len}_end_motifs.csv")
    with open(outpath, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(
            ["motif_5p", "count_5p", "motif_3p", "count_3p"]
        )
        for motif in sorted_motifs:
            count_5p: int = motifs_5p[motif]
            count_3p: int = motifs_3p[motif]
            writer.writerow(
                [motif, count_5p, motif, count_3p]
            )

    logger.info(
        "exported {} unique k-mer motifs to {}".format(
            len(sorted_motifs), outpath
        )
    )
