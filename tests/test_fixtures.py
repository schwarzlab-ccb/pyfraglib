# This file is part of `pyfraglib`, a software suite to calculate fragmentomics
# features from cfDNA and perform downstream analyses.
#
# Copyright (C) 2026 Daniel Schütte, daniel.schuette@iccb-cologne.org
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version. This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
# more details. You should have received a copy of the GNU General Public
# License along with this program. If not, see <https://www.gnu.org/licenses/>.
import tempfile
import os

from unittest.mock import Mock
from pyfraglib.fragment import Fragment, FragmentList


class MockAlignedSegment:
    """Mock pysam.AlignedSegment for testing purposes."""

    def __init__(
        self,
        query_name: str = "test_read",
        reference_start: int = 100,
        reference_end: int = 250,
        reference_name: str = "1",
        query_sequence: str | None = "ATCGATCGATCG",
        is_forward: bool = True,
        is_read1: bool = True,
        is_read2: bool = False,
        is_reverse: bool = False,
        is_unmapped: bool = False,
        is_duplicate: bool = False,
        mapping_quality: int = 30,
        template_length: int = 150
    ) -> None:
        self.query_name = query_name
        self.reference_start = reference_start
        self.reference_end = reference_end
        self.reference_name = reference_name
        self.query_sequence = query_sequence
        self.is_forward = is_forward
        self.is_read1 = is_read1
        self.is_read2 = is_read2
        self.is_reverse = is_reverse
        self.is_unmapped = is_unmapped
        self.is_duplicate = is_duplicate
        self.mapping_quality = mapping_quality
        self.template_length = template_length

    def get_reference_positions(self, full_length: bool = False) -> list[int]:
        """Mock method to return reference positions."""
        return list(range(self.reference_start, self.reference_end))

    def get_tag(self, tag: str) -> str:
        """Mock method for getting BAM tags."""
        if tag in ["XI", "XJ"]:
            return "mock_value"
        raise KeyError(f"Tag {tag} not found")


def create_mock_fragment(
    start_pos: int = 100,
    end_pos: int = 250,
    chrom: str = "1",
    length: int = 150,
    end5p: str = "ATCG",
    end3p: str = "CGAT",
    is_bogus: bool = False,
    is_mutated: bool | None = None
) -> Fragment:
    """Create a mock Fragment for testing."""
    fragment = Fragment.__new__(Fragment)
    fragment.start_pos = start_pos
    fragment.end_pos = end_pos
    fragment.chrom = chrom
    fragment.length = length
    fragment.end5p = end5p
    fragment.end3p = end3p
    fragment.is_bogus = is_bogus
    fragment.is_mutated = is_mutated
    fragment.is_forward = True
    return fragment


def create_test_fragment_list(num_fragments: int = 10) -> FragmentList:
    """Create a FragmentList with test fragments."""
    fragment_list = FragmentList()

    for i in range(num_fragments):
        fragment = create_mock_fragment(
            start_pos=100 + i * 200,
            end_pos=250 + i * 200,
            chrom=str((i % 22) + 1),
            length=150 + i * 10,
            is_bogus=(i % 5 == 0),
            is_mutated=(i % 3 == 0)
        )
        fragment_list.append(fragment)

    return fragment_list


def create_temp_frag_file(fragments: FragmentList) -> str:
    """Create a temporary .frag file for testing."""
    tmp_dir = tempfile.mkdtemp()
    name = "test_fragments"
    fragments.to_frag_file(name, tmp_dir)
    return os.path.join(tmp_dir, name + ".frag")


def create_mock_bam_header() -> dict[str, object]:
    """Create a mock BAM header for testing."""
    return {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [
            {"SN": "1", "LN": 249250621},
            {"SN": "2", "LN": 243199373},
            {"SN": "X", "LN": 155270560}
        ]
    }


def create_mock_variant_record(
    contig: str = "1",
    start: int = 125,
    stop: int = 126,
    ref: str = "A",
    alts: tuple[str, ...] = ("T",),
    rlen: int = 1
) -> Mock:
    """Create a mock VCF variant record."""
    variant = Mock()
    variant.contig = contig
    variant.start = start
    variant.stop = stop
    variant.ref = ref
    variant.alts = alts
    variant.rlen = rlen
    return variant


def create_mock_vcf_file(variants: list[Mock]) -> Mock:
    """Create a mock VCF file with specified variants."""
    vcf_file = Mock()
    vcf_file.filename = b"test.vcf"
    vcf_file.__iter__ = Mock(return_value=iter(variants))
    return vcf_file


def create_test_bam_with_reads(reads: list[MockAlignedSegment]) -> Mock:
    """Create a mock BAM file with specified reads."""
    bam_file = Mock()
    bam_file.header.to_dict.return_value = {"SQ": [{"SN": "chr1"}]}
    bam_file.fetch.return_value = reads
    bam_file.has_index.return_value = True
    return bam_file


def cleanup_temp_file(filepath: str) -> None:
    """Clean up a temporary test file and, if it lives inside a
    ``tempfile.mkdtemp`` directory, remove that directory too."""
    if os.path.exists(filepath):
        os.unlink(filepath)
    parent = os.path.dirname(filepath)
    if parent and parent.startswith(tempfile.gettempdir()) and \
            os.path.isdir(parent) and not os.listdir(parent):
        os.rmdir(parent)
