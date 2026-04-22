"""
Fragment File I/O and Serialization
====================================

This module provides efficient I/O operations for fragment data using a
columnar on-disk format. The FragFile class enables streaming access to
serialized fragment data stored in Apache Parquet files.

File Format
-----------
Fragment files (``.frag``) use Apache Parquet (via :mod:`pyarrow`) with one
row per fragment and the following columns:

.. list-table::
   :header-rows: 1
   :widths: 20 15 65

   * - Column
     - Arrow dtype
     - Description
   * - ``start_pos``
     - int64
     - First aligned position (0-based).
   * - ``end_pos``
     - int64
     - One-past-the-last aligned position.
   * - ``length``
     - int64
     - Fragment length in base pairs.
   * - ``chrom``
     - string
     - Contig name (dictionary-encoded on disk).
   * - ``is_forward``
     - bool
     - Strand / orientation flag.
   * - ``end5p``
     - string
     - 5' end motif (up to ``MAX_KMER_LEN`` nt).
   * - ``end3p``
     - string
     - 3' end motif (up to ``MAX_KMER_LEN`` nt).
   * - ``is_bogus``
     - bool
     - Quality flag; ``True`` marks excluded fragments.
   * - ``is_mutated``
     - bool
     - Mutation-carrying flag; nullable.

Example Usage
-------------
Basic fragment file operations::

    from pyfraglib.fragfile import FragFile
    from pyfraglib.fragment import Fragment

    # Create fragments and save to file
    fragments = Fragment.from_bam("sample.bam", "variants.vcf")
    fragments.to_frag_file("sample", "output/")

    # Load fragments from file (streaming)
    with FragFile("output/sample.frag") as fragfile:
        for fragment in fragfile:
            print(f"{fragment.chrom}:{fragment.start_pos}-"
                  f"{fragment.end_pos}")

    # Load all fragments into memory
    with FragFile("output/sample.frag") as fragfile:
        all_fragments = fragfile.get_fragment_list()
        print(f"Loaded {len(all_fragments)} fragments")

License
-------
This file is part of ``pyfraglib``, a software suite to calculate
fragmentomics features from cfDNA and perform downstream analyses.

Copyright (C) 2026 Daniel Schütte, daniel.schuette@iccb-cologne.org

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version. This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
more details. You should have received a copy of the GNU General Public
License along with this program. If not, see <https://www.gnu.org/licenses/>.
"""
from typing import Generator

import pyarrow as pa
import pyarrow.parquet as pq

from pyfraglib import Fragment, FragmentList, FRAG_SCHEMA


def fragment_from_row(
    start_pos: int, end_pos: int, length: int, chrom: str,
    is_forward: bool, end5p: str, end3p: str, is_bogus: bool,
    is_mutated: bool | None,
) -> Fragment:
    """
    Reconstruct a :class:`pyfraglib.fragment.Fragment` from a row of raw
    field values read from disk.

    Bypasses ``Fragment.__init__`` (which expects ``pysam.AlignedSegment``
    objects) by using ``Fragment.__new__`` followed by explicit attribute
    assignment. ``Fragment`` uses ``@dataclass(slots=True)``, so only the
    declared fields can be set; unknown fields raise ``AttributeError``.

    Args:
        start_pos: First aligned position (0-based).
        end_pos: One-past-the-last aligned position.
        length: Fragment length in base pairs.
        chrom: Contig name.
        is_forward: Strand / orientation flag.
        end5p: 5' end motif.
        end3p: 3' end motif.
        is_bogus: Quality flag.
        is_mutated: Mutation flag (may be ``None``).

    Returns:
        Fragment: A fully populated ``Fragment`` instance.
    """
    fragment: Fragment = Fragment.__new__(Fragment)
    fragment.start_pos = int(start_pos)
    fragment.end_pos = int(end_pos)
    fragment.length = int(length)
    fragment.chrom = str(chrom)
    fragment.is_forward = bool(is_forward)
    fragment.end5p = str(end5p)
    fragment.end3p = str(end3p)
    fragment.is_bogus = bool(is_bogus)
    fragment.is_mutated = None if is_mutated is None else bool(is_mutated)
    return fragment


class FragFile:
    """
    Reader for fragment files (``.frag``) stored as Apache Parquet.

    Provides streaming and bulk access to serialized fragment data with
    context-manager semantics for automatic cleanup. The underlying
    storage is a :class:`pyarrow.parquet.ParquetFile` and reads go through
    :meth:`pyarrow.parquet.ParquetFile.iter_batches` for memory-efficient
    iteration.

    Example:
        Basic usage with context manager::

            from pyfraglib.fragfile import FragFile

            # Streaming iteration (memory efficient)
            with FragFile("sample.frag") as fragfile:
                for fragment in fragfile:
                    print(f"{fragment.chrom}:{fragment.start_pos}")

            # Bulk loading for analysis
            with FragFile("sample.frag") as fragfile:
                fragments = fragfile.get_fragment_list()
                print(f"Loaded {len(fragments)} fragments")

    See Also:
        * :meth:`pyfraglib.fragment.FragmentList.to_frag_file` - Writing
          fragment files.
        * :class:`pyfraglib.fragment.Fragment` - Individual fragment objects.
        * :class:`pyfraglib.fragment.FragmentList` - Fragment collections.
        * :data:`FRAG_SCHEMA` - On-disk column schema.
    """

    #: Batch size used by :meth:`__iter__` when pulling rows from Parquet.
    ITER_BATCH_SIZE: int = 65_536

    def __init__(self, path: str) -> None:
        """
        Initialize a FragFile reader for the specified fragment file.

        Opens the Parquet file lazily — metadata is read but fragment rows
        are not loaded until iteration or bulk loading is invoked.

        Args:
            path: Path to a ``.frag`` Parquet file produced by
                :meth:`FragmentList.to_frag_file`.

        Raises:
            FileNotFoundError: If the specified file does not exist.
            pyarrow.lib.ArrowInvalid: If the file is not a valid Parquet
                file or does not carry the expected fragment columns.
            PermissionError: If the file cannot be opened due to
                permissions.

        Note:
            * Column presence is validated on open; extra columns beyond
              the fragment schema are silently ignored.
            * Use a context manager for automatic cleanup.
        """
        self.__path: str = path
        self.__parquet: pq.ParquetFile | None = pq.ParquetFile(path)
        self.__closed: bool = False

        required: set[str] = {f.name for f in FRAG_SCHEMA}
        present: set[str] = set(self.__parquet.schema_arrow.names)
        missing: set[str] = required - present
        if missing:
            self.close()
            raise pa.lib.ArrowInvalid(
                f"{path!r} is missing fragment columns: {sorted(missing)}"
            )

    def __iter__(self) -> Generator[Fragment, None, None]:
        """
        Iterate through fragments in the file using streaming row batches.

        Yields fragments in the order they appear on disk, pulling
        ``ITER_BATCH_SIZE``-row batches from Parquet at a time so that
        memory usage stays bounded regardless of file size. Each batch is
        converted to Python via ``to_pylist`` and the rows are then
        reconstructed into :class:`Fragment` objects.

        Yields:
            Fragment: Individual Fragment objects from the file.

        Raises:
            ValueError: If called after :meth:`close`.
            pyarrow.lib.ArrowInvalid: If the Parquet stream is corrupted.

        Note:
            * Iteration is single-pass per call; calling ``iter(fragfile)``
              again starts a fresh read from the first row group.
            * Column selection is fixed to the full fragment schema.

        See Also:
            * :meth:`get_fragment_list` - Bulk loading all fragments into
              memory.
            * :class:`pyfraglib.fragment.Fragment` - Fragment object
              structure.
        """
        if self.__closed or self.__parquet is None:
            raise ValueError("FragFile is closed")

        columns: list[str] = [f.name for f in FRAG_SCHEMA]
        batch: pa.RecordBatch
        for batch in self.__parquet.iter_batches(
            batch_size=self.ITER_BATCH_SIZE, columns=columns,
        ):
            start_cols = batch.column("start_pos").to_pylist()
            end_cols = batch.column("end_pos").to_pylist()
            length_cols = batch.column("length").to_pylist()
            chrom_cols = batch.column("chrom").to_pylist()
            forward_cols = batch.column("is_forward").to_pylist()
            e5_cols = batch.column("end5p").to_pylist()
            e3_cols = batch.column("end3p").to_pylist()
            bogus_cols = batch.column("is_bogus").to_pylist()
            mutated_cols = batch.column("is_mutated").to_pylist()
            for s, e, l, c, fw, e5, e3, b, m in zip(
                start_cols, end_cols, length_cols, chrom_cols,
                forward_cols, e5_cols, e3_cols, bogus_cols, mutated_cols,
            ):
                yield fragment_from_row(s, e, l, c, fw, e5, e3, b, m)

    def get_fragment_list(self) -> FragmentList:
        """
        Load all fragments from the file into a :class:`FragmentList`.

        Reads the entire Parquet table in one shot (more efficient than
        streaming for bulk loads) and loads each row as a
        ``Fragment``. Use :meth:`__iter__` for very large files.

        Returns:
            FragmentList: A FragmentList containing every fragment from
            the file.

        Raises:
            ValueError: If called after :meth:`close`.
            MemoryError: If the file is too large to fit in memory.

        Warning:
            This method loads all fragments into memory simultaneously.
            Prefer streaming iteration for memory-constrained environments.

        See Also:
            * :meth:`__iter__` - Memory-efficient streaming iteration.
            * :class:`pyfraglib.fragment.FragmentList` - Fragment
              collection class.
        """
        if self.__closed or self.__parquet is None:
            raise ValueError("FragFile is closed")

        fragment_list: FragmentList = FragmentList()
        columns: list[str] = [f.name for f in FRAG_SCHEMA]
        tbl: pa.Table = self.__parquet.read(columns=columns)

        start_cols = tbl.column("start_pos").to_pylist()
        end_cols = tbl.column("end_pos").to_pylist()
        length_cols = tbl.column("length").to_pylist()
        chrom_cols = tbl.column("chrom").to_pylist()
        forward_cols = tbl.column("is_forward").to_pylist()
        e5_cols = tbl.column("end5p").to_pylist()
        e3_cols = tbl.column("end3p").to_pylist()
        bogus_cols = tbl.column("is_bogus").to_pylist()
        mutated_cols = tbl.column("is_mutated").to_pylist()
        for s, e, l, c, fw, e5, e3, b, m in zip(
            start_cols, end_cols, length_cols, chrom_cols,
            forward_cols, e5_cols, e3_cols, bogus_cols, mutated_cols,
        ):
            fragment_list.append(
                fragment_from_row(s, e, l, c, fw, e5, e3, b, m)
            )
        return fragment_list

    def close(self) -> None:
        """
        Close the fragment file and release system resources.

        Drops the underlying :class:`pyarrow.parquet.ParquetFile` handle.
        Called automatically by the context manager and destructor, but
        may be called manually. Subsequent reads raise ``ValueError``.
        """
        if self.__closed:
            return
        self.__closed = True
        self.__parquet = None  # hint to GC

    @property
    def closed(self) -> bool:
        """Whether :meth:`close` has been invoked on this file."""
        return self.__closed

    def __enter__(self) -> "FragFile":
        """Context-manager entry — returns self."""
        return self

    def __exit__(
        self, exc_type: object, exc_val: object, exc_tb: object
    ) -> None:
        """Context-manager exit — ensures the file is closed."""
        self.close()

    def __del__(self) -> None:
        """Destructor — closes the file if still open."""
        try:
            self.close()
        except Exception:
            pass
