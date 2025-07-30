"""
Fragment File I/O and Serialization
====================================

This module provides efficient I/O operations for fragment data using
compressed serialization. The FragFile class enables access to serialized
fragment data for convenient iteration.

File Format
-----------
Fragment files (.frag) use the following format:

- **Compression**: Gzip compression (typically 70-80% size reduction)
- **Serialization**: Python pickle protocol for object storage
- **Streaming iteration**: Space-constant streaming support for large files

Example Usage
-------------
Basic fragment file operations:

.. code-block:: python

    from pyfraglib.fragfile import FragFile
    from pyfraglib.fragment import Fragment

    # Create fragments and save to file
    fragments = Fragment.from_bam("sample.bam", "variants.vcf")
    fragments.to_frag_file("sample", "output/")

    # Load fragments from file (streaming!)
    with FragFile("output/sample.frag") as fragfile:
        for fragment in fragfile:
            print(f"Fragment: {fragment.chrom}:{fragment.start_pos}-"
                  f"{fragment.end_pos}")

    # Load all fragments into memory
    with FragFile("output/sample.frag") as fragfile:
        all_fragments = fragfile.get_fragment_list()
        print(f"Loaded {len(all_fragments)} fragments")

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
import pickle
import gzip

from typing import Generator
from pyfraglib import Fragment, FragmentList


class FragFile:
    """
    Efficient reader for compressed fragment files (.frag format).

    Provides access to serialized fragment data stored in compressed pickle
    format. Supports automatic file management and cleanup using the Python
    context manager interface.

    Example:
        Basic usage with context manager::

            from pyfraglib.fragfile import FragFile

            # Streaming iteration (memory efficient)
            with FragFile("sample.frag") as fragfile:
                for fragment in fragfile:
                    print(f"Fragment: {fragment.chrom}:{fragment.start_pos}")

            # Bulk loading for analysis
            with FragFile("sample.frag") as fragfile:
                fragments = fragfile.get_fragment_list()
                print(f"Loaded {len(fragments)} fragments")

    See Also:
        * :meth:`pyfraglib.fragment.FragmentList.to_frag_file` - Writing
          fragment files
        * :class:`pyfraglib.fragment.Fragment` - Individual fragment objects
        * :class:`pyfraglib.fragment.FragmentList` - Fragment collections
    """

    def __init__(self, path: str) -> None:
        """
        Initialize FragFile reader for the specified fragment file.

        Opens the fragment file for reading with gzip decompression.

        Args:
            path: Path to the fragment file (.frag format). Must be a valid
                file path to a gzip-compressed pickle file containing Fragment
                objects.

        Raises:
            FileNotFoundError: If the specified file does not exist.
            gzip.BadGzipFile: If the file is not a valid gzip file.
            PermissionError: If the file cannot be opened due to permissions.

        Note:
            * No validation of file contents is performed during initialization
            * Use context manager for automatic cleanup
        """
        self.__path: str = path
        self.__file: gzip.GzipFile = gzip.open(self.__path, "rb")

    def __iter__(self) -> Generator[Fragment, None, None]:
        """Iterate through fragments in the file using streaming I/O.

        Provides memory-efficient iteration through all fragments in the file
        without loading the entire file into memory. Each fragment is
        deserialized on-demand as the iterator is consumed.

        Yields:
            Fragment: Individual Fragment objects from the file in the order
            they were originally serialized.

        Raises:
            pickle.UnpicklingError: If fragment deserialization fails due to
                corrupted data or version incompatibility.
            gzip.BadGzipFile: If the file format is corrupted.

        Note:
            * Iteration is single-pass - file position advances with each yield
            * Memory usage remains constant regardless of file size
            * Iteration stops automatically when EOF is reached
            * File position is not reset between iterations

        Example:
            Memory-efficient processing of large fragment files::

                from pyfraglib.fragfile import FragFile

                # Process fragments without loading entire file
                fragment_count = 0
                mutated_count = 0

                with FragFile("large_sample.frag") as fragfile:
                    for fragment in fragfile:
                        fragment_count += 1
                        if fragment.is_mutated:
                            mutated_count += 1

                print(f"Processed {fragment_count} fragments")

        See Also:
            * :meth:`get_fragment_list` - Bulk loading all fragments into
              memory
            * :class:`pyfraglib.fragment.Fragment` - Fragment object structure
        """
        while True:
            try:
                frag: Fragment = pickle.load(self.__file)
                yield frag
            except EOFError:
                return None

    def get_fragment_list(self) -> FragmentList:
        """
        Load all fragments from the file into a FragmentList object.

        Reads the entire fragment file and creates a FragmentList containing
        all fragments. This method is convenient for analysis workflows that
        require random access to fragments or need to process all fragments
        multiple times.

        Returns:
            FragmentList: A FragmentList object containing all fragments from
            the file, ready for analysis and manipulation.

        Raises:
            pickle.UnpicklingError: If any fragment deserialization fails.
            gzip.BadGzipFile: If the file format is corrupted.
            MemoryError: If the file is too large to fit in available memory.

        Warning:
            This method loads all fragments into memory simultaneously, which
            can consume significant memory for large files. Consider using
            streaming iteration for memory-constrained environments.

        See Also:
            * :meth:`__iter__` - Memory-efficient streaming iteration
            * :class:`pyfraglib.fragment.FragmentList` - Fragment collection
              class
            * :mod:`pyfraglib.lengths` - Fragment length analysis functions
        """
        fragment: Fragment
        fragment_list: FragmentList = FragmentList()
        for fragment in self:
            fragment_list.append(fragment)
        return fragment_list

    def close(self) -> None:
        """
        Close the fragment file and release system resources.

        Explicitly closes the underlying gzip file handle and releases any
        associated system resources. This method is called automatically by
        the context manager and destructor, but can be called manually when
        needed.
        """
        self.__file.close()

    def __enter__(self) -> "FragFile":
        """
        Context manager entry - returns self for use in with statements.

        Enables FragFile to be used as a context manager, ensuring proper
        resource cleanup even if exceptions occur during processing.

        Returns:
            FragFile: Self reference for use in with statement.
        """
        return self

    def __exit__(
        self, exc_type: object, exc_val: object, exc_tb: object
    ) -> None:
        """
        Context manager exit - ensures file is closed.

        Automatically closes the file when exiting the context manager,
        regardless of whether an exception occurred during processing.

        Args:
            exc_type: Exception type (if any) that caused context exit.
            exc_val: Exception value (if any) that caused context exit.
            exc_tb: Exception traceback (if any) that caused context exit.
        """
        self.close()

    def __del__(self) -> None:
        """
        Destructor that ensures file closing when object is garbage collected.

        Provides automatic cleanup of file resources when the FragFile object
        is destroyed. This serves as a safety mechanism to prevent resource
        leaks, but explicit cleanup using context managers is preferred.
        """
        self.close()
