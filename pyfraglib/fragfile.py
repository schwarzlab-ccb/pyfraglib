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
import pickle
import gzip

from typing import Generator
from pyfraglib import Fragment, FragmentList


class FragFile:
    def __init__(self, path: str) -> None:
        self.__path: str = path
        self.__file: gzip.GzipFile = gzip.open(self.__path, "rb")

    def __iter__(self) -> Generator[Fragment, None, None]:
        while True:
            try:
                frag: Fragment = pickle.load(self.__file)
                yield frag
            except EOFError:
                return None

    def get_fragment_list(self) -> FragmentList:
        fragment: Fragment
        fragment_list: FragmentList = FragmentList()
        for fragment in self:
            fragment_list.append(fragment)
        return fragment_list

    def close(self) -> None:
        self.__file.close()

    def __del__(self) -> None:
        self.close()
