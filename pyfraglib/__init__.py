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

from pyfraglib.core import fail, PyfragManager
from pyfraglib.math import fit_gmm, plot_gmm
from pyfraglib.fragment import Fragment, FragmentList, FragmentCollection
from pyfraglib.fragfile import FragFile
from pyfraglib.stats import fragments_per_chromosome_barplot
from pyfraglib.lengths import fragment_length_plot

__version__ = "0.1.0"
__repo_url__ = "https://bitbucket.org/schwarzlab/project-lymphoma-cfdna/"
__all__ = ["Fragment", "FragmentList", "FragmentCollection", "fail",
           "FragFile",
           "fragments_per_chromosome_barplot",
           "fragment_length_plot",
           "fit_gmm", "plot_gmm"]

logging.basicConfig(level=logging.NOTSET)
logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)
logging.getLogger('matplotlib.pyplot').setLevel(logging.ERROR)

PyfragManager.register(  # type: ignore
    "FragmentCollection", FragmentCollection)
