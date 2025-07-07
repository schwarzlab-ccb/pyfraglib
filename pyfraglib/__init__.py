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

from pyfraglib.core import get_logger, fail, PyfragManager
from pyfraglib.math import fit_gmm, plot_gmm
from pyfraglib.fragment import Fragment, FragmentList, FragmentCollection
from pyfraglib.fragfile import FragFile
from pyfraglib.stats import fragments_per_chromosome_barplot
from pyfraglib.lengths import fragment_length_plot

__version__ = "0.4.3"
__repo_url__ = "https://bitbucket.org/schwarzlab/project-lymphoma-cfdna/"
__all__ = [
    "Fragment", "FragmentList", "FragmentCollection", "fail",
    "get_logger", "FragFile", "fragments_per_chromosome_barplot",
    "fragment_length_plot", "fit_gmm", "plot_gmm"
]


# @NOTE(ds): This filter changes the level of every record to `level' and only
# logs the record if that new level is more severe than the log level of
# the logger that emitted the record.
class FixedLogLevelFilter(logging.Filter):
    def __init__(self, level: int):
        super().__init__()
        self.level: int = level

    def filter(self, record: logging.LogRecord) -> bool:
        this_logger_name: str = record.name
        this_logger: logging.Logger = logging.getLogger(this_logger_name)

        record.levelno = self.level
        record.levelname = logging.getLevelName(self.level)
        return record.levelno >= this_logger.level


logging.basicConfig(
    level=logging.NOTSET,
    format='[%(asctime)s %(levelname)-8s %(name)-9s %(process)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logging.captureWarnings(True)

# @NOTE(ds): The `matplotlib' logging messages are just annoying. We turn them
# off completely.
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
logging.getLogger("matplotlib.pyplot").setLevel(logging.ERROR)
logging.getLogger("PIL.PngImagePlugin").setLevel(logging.ERROR)

# @NOTE(ds): Not sure if the following log messages could be interesting for
# debugging.
logging.getLogger("asyncio").setLevel(logging.DEBUG)

# @NOTE(ds): `py.warnings' logs everything as a WARNING. We don't want that
# so we apply a filter to log everything as DEBUG.
py_warnings_logger: logging.Logger = logging.getLogger("py.warnings")
log_filter: FixedLogLevelFilter = FixedLogLevelFilter(logging.DEBUG)
py_warnings_logger.addFilter(log_filter)

PyfragManager.register(  # type: ignore
    "FragmentCollection", FragmentCollection
)
