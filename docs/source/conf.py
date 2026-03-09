# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html
#
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
import os
import sys

sys.path.insert(0, os.path.abspath("../../"))

# General configuration information:
# https://www.sphinx-doc.org/en/master/usage/configuration.html
project = "pyfraglib"
copyright = "2025, Daniel Schütte"
author = "Daniel Schütte"
release = "0.5.1"
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.todo",
    "sphinx.ext.coverage",
    "sphinx.ext.doctest",
    "sphinx.ext.githubpages",
    "sphinx_autodoc_typehints",
    "myst_parser",
    "sphinx_rtd_theme",
    "rst2pdf.pdfbuilder",
]

napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = False
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_type_aliases = None

autodoc_default_options = {
    "members": True,
    "member-order": "bysource",
    "special-members": "__init__",
    "undoc-members": True,
    "exclude-members": "__weakref__"
}

autodoc_type_aliases = {
    "SequenceContextGenerator": "pyfraglib.simulator.SequenceContextGenerator"
}
suppress_warnings = ["ref.python"]
autodoc_mock_imports: list[object] = []  # for missing modules
autodoc_warningiserror = False  # don't show warnings about missing modules
keep_warnings = True  # keep going on warnings

autosummary_generate = True
autosummary_generate_overwrite = True

typehints_defaults = "comma"
typehints_use_signature = True
typehints_use_signature_return = True
always_document_param_types = True

intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "matplotlib": ("https://matplotlib.org/stable/", None),
    "pandas": ("https://pandas.pydata.org/docs/", None),
    "pysam": ("https://pysam.readthedocs.io/en/latest/", None),
}

myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "dollarmath",
    "html_admonition",
    "html_image",
    "linkify",
    "replacements",
    "smartquotes",
    "substitution",
    "tasklist",
]

templates_path = ["_templates"]
exclude_patterns: list[object] = []

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_theme_options = {
    "logo_only": False,
    "prev_next_buttons_location": "bottom",
    "style_external_links": False,
    "vcs_pageview_mode": "",
    "style_nav_header_background": "#343131",
    "collapse_navigation": True,
    "sticky_navigation": True,
    "navigation_depth": 4,
    "includehidden": True,
    "titles_only": False
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_logo = "../imgs/pyfraglib_logo.jpg"
html_favicon = "../imgs/pyfraglib_logo.jpg"
todo_include_todos = True

source_suffix = {
    ".rst": None,
    ".md": "myst_parser",
}
master_doc = "index"  # the master toctree document
html_use_modindex = True
html_use_index = True
html_search_language = "en"
html_show_sourcelink = True
html_copy_source = True

latex_engine = "pdflatex"  # xelatex
latex_elements = {
    "papersize": "letterpaper",  # or a4paper
    "pointsize": "10pt",
    "preamble": r"""
\usepackage{newunicodechar}
\usepackage{amsmath,amsfonts,amssymb,amsthm}
\usepackage{graphicx}
\usepackage{hyperref}
\usepackage{xcolor}
\definecolor{VerbatimColor}{rgb}{0.95,0.95,0.95}
\definecolor{VerbatimBorderColor}{rgb}{0.8,0.8,0.8}
% Disable footnotes for hyperlinks
\hypersetup{
    colorlinks=true,
    linkcolor=blue,
    urlcolor=blue,
    citecolor=blue,
    filecolor=blue
}
""",
    "figure_align": "htbp",
    "geometry": r"\usepackage[margin=1in]{geometry}",
}

latex_documents = [
    (master_doc, "pyfraglib_documentation.tex", "pyfraglib Documentation",
     "Daniel Schütte", "manual"),
]
latex_logo = "../imgs/pyfraglib_logo.jpg"
latex_use_parts = False
latex_show_pagerefs = False
latex_show_urls = "no"
latex_appendices: list[object] = []
latex_domain_indices = True

man_pages = [
    (master_doc, "pyfraglib", "pyfraglib Documentation",
     [author], 1)
]
texinfo_documents = [
    (master_doc, "pyfraglib", "pyfraglib Documentation",
     author, "pyfraglib", "Python library for cfDNA fragmentomics analysis.",
     "Miscellaneous"),
]

epub_title = project
epub_author = author
epub_publisher = author
epub_copyright = copyright
epub_identifier = "pyfraglib"
epub_uid = "pyfraglib"
epub_exclude_files = ["search.html"]
