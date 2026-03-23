"""Sphinx configuration for the mcb documentation."""

from __future__ import annotations

import os
import sys

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

sys.path.insert(0, os.path.abspath("../src"))

from mcb import __version__

project = "mcb"
author = "mcb contributors"
release = __version__

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
]

source_suffix = {
    ".md": "markdown",
}
master_doc = "index"
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

autosummary_generate = False
autodoc_typehints = "description"
autodoc_member_order = "bysource"
napoleon_google_docstring = False
napoleon_numpy_docstring = True

myst_enable_extensions = [
    "amsmath",
    "dollarmath",
]
myst_heading_anchors = 3

html_theme = "sphinx_rtd_theme"
html_title = f"mcb {release} documentation"
html_static_path = []
