from __future__ import annotations

import importlib.util
from pathlib import Path

_SPEC = importlib.util.spec_from_file_location(
    "make_sample_size_figures",
    Path(__file__).resolve().parents[1] / "tools" / "make_sample_size_figures.py",
)
figures = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(figures)


def test_figures_still_show_what_the_guide_narrates(tmp_path):
    # Generating asserts that the cut line matches constrained_mcb and that each narrated outcome holds.
    figures.main(tmp_path)

    assert sorted(path.name for path in tmp_path.iterdir()) == [
        "sample_size_luck_in_one_test.svg",
        "sample_size_luck_needed.svg",
        "sample_size_many_tests.svg",
        "sample_size_more_users.svg",
        "sample_size_one_test.svg",
    ]
