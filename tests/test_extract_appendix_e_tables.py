from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

_SPEC = importlib.util.spec_from_file_location(
    "extract_appendix_e_tables",
    Path(__file__).resolve().parents[1] / "tools" / "extract_appendix_e_tables.py",
)
extractor = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(extractor)

TABLE = "|  ν | k  |   |\n| --- | --- | --- |\n|   | 2 | 3 |\n|  3 | 2.353 | 2.989 |\n|  ∞ | 1.645 | 1.916 |\n"


@pytest.mark.parametrize(
    "heading",
    [
        "Table E.2 Values of  $d$  for constrained MCB and one-sided MCC $\\alpha = 0.05$\n",
        "Table E.2 Values of d for constrained MCB and one-sided MCC\n\nα = 0.05\n\n",
        "Table E.2 Values of d for constrained MCB and one-sided MCC\n\n$$\\alpha = 0.05$$\n\n",
    ],
)
def test_alpha_is_parsed_inline_or_on_its_own_line(heading: str) -> None:
    table = extractor.parse_appendix_e_tables(heading + TABLE)["E.2"]

    assert table["alpha"] == 0.05
    assert table["k_values"] == [2, 3]
    assert table["matrix"] == [[2.353, 2.989], [1.645, 1.916]]
