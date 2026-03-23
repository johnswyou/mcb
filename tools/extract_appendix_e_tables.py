#!/usr/bin/env python3
"""Extract Appendix E tables from Hsu's OCR markdown into JSON test fixtures."""

from __future__ import annotations

import argparse
import json
import math
import re
from collections import OrderedDict
from pathlib import Path

DEFAULT_OUTPUT_PATH = Path(__file__).resolve().parents[1] / "tests" / "data" / "appendix_e_tables.json"


def parse_pipe_row(line: str) -> list[str]:
    return [p.strip() for p in line.strip().strip("|").split("|")]


def parse_table_block(lines: list[str], start_idx: int):
    i = start_idx + 1
    while i < len(lines) and not lines[i].startswith("|"):
        i += 1
    if i + 2 >= len(lines):
        raise ValueError(f"No markdown table found after line {start_idx+1}")
    header1 = lines[i]
    header2 = lines[i + 1]
    header3 = lines[i + 2]
    data_lines = []
    j = i + 3
    while j < len(lines) and lines[j].startswith("|"):
        data_lines.append(lines[j])
        j += 1
    return header1, header2, header3, data_lines, j


def parse_appendix_e_tables(markdown_text: str) -> OrderedDict[str, dict]:
    lines = markdown_text.splitlines()
    heading_indices = [i for i, l in enumerate(lines) if re.match(r"Table E\.\d+", l)]
    tables = OrderedDict()

    for idx in heading_indices:
        line = lines[idx]
        m = re.match(
            r"Table (E\.\d+)\s+Values of\s+(.*?)\s+for\s+(.*?)(?:\s+\$\\alpha = ([0-9.]+)\$)?\s*$",
            line,
        )
        if not m:
            raise ValueError(f"Could not parse Appendix E heading: {line!r}")

        table_id, quantity, method, alpha = m.group(1), m.group(2).strip(), m.group(3).strip(), m.group(4)
        _, _, header_row, data_lines, _ = parse_table_block(lines, idx)
        k_values = [int(x) for x in parse_pipe_row(header_row)[1:]]

        record = tables.setdefault(
            table_id,
            {
                "table_id": table_id,
                "title": "",
                "quantity": quantity,
                "method": method,
                "alpha": None,
                "nu_values": [],
                "k_values": [],
                "matrix_by_nu": OrderedDict(),
                "source_heading_lines": [],
            },
        )
        record["source_heading_lines"].append(line)
        if alpha is not None:
            record["alpha"] = float(alpha)

        for row_line in data_lines:
            row = parse_pipe_row(row_line)
            nu_label = row[0]
            nu_value = math.inf if "∞" in nu_label else int(nu_label)
            vals = [float(x) for x in row[1:]]

            if nu_value not in record["matrix_by_nu"]:
                record["matrix_by_nu"][nu_value] = {}
                record["nu_values"].append(nu_value)

            for k, v in zip(k_values, vals):
                record["matrix_by_nu"][nu_value][k] = v
            for k in k_values:
                if k not in record["k_values"]:
                    record["k_values"].append(k)

    for record in tables.values():
        alpha_str = f"alpha = {record['alpha']:.2f}" if record["alpha"] is not None else "alpha unknown"
        record["title"] = f"Values of {record['quantity']} for {record['method']} ({alpha_str})"
        record["matrix"] = [
            [record["matrix_by_nu"][nu][k] for k in record["k_values"]]
            for nu in record["nu_values"]
        ]
        del record["matrix_by_nu"]

    return tables


def write_json(tables, path: Path) -> None:
    serializable = {}
    for table_id, rec in tables.items():
        serializable[table_id] = {
            "table_id": table_id,
            "title": rec["title"],
            "quantity": rec["quantity"],
            "method": rec["method"],
            "alpha": rec["alpha"],
            "nu_values": ["inf" if v == math.inf else v for v in rec["nu_values"]],
            "k_values": rec["k_values"],
            "matrix": rec["matrix"],
            "source_heading_lines": rec["source_heading_lines"],
        }
    path.write_text(json.dumps(serializable, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("markdown_path", type=Path)
    parser.add_argument("--json-out", type=Path, default=DEFAULT_OUTPUT_PATH)
    args = parser.parse_args()

    tables = parse_appendix_e_tables(args.markdown_path.read_text())
    write_json(tables, args.json_out)
    print(f"Wrote {len(tables)} Appendix E tables to {args.json_out}.")


if __name__ == "__main__":
    main()
