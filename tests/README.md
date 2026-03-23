
The test suite uses `pytest` and is organized around the installable `mcb` package.

Contents:

* `data/appendix_e_tables.json` — canonical Appendix E fixture data used by the regression tests
* `test_appendix_e_critical_values.py` — balanced-case checks for Tables E.2, E.3, E.5, and E.6
* `test_core.py` — summary validation, synthetic result cases, and output formatting checks
* `test_plotting.py` — smoke test for the optional plotting API

The Appendix E fixture can be regenerated from the OCR markdown with:

```bash
python3 tools/extract_appendix_e_tables.py multiple_comparisons_hsu_1996/markdown.md
```

Run the tests with:

```bash
python3 -m pip install -e '.[test]'
python3 -m pytest
```
