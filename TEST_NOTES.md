# Verification notes

Last verified in the packaging environment:

- `python -m compileall -q .` — passed.
- `python -m pytest -q tests` — **61 passed**, with one Biopython deprecation warning from an existing FASTA error-handling test.
- New algorithm tests cover input validation, LCS reconstruction, Levenshtein distance and edit operations, Smith–Waterman local alignment, and recursive-vs-DP agreement.
- The Streamlit UI was not launched in the verification environment because Streamlit was not installed there. Install the `web` extra and run `streamlit run app.py` to verify the interactive interface on your machine.

To repeat the tests from the project root:

```bash
python -m pip install -e ".[dev,web]"
python -m pytest -q tests
```
