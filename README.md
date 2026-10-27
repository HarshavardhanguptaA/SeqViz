# DNA Sequence Alignment Tool

An educational bioinformatics application for comparing DNA sequences using dynamic programming. It includes a Streamlit web interface and retains the upstream SeqViz command-line alignment tools.

## Features

- **Longest Common Subsequence (LCS):** recover a common subsequence and inspect its dynamic-programming matrix.
- **Levenshtein edit distance:** calculate the minimum number of insertions, deletions, and substitutions; show the operation trace and DP matrix.
- **Smith–Waterman local alignment:** find high-scoring local regions using configurable match, mismatch, and gap scores.
- **Alignment visualization:** color-code matches, substitutions, and gaps.
- **Mutation simulator:** apply a substitution, insertion, or deletion to sequence 1 and compare the resulting edit distance.
- **Similarity percentage:** report a normalized edit-distance score.
- **Algorithm benchmark:** compare dynamic programming with naive recursive edit distance on short inputs.
- **Existing CLI workflow:** use Needleman–Wunsch global alignment, Smith–Waterman local alignment, FASTA input, and heatmap output through the `seqviz` command.

## Requirements

- Python 3.11 or newer
- Internet access for initial package installation

## Install and run (Windows PowerShell)

From the extracted project root:

```powershell
python -m venv .venv
.venv\\Scripts\\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -e ".[dev,web]"
streamlit run app.py
```

If PowerShell blocks virtual-environment activation, run the executable directly instead:

```powershell
.venv\\Scripts\\python.exe -m pip install -e ".[dev,web]"
.venv\\Scripts\\python.exe -m streamlit run app.py
```

The app opens in your browser. Enter two sequences in the sidebar and use the tabs to explore comparison results, local alignment, mutations, and performance.

### Install and run (macOS / Linux)

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e ".[dev,web]"
streamlit run app.py
```

## Command-line interface

The original `seqviz` CLI remains available after installing the project.

```bash
seqviz --help
seqviz --seq1 GATTACA --seq2 GCATGCU --algorithm nw
seqviz --seq1 ACACACTA --seq2 AGCACACA --algorithm sw
```

For FASTA input and heatmap output, use the example files under `examples/`. Consult `seqviz --help` for the complete options list.

## Run tests

```bash
python -m pytest
```

The tests cover DNA validation, LCS reconstruction, Levenshtein distance, Smith–Waterman alignment, and agreement between recursive and dynamic-programming edit distance on short strings.

## Algorithms and complexity

Let `m` and `n` be the lengths of the two sequences.

| Algorithm | Time | Space | Purpose |
|---|---:|---:|---|
| LCS dynamic programming | O(mn) | O(mn) | Longest shared subsequence |
| Levenshtein dynamic programming | O(mn) | O(mn) | Minimum insertion/deletion/substitution count |
| Smith–Waterman | O(mn) | O(mn) | Best local alignment under a scoring scheme |
| Naive recursive edit distance | Exponential in the worst case | O(m+n) recursion depth | Small-input comparison baseline |

The naive-recursion tab limits sequence length because the method grows exponentially and can become very slow.

## Input and interpretation notes

- The interactive interface accepts `A`, `C`, `G`, `T`, and `N` and ignores whitespace.
- The interactive sequence length is capped at 500 bases to keep matrix rendering manageable.
- The similarity score is defined as `100 × (1 − edit_distance / max(lengths))`. It is a normalized edit-distance score, not a clinical diagnosis or a universal measure of evolutionary similarity.
- Smith–Waterman scores depend on the selected scoring parameters.
- This is an educational project and is not intended for clinical decision-making or validated genomic analysis.

## Project structure

```text
.
├── app.py                         # Streamlit interface
├── pyproject.toml                 # Package metadata and dependencies
├── requirements-app.txt           # Web app dependencies (alternative install)
├── src/
│   └── seqviz/
│       ├── dna_algorithms.py      # LCS, edit distance, Smith–Waterman
│       ├── algorithms/            # Existing alignment implementations
│       ├── visualise/              # Existing terminal and heatmap output
│       ├── cli.py                 # Command-line entry point
│       └── fasta_utils.py         # FASTA helpers
├── tests/                          # Unit and regression tests
├── examples/                       # Sample FASTA sequences
└── LICENSE
```

## Credits and license

This repository extends **SeqViz**, originally published by Ryan Young / BarefootRobber. The upstream implementation, example data, and MIT license are retained. The interactive Streamlit interface and the additional educational algorithms/tests are extensions in this project. Please retain the `LICENSE` file and this attribution when redistributing the repository.

Upstream project: https://github.com/BarefootRobber/SeqViz

The software is provided under the terms in `LICENSE`. Review the license and any third-party data/source terms before redistribution.
