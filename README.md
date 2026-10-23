# Circumcoronene and Face-Stellated Graph Indices

Computational code and reproducibility data for circumcoronene graphs (`Cor`, or H_n) and their face-stellated counterparts (`St`, or T_n). The scripts construct both graph families, compute 19 graph descriptors, and compare graph-based calculations with coordinate formulas and symbolic identities.

This repository intentionally excludes the manuscript, LaTeX files, and rendered figures. Figure-generation code is included and can generate the figures locally. No manuscript or external dataset is needed to run the code.

## Contents

| File | Purpose |
| --- | --- |
| `verify.py` | Graph construction, all-source BFS, index calculations, and numerical verification |
| `symbolic_checks.py` | Symbolic formulas and algebraic certificates |
| `make_figures.py` | Optional figure generation from code and CSV data |
| `requirements.txt` | Pinned Python dependencies |
| `indices.csv` | Descriptor values for both families, n = 1 through 15 |
| `distance_partitions.csv` | Counts of unordered vertex pairs at each positive distance |
| `verification.csv` | Formula-versus-direct-calculation comparisons |
| `old_vertex_comparison.csv` | Changes in distances between original vertices after face stellation |
| `verification_summary.json` | Verification counts, status, and recorded runtime versions |
| `DATA_DICTIONARY.md` | Column meanings and index conventions |
| `.gitignore` | Excludes virtual environments, caches, manuscript files, and generated figures |
| `manifest_sha256.json` | SHA-256 checksums of the distributed files, excluding the manifest itself |

## Setup

Tested with Python 3.13.3. Install the pinned versions in an isolated environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

On Windows PowerShell, use `.venv\Scripts\Activate.ps1` instead of the `source` command.

## Reproduce Results

Run from the repository directory:

```bash
python verify.py --max-n 15
```

This rebuilds both families for n = 1 through 15 and overwrites the four CSV files and `verification_summary.json`. Run without Python's `-O` option: the verification checks use assertions. A failed check raises an exception and the process exits unsuccessfully.

The distributed run covers 30 graphs and 19 descriptors per graph, with 405 formula comparisons and 24 symbolic certificates. It also checks coordinate distances against BFS for 10,049,712 unordered vertex pairs. Integer and rational distance calculations are exact; radical-valued degree indices and the logarithmic product representation use floating-point arithmetic. Degree formula comparisons allow relative tolerance 1e-12 and absolute tolerance 1e-10.

For symbolic checks only:

```bash
python symbolic_checks.py
```

A smaller `--max-n` can be used for a quick check, but replaces the distributed datasets with that smaller range. Run the full command again to restore the n = 1 through 15 data. Larger values require more memory and time because distance matrices grow with the square of the vertex count. Runtime-version fields in the summary can differ on another machine.

## Optional Figures

After the full verification run:

```bash
python make_figures.py
```

This creates `graphs.pdf`, `normalized_indices.pdf`, `parity.pdf`, and `distance_distribution.pdf`. These files are not distributed and are ignored by Git. The distance-distribution plot requires n = 5 in the CSV data. No graphical display is required.

## Interpretation

The finite computations check the implemented formulas over the stated range; they are not a proof for all n. Symbolic certificates check the algebra of selected expressions. Not every descriptor is represented by a polynomial closed form: some are evaluated from distance distributions or finite sums. The work concerns graph descriptors, not experimental chemical measurements or validated predictions of chemical properties.

## Uploading to Git

Upload all extracted repository files, including `.gitignore`. Do not upload the ZIP itself as the repository source. Keep the manuscript, rendered figures, virtual environment, and cache files outside version control. No remote repository is configured by this package.

No software license has been selected in this package. The rights holders should choose a license before authorizing reuse; public visibility alone does not grant a reuse license.
