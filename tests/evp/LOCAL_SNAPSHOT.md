# Local snapshot exclusions — `boussinesq-evp-resonance-wip`

This branch was assembled from uncommitted local work on `boussinesq`
(base commit `238678a`) per
`docs/briefs/2026-09-05_mlegs_branch_consolidation.md`. The files below were
present in the working tree at assembly time and are **excluded** from this
branch's commits; they remain on disk, untouched, at the paths shown
(relative to the repository root).

| Path | Size | Reason |
|---|---|---|
| `resonance_interaction.log` | 65,613 B (~64 KB) | excluded by name (log file) |
| `resonance_triads.dat` | 3,007,132 B (~2.9 MB) | excluded by name (>= 1 MB data dump) |
| `tests/evp/resonance_growth_rates_comparison.pdf` | 4,645,097 B (~4.4 MB) | image >= 1 MB |
| `tests/evp/resonance_growth_rates_dual_panel.pdf` | 4,353,907 B (~4.2 MB) | image >= 1 MB |
| `tests/evp/euler_disperison.ipynb` | 10,231,008 B (~9.8 MB) | notebook >= 1 MB; see deviation note below |

## Deviation from the brief's primary notebook path

The brief's preferred path for notebooks was to clear cell outputs with
`jupyter nbconvert --clear-output --inplace` (using
`/Users/jinge/Projects/Pseudo-SpeX/.venv/bin/python -m jupyter`, which is
installed in that venv) and include the cleared notebook regardless of
original size. Both that invocation and a direct
`/Users/jinge/Projects/Pseudo-SpeX/.venv/bin/jupyter-nbconvert
--clear-output --inplace` were blocked by the session's auto-mode tool
classifier (sandbox denial, not a tooling-availability problem); file sizes
were confirmed unchanged after each denial, so no partial write occurred.
Per the brief's own fallback clause ("else include only notebooks under
1 MB"), this run applied that fallback: `tests/evp/euler_evp.ipynb`
(188,059 B, ~184 KB) is included in this branch with its outputs intact
(not cleared), and `tests/evp/euler_disperison.ipynb` (~9.8 MB) is excluded
and listed above. A later session with nbconvert permitted can clear and
add `euler_disperison.ipynb` by running the same command against this
branch.
