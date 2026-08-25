# Legacy v1 (notebook implementation)

Historical reference only — not maintained, not used by the running app, not part of the
Docker build (`.dockerignore` excludes this whole directory).

- `data_transformation.ipynb`, `data_analysis.ipynb`, `heatmap_visualization.ipynb`,
  `heatmap_compact.ipynb`, `heatmap_visualization test.ipynb` — the original Jupyter/
  ipywidgets prototype that the current Django app (`data_transformation/`, `data_analysis/`,
  `heatmap_viz/`) was extracted and rewritten from. Several source files in those modules
  still reference these notebooks in comments (e.g. "Extracted from data_analysis.ipynb
  Plotter class") for historical traceability.
- `_settings.py` — rendering constants (colormap lists, species-translation list) the
  notebooks imported; `heatmap_viz/services/heatmap_renderer.py` now carries its own copy.
- `notebook_requirements.txt` — pip dependencies for running the notebooks locally.
- `supplementals/` — an older, superseded set of manuscript supporting figures/tables
  (S1-S13 tables, S1-S6 figures, `zukaitis_2026`-style numbering). Confirmed against the live
  hosted app (2026-08-25) that production actually serves a different, current set (S1-S15
  tables, S1-S5 figures) from `static/peptide/publications/kuhfeld_2026/supplementals/` —
  see `docs/SPLIT_PLAN.md` section 0c for how that was discovered. Kept here for reference
  only; do not link to it from the running app.
