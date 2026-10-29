# Nucleic-acid campaign

39 structures with DNA or RNA (`scripts/na_validate.py`), each drawn and judged against its 3D structure and, with
`--require-rnapolis`, against RNApolis's base pairs. Limits: `docs/validation/na/report.md`.

- `report.md`, `report.json`, `<ID>.png` (contact sheets) and `<ID>-figure.svg`: the run at the paper's release,
  v0.2.0 (code fingerprint 7879eea05bd9301c), made by `M2_RERUN=1 scripts/paper/run_all.sh v0.2.0 2f53b93` on
  zyla-office: 34 of 39 pass, the count of plan 15's run (1d38e45). Against that run, every entry gains
  `stems_nearer_not_applicable` (NA fix F7), and 30QM, still failing `stems_nearer_frac`, no longer fails spearman,
  procrustes and the angles (NA fix F2, each nucleic acid unrolled with the protomer holding it;
  `docs/validation/na-fixes/`). The figures are byte for byte those committed before.
- `report-step*.json/md`, `run*.log`: the campaign's earlier steps (Task 13 and plan 15), kept as they are.
- `review.md`: the campaign's phase A review (26 of 39 then), kept as it is.
