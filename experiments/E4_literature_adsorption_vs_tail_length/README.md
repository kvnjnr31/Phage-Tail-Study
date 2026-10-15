# E4 — Literature adsorption rates vs tail length and family

**Question.** What does the measured adsorption data say about tail length? The simulations (E1–E3) need to
explain it.

**Analyses** (all on the cleaned tables in `data/`: `adsorption_data.xlsx`, `phage_characteristics.xlsx`,
`bacteria_characteristics.xlsx`).
1. **Bootstrap regression**: z-scored rate ~ bacterial size + tail + capsid over Study-group means, resampling the
   measurements and drawing family-typical sizes for phage without their own values. Port of
   `adsorption_model_bootstrap.m`.
2. **Family and host comparisons** (rank-sum tests, boxplots).
3. **Kinetic predictors** per group: ballistic and run-and-tumble diffusive limits vs the data.
4. **Theory.** A map of k_diff over capsid diameter × tail length, and Moldovan's lambda temperature series
   regressed on k_diff.

**Run.** `python run.py --trials 2000`

**Violin figure.** `python violin_figure.py` draws the four-panel figure (tail/capsid sizes, rates by family,
family × host, host) from the cleaned data and writes Mann–Whitney tests to `results/violin_tests.csv`. It ports
the 2025 `data_violinplots.py` (kept in `archive/violin_study_2025/`); the plotting helpers are `tiger.analysis.violin`.

**Results (Oct 2026).** Myoviridae vs Siphoviridae p = 0.018 and Myoviridae vs Podoviridae p = 0.011 (Mann–Whitney);
Podoviridae vs Siphoviridae is not significant. The Abedon 2023 T4 row was relabelled Myoviridae during cleaning
(see the `changes` sheet in `data/adsorption_data.xlsx`).

**Legacy.** `legacy_adsorption_studies/adsorption_stats_from_csv.py` is the 2025 boxplot script, which reads
`adsorption_data_original.csv`, now in `archive/data_legacy/`.
