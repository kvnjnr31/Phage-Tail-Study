# 2025 violin-plot study (archived)

The original contents of this repository (Nov 2025): the four-panel violin figure of literature adsorption
rates, its plotting helpers, the Word report and the templates it was built from.

| file | what it is |
|---|---|
| `data_violinplots.py`, `violin.py`, `violinplot.py` | 2025 figure script and helpers. They read the raw 2026 spreadsheet (`archive/data_legacy/adsorption_data_2026.xlsx`) with hard-coded tail and head arrays. |
| `figure_violin_panels.png` | The figure they produced (manuscript Fig. 4). |
| `tails_heads_violin_template.py`, `figure_violin_tails_heads.png` | Template for panel A. |
| `phage_tail_vs_adsorption.docx`, `phage_tail_vs_adsorption_template.docx` | Report draft and its template, superseded by `manuscript/`. |
| `adsorption_data_template.xlsx` | Column template used to clean `data/adsorption_data.xlsx`. |

**Replaced by** `experiments/E4_literature_adsorption_vs_tail_length/violin_figure.py` and
`tiger/analysis/violin.py`, which draw the same figure from the cleaned tables in `data/`.

Two duplicates were dropped when the repo was consolidated: `deliverables/adsorption_data.xlsx` (identical to
`archive/data_legacy/adsorption_data_2026.xlsx`) and `references/adsorption_data_urls.csv` (identical to
`archive/data_legacy/adsorption_data_2026_urls.csv`).
