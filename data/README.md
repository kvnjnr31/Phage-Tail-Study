# data — TIGER working data

| file | contents |
|---|---|
| `adsorption_data.xlsx` | **The working adsorption dataset.** It follows the 2026 template and has four sheets. **adsorption** has 66 rates in 8 columns: Reference, Host, Tail type, Phage, Study, Bac type, Adsorption Rate d(ml/h), Temp deg C. **codes** explains the Study and Bac type codes. **data_issues** lists rows to check against the papers. **changes** records what was edited. |
| `adsorption_references.xlsx` | One row per reference: year, title, URL, whether the dataset cites it, and notes. |
| `phage_characteristics.xlsx` | **phage** has one row per phage in the adsorption data: family, hosts, capsid and tail dimensions (nm), the basis of each value (measured, assumed, or family typical) and its source. **literature_values** lists every measurement on file with its reference. **family_typical** holds the per-family pools of tail and capsid values, used for phage without their own values. |
| `bacteria_characteristics.xlsx` | **hosts** has one row per host: diameter, length, equivalent spherical diameter, swim speed and run time, with sources and checks. **ecoli_motility_vs_T** holds Maeda 1976 E. coli swim speed and tumble rate against temperature. |
| `binding_force_curve_gp37_LPS.csv` | Force–distance curve for the T4 gp37 adhesin against E. coli LPS (Brzozowska 2019), digitised. |
| `schwartz1976_receptor_density_vs_adsorption.csv` | Lambda adsorption rate against receptor density (Schwartz 1976), digitised. |

Load these from Python with the `tiger.analysis` loaders: `load_adsorption`, `load_phage_characteristics`,
`load_bacteria_characteristics`, `maeda_arrays` and `load_force_curve`.

The older versions are in `../archive/data_legacy`. That folder has the 2025 adsorption sheet, which includes
Schwartz, Chattopadhyay, Knezevic and Waterbury rows that the 2026 template dropped. It also has the bootstrap
sheet, the original CSVs, and the 2025 phage characteristics sheet. The unlabelled digitised curves `sipho` and
`podo_delta`, which nothing used, are there as CSV files.
