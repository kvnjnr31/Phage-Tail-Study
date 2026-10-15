"""E4 - Literature adsorption rates vs tail length and family (the data the model must explain).

Uses the cleaned data in ../../data:
  adsorption_data.xlsx, phage_characteristics.xlsx, bacteria_characteristics.xlsx

1. Tail-length bootstrap regression  A ~ bacterial size (B) + tail (T) + capsid (H),
   z-scored over Study-group means, resampling measurements and drawing family-typical
   sizes for phage without their own values.
2. Family / host comparisons (rank-sum tests, boxplots).
3. Kinetic predictors per group: ballistic and run-and-tumble diffusive limits vs data.
4. Theory: k_diff over capsid diameter x tail length; Moldovan lambda series vs k_diff.

    python run.py --trials 2000

Outputs (results/): bootstrap.csv, family_tests.csv, kernel_predictors.csv, *.png
"""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from tiger.analysis import (family_ranksum, kdiff_size_map, kernel_predictors, linear_model_moldovan,
                            load_adsorption, load_bacteria_characteristics, load_phage_characteristics,
                            maeda_arrays, plot_family_boxplots, tail_regression_bootstrap)

DATA = pathlib.Path(__file__).resolve().parents[2] / "data"
ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--trials", type=int, default=2000)
ap.add_argument("--data", default=str(DATA))
ap.add_argument("--out", default=str(pathlib.Path(__file__).parent / "results"))
a = ap.parse_args()
D = pathlib.Path(a.data)
out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True)

ads = load_adsorption(D / "adsorption_data.xlsx")
phage, pools = load_phage_characteristics(D / "phage_characteristics.xlsx")
hosts, _ = load_bacteria_characteristics(D / "bacteria_characteristics.xlsx")
print(f"{len(ads)} adsorption measurements, {ads.phage.nunique()} phage, {ads.reference.nunique()} references")

# 1. bootstrap on tail length
bt = tail_regression_bootstrap(ads, phage, pools, hosts, n_trials=a.trials, seed=0)
bt.to_csv(out / "bootstrap.csv", index=False)
fig, axs = plt.subplots(1, 4, figsize=(16, 3.8))
for ax, col in zip(axs, ["R2", "b_B", "b_T", "b_H"]):
    ax.hist(bt[col], 25, color="0.6", ec="k"); ax.set_title(col)
fig.suptitle("Bootstrap: z-scored adsorption ~ bacterium (B) + tail (T) + capsid (H), Study-group means")
fig.tight_layout(); fig.savefig(out / "bootstrap_regression.png", dpi=150)
print("bootstrap medians:\n" + bt.median().round(3).to_string())

# 2. family comparisons
rows = [dict(comparison=k, statistic=v["statistic"], p=v["p"], n=v["n"]) for k, v in family_ranksum(ads).items()]
pd.DataFrame(rows).to_csv(out / "family_tests.csv", index=False)
print(pd.DataFrame(rows).round(4).to_string())
fig, ax = plt.subplots(figsize=(8, 5)); plot_family_boxplots(ads, ax)
fig.tight_layout(); fig.savefig(out / "families.png", dpi=150)

# 3. kinetic predictors per group
kp = kernel_predictors(ads, phage, pools, hosts)
kp.to_csv(out / "kernel_predictors.csv", index=False)
fig, ax = plt.subplots(figsize=(6, 5))
ax.loglog(kp.J_ball, kp.rate_median, "o", label="ballistic pi Rc^2 v")
ax.loglog(kp.J_diff, kp.rate_median, "s", label="run-tumble 4 pi Rc D_B")
lim = [min(kp.rate_median.min(), kp.J_ball.min()) / 2, max(kp.J_diff.max(), kp.rate_median.max()) * 2]
ax.plot(lim, lim, "k--", lw=1)
ax.set(xlabel="predicted (mL/h)", ylabel="measured median (mL/h)"); ax.legend()
fig.tight_layout(); fig.savefig(out / "kernel_predictors.png", dpi=150)

# 4. theory
X, Y, K = kdiff_size_map()
fig, ax = plt.subplots(figsize=(6, 5))
c = ax.pcolormesh(X, Y, K, norm=matplotlib.colors.LogNorm(), shading="auto")
fig.colorbar(c, label="k_diff (mL/hr)"); ax.set(xlabel="capsid diameter (nm)", ylabel="tail length (nm)")
fig.tight_layout(); fig.savefig(out / "kdiff_capsid_tail_map.png", dpi=150)
speed, _ = maeda_arrays(D / "bacteria_characteristics.xlsx")
for w in ("maltose", "glucose"):
    print(f"Moldovan {w}: k_exp ~ k_diff\n" + linear_model_moldovan(speed, w)["fit"].summary())
print(f"results in {out}")
