"""Adsorption-rate data: loading, family comparisons, tail-length regression.

Reads the cleaned files in ``code/data``:

* ``adsorption_data.xlsx``          sheet "adsorption" (one row per measured rate)
* ``phage_characteristics.xlsx``    sheets "phage" (one row per phage), "family_typical"
* ``bacteria_characteristics.xlsx`` sheets "hosts" (one row per host), "ecoli_motility_vs_T"
* ``binding_force_curve_gp37_LPS.csv``

Ports Phage_Bacteria_2/Data: adsorption_model_bootstrap.m (tail-length bootstrap),
data_boxplots.m (rank-sum tests), model_scratch.m / model_fun.m (ballistic vs
diffusive predictors), AdsorptionRatesViolinPlots.m and the T2/T4/T6 block of
LinearModelingMoldovan.m.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .stats import ols

DATA_DIR = Path(__file__).resolve().parents[2] / "data"
HOST_OF_CODE = {1: "Enterobacteriaceae", 2: "Synechococcaceae", 3: "Pseudomonadaceae", 4: "Vibrionaceae"}
FAMILIES = ("Siphoviridae", "Myoviridae", "Podoviridae")


def _p(path, name):
    return Path(path) if path is not None else DATA_DIR / name


# ----------------------------------------------------------------------------- loaders
def load_adsorption(path=None) -> pd.DataFrame:
    """Tidy adsorption table: reference, host, family, phage, study, host_code, rate (mL/h), temp_C."""
    raw = pd.read_excel(_p(path, "adsorption_data.xlsx"), sheet_name="adsorption")
    df = pd.DataFrame({
        "reference": raw["Reference"].astype(str).str.strip(), "host": raw["Host"].astype(str).str.strip(),
        "family": raw["Tail type"].astype(str).str.strip(), "phage": raw["Phage"].astype(str).str.strip(),
        "study": pd.to_numeric(raw["Study"], errors="coerce"),
        "host_code": pd.to_numeric(raw["Bac type"], errors="coerce"),
        "rate": pd.to_numeric(raw["Adsorption Rate d(ml/h)"], errors="coerce"),
        "temp_C": pd.to_numeric(raw["Temp deg C"], errors="coerce"),
    })
    return df[df.rate.notna()].reset_index(drop=True)


def load_phage_characteristics(path=None):
    """Return (phage table, family pools dict family -> {"tail": array, "capsid": array}) in nm."""
    p = _p(path, "phage_characteristics.xlsx")
    ph = pd.read_excel(p, sheet_name="phage")
    fam = pd.read_excel(p, sheet_name="family_typical")
    pools = {}
    for _, r in fam.iterrows():
        vals = np.array([float(v) for v in str(r["values_nm"]).split(",")])
        pools.setdefault(r["family"], {})["tail" if "tail" in r["quantity"] else "capsid"] = vals
    return ph, pools


def load_bacteria_characteristics(path=None):
    """Return (hosts table, Maeda E. coli table with temperature, swim speed, tumble rate)."""
    p = _p(path, "bacteria_characteristics.xlsx")
    return pd.read_excel(p, sheet_name="hosts"), pd.read_excel(p, sheet_name="ecoli_motility_vs_T")


def maeda_arrays(path=None):
    """(speed, tumble) as (T_C, value) arrays, the layout of the old MaedaSpeed/MaedaTumble.mat."""
    _, m = load_bacteria_characteristics(path)
    return (m[["Temp deg C (speed)", "Swim speed (um/s)"]].to_numpy(float),
            m[["Temp deg C (tumble)", "Tumble rate (1/s)"]].to_numpy(float))


def load_force_curve(path=None) -> np.ndarray:
    """gp37-LPS force curve as an (N, 2) array: distance (nm), force (pN)."""
    return pd.read_csv(_p(path, "binding_force_curve_gp37_LPS.csv")).to_numpy(float)


def host_table(hosts: pd.DataFrame) -> pd.DataFrame:
    """One row per Bac type code (first listed representative species)."""
    return hosts.drop_duplicates("Bac type").set_index("Bac type")


# ----------------------------------------------------------------------------- family comparisons
def family_ranksum(df: pd.DataFrame):
    """Wilcoxon rank-sum tests between tail families (data_boxplots.m)."""
    from scipy.stats import ranksums

    g = {f: df.loc[df.family == f, "rate"].to_numpy() for f in FAMILIES}
    out = {}
    for a, b in (("Siphoviridae", "Myoviridae"), ("Siphoviridae", "Podoviridae"), ("Myoviridae", "Podoviridae")):
        r = ranksums(g[a], g[b])
        out[f"{a} vs {b}"] = dict(statistic=r.statistic, p=r.pvalue, n=(len(g[a]), len(g[b])))
    return out


def plot_family_boxplots(df: pd.DataFrame, ax=None):
    """Grouped log-scale boxplots by family and host (data_boxplots.m)."""
    import matplotlib.pyplot as plt

    ax = ax or plt.gca()
    colors = {"Podoviridae": "tab:red", "Myoviridae": "tab:green", "Siphoviridae": "tab:blue"}
    data, pos, labels, cols = [], [], [], []
    for i, fam in enumerate(("Podoviridae", "Myoviridae", "Siphoviridae")):
        hosts = sorted(df.loc[df.family == fam, "host_code"].dropna().unique())
        for j, h in enumerate(hosts):
            v = df.loc[(df.family == fam) & (df.host_code == h), "rate"].to_numpy()
            data.append(v); pos.append(i + 1 + 0.22 * j); labels.append(HOST_OF_CODE.get(int(h), str(h)))
            cols.append(colors[fam])
    bp = ax.boxplot(data, positions=pos, widths=0.18, patch_artist=True)
    for patch, c in zip(bp["boxes"], cols):
        patch.set_facecolor(c); patch.set_alpha(0.5)
    ax.set_xticks(pos); ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_yscale("log"); ax.set_ylabel("adsorption rate (mL/h)")
    for fam, c in colors.items():
        ax.plot([], [], "s", color=c, alpha=0.5, label=fam)
    ax.legend(fontsize=8)
    return ax


# ----------------------------------------------------------------------------- tail-length regression
def _row_sizes(ads, phage, pools, rng):
    """Tail and capsid (nm) for each adsorption row: the phage's own value when it is
    measured/assumed from a named phage, otherwise a draw from its family pool
    (as adsorption_model_bootstrap.m did for unnamed isolates)."""
    ph = phage.set_index("Phage")
    tails = np.empty(len(ads)); caps = np.empty(len(ads))
    for i, (p, fam) in enumerate(zip(ads.phage, ads.family)):
        own = p in ph.index and not str(ph.loc[p, "Value basis"]).startswith("family typical")
        if own:
            tails[i] = ph.loc[p, "Tail length (nm)"]; caps[i] = ph.loc[p, "Capsid diameter (nm)"]
        else:
            tails[i] = rng.choice(pools[fam]["tail"]); caps[i] = rng.choice(pools[fam]["capsid"])
    return tails, caps


def tail_regression_bootstrap(ads=None, phage=None, pools=None, hosts=None, n_trials=2000, seed=None,
                              group="study"):
    """Bootstrap of the z-scored regression  A ~ B + T + H  over group means.

    Each trial resamples the measurements within each group (default: Study code,
    i.e. family x host), assigns tail length T and capsid H per row (own value or
    family-pool draw), takes the host's equivalent spherical diameter B, averages
    per group, z-scores and fits OLS. Returns a DataFrame of R2, coefficients and
    p-values per trial. Generalises adsorption_model_bootstrap.m to the cleaned data.
    """
    if ads is None:
        ads = load_adsorption()
    if phage is None or pools is None:
        phage, pools = load_phage_characteristics()
    if hosts is None:
        hosts, _ = load_bacteria_characteristics()
    H = host_table(hosts)
    B_of = {k: H.loc[k, "Equivalent spherical diameter (um)"] for k in H.index}
    rng = np.random.default_rng(seed)
    groups = [g for _, g in ads.groupby(group)]
    rows = []
    for _ in range(n_trials):
        data = []
        for g in groups:
            s = g.iloc[rng.integers(0, len(g), len(g))]
            t, c = _row_sizes(s, phage, pools, rng)
            data.append([np.mean([B_of[k] for k in s.host_code]), t.mean(), c.mean(), s.rate.mean()])
        X = np.array(data)
        sd = X.std(0, ddof=1); sd[sd == 0] = 1
        Z = (X - X.mean(0)) / sd
        fit = ols(Z[:, 3], Z[:, :3], names=["B", "T", "H"])
        rows.append(dict(R2=fit.r2, b_B=fit.beta[1], b_T=fit.beta[2], b_H=fit.beta[3],
                         p_B=fit.pvalues[1], p_T=fit.pvalues[2], p_H=fit.pvalues[3]))
    return pd.DataFrame(rows)


def kernel_predictors(ads=None, phage=None, pools=None, hosts=None, seed=0):
    """model_scratch.m: per Study group, J_ball = pi Rc^2 v and J_diff = 4 pi Rc D_B (mL/h).

    Rc = tail + capsid/2 + ESD/2 (median tail/capsid of the group's rows); v and run
    time from the host table; D_B = v^2 tau / 3. J_diff falls back to J_ball when the
    run time is missing or infinite.
    """
    if ads is None:
        ads = load_adsorption()
    if phage is None or pools is None:
        phage, pools = load_phage_characteristics()
    if hosts is None:
        hosts, _ = load_bacteria_characteristics()
    H = host_table(hosts)
    rng = np.random.default_rng(seed)
    out = []
    for g, d in ads.groupby("study"):
        t, c = _row_sizes(d, phage, pools, rng)
        h = H.loc[int(d.host_code.iloc[0])]
        Rc = (np.median(t) + np.median(c) / 2) * 1e-9 + h["Equivalent spherical diameter (um)"] / 2 * 1e-6
        v = h["Swim speed (um/s)"] * 1e-6
        tau = h["Run time (s)"]
        Jb = np.pi * Rc ** 2 * v * 1e6 * 3600
        D = v * v * tau / 3
        Jd = Jb if not np.isfinite(D) else 4 * np.pi * Rc * D * 1e6 * 3600
        out.append(dict(study=g, family=d.family.iloc[0], host=d.host.iloc[0], Rc_um=Rc * 1e6,
                        J_ball=Jb, J_diff=Jd, rate_mean=d.rate.mean(), rate_median=d.rate.median(), n=len(d)))
    return pd.DataFrame(out)


def host_rates_log10(ads=None):
    """AdsorptionRatesViolinPlots.m: log10 adsorption rates by host family."""
    ads = load_adsorption() if ads is None else ads
    return {h: np.log10(d.rate.to_numpy()) for h, d in ads.groupby("host")}


def t_even_rates(ads=None):
    """T2/T4/T6 adsorption rates on E. coli (LinearModelingMoldovan.m block)."""
    ads = load_adsorption() if ads is None else ads
    return {p: ads.loc[ads.phage == p, "rate"].to_numpy() for p in ("T4", "T2", "T6")}
