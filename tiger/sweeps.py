"""Parameter sweeps for the tail-length question.

* :func:`tail_length_sweep` runs the receptor model (optionally inside a flow:
  static cell, Stokes swimmer or TIGER LBM field) over tail lengths and
  replicates and summarises binding characteristics.
* :func:`encounter_vs_tail` gives the encounter kernel as a function of tail
  length: analytic diffusive/ballistic limits with capture radius
  Rc = Rb + Rp + L, plus the hydro simulation with the same Rc.
"""
from __future__ import annotations

from dataclasses import replace

import numpy as np
import pandas as pd

from .model.simulation import ReceptorSimParams, run_receptor_simulation


def summarize(res) -> dict:
    """Binding characteristics of one receptor-model run."""
    P = res.params["P"]
    T = res.params["nsteps"] * res.params["dt"]
    ta, tb = res.t_attach, res.t_bind
    att, bnd = np.isfinite(ta), np.isfinite(tb)
    both = att & bnd
    # Kaplan-Meier-free simple rates: events per phage-second at risk (censored at T)
    exp_att = np.where(att, ta, T).sum()
    exp_bnd = np.where(bnd, tb, T).sum()
    return dict(
        frac_attached=att.mean(), frac_bound=bnd.mean(),
        frac_bound_given_attached=(both.sum() / att.sum()) if att.any() else np.nan,
        median_t_attach=np.nanmedian(ta) if att.any() else np.nan,
        median_t_bind=np.nanmedian(tb) if bnd.any() else np.nan,
        median_search_time=np.nanmedian((tb - ta)[both & (tb >= ta)]) if both.any() else np.nan,
        rate_attach=att.sum() / exp_att, rate_bind=bnd.sum() / exp_bnd,
        mean_kernel=np.nanmean(res.kernel), n_phage=P, T=T,
    )


def survival_curve(times: np.ndarray, t_grid: np.ndarray) -> np.ndarray:
    """Fraction of phage that have NOT yet had the event by each time in t_grid."""
    t = np.where(np.isfinite(times), times, np.inf)
    return (t[None, :] > t_grid[:, None]).mean(axis=1)


def tail_length_sweep(tails, reps: int = 3, base: ReceptorSimParams | None = None, flow_factory=None,
                      seed: int = 0, progress: bool = True, keep_runs: bool = False):
    """Run the receptor model for each tail length (um) and replicate.

    ``flow_factory`` is None (static cell) or a callable ``tail -> flow`` so the
    flow can depend on the case. Tail mutation is switched off so each run keeps
    its tail length. Returns (summary DataFrame, survival DataFrame[, runs]).
    """
    base = base or ReceptorSimParams()
    base = replace(base, mutate_tails=False)
    rows, surv, runs = [], [], []
    t_grid = base.dt * np.arange(1, base.nsteps + 1)
    for i, L in enumerate(tails):
        for r in range(reps):
            p = replace(base, Ltail=float(L), seed=seed + 1000 * i + r)
            flow = flow_factory(L) if flow_factory else None
            res = run_receptor_simulation(p, flow=flow, warn_overlap=False)
            s = summarize(res); s.update(tail=L, rep=r)
            rows.append(s)
            surv.append(pd.DataFrame(dict(tail=L, rep=r, t=t_grid,
                                          S_attach=survival_curve(res.t_attach, t_grid),
                                          S_bind=survival_curve(res.t_bind, t_grid))))
            if keep_runs:
                runs.append(res)
            if progress:
                print(f"tail {L:.3f} um rep {r}: bound {s['frac_bound']:.2f}, attached {s['frac_attached']:.2f}")
    out = (pd.DataFrame(rows), pd.concat(surv, ignore_index=True))
    return out + (runs,) if keep_runs else out


def encounter_vs_tail(tails, Rb=1.0, Rp=0.05, D0=None, v=30.0, simulate=True, P=20_000, nsteps=300,
                      seed=0, **hydro_kw):
    """Encounter kernel vs tail length with capture radius Rc = Rb + Rp + L (um^3/s).

    Analytic: k_diff = 4 pi D0 Rc, k_ball = pi Rc^2 v. Simulation:
    tiger.encounter.encounter_rates_hydro with Rcrit = Rc (exclusion radius Rb so
    diffusion can reach the capture shell).
    """
    from .encounter.hydro import encounter_rates_hydro
    from .model.surface_walk import stokes_einstein_D

    D0 = stokes_einstein_D(Rp) if D0 is None else D0
    rows = []
    for i, L in enumerate(tails):
        Rc = Rb + Rp + L
        row = dict(tail=L, Rc=Rc, k_diff=4 * np.pi * D0 * Rc, k_ball=np.pi * Rc ** 2 * v)
        if simulate:
            h = encounter_rates_hydro(P=P, nsteps=nsteps, v=v, Rb=Rb, Rp=Rp + L, D0=D0,
                                      exclusion_radius=Rb, seed=seed + i, **hydro_kw)
            row["k_sim"] = np.nanmean(h["kernel"][nsteps // 5:])
        rows.append(row)
    return pd.DataFrame(rows)
