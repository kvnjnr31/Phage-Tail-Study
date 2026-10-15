"""Receptor-level phage search simulation (driver).

Port of encounter_rates_main.m (Phage_Bacteria_2/Simulation): phage start in
a thin shell around a capsule-shaped E. coli covered by a helical strip of
receptors, jump through the bulk (CTRW), attach on contact, random-walk on
the surface, and bind receptors that fall inside their tail cone.

Each step, in the MATLAB order:
  surface_step  -> overlap check
  [flow advection of free phage]      (new: optional, see tiger.flow)
  bulk_jumps    -> overlap check
  binding_fun   -> metrics
  boundarys_receptor -> overlap check

Additions over the MATLAB driver (all optional, defaults reproduce it):
* ``flow``: an object with velocity(pos, t) (um/s) that advects detached phage,
  e.g. the Stokes flow of a swimming cell or the TIGER LBM field;
* per-phage event times: first attachment and first receptor binding;
* ``record_every``: snapshots for rendering (tiger.viz.render_tiger_movie);
* ``mutate_tails``: the MATLAB boundary routine mutates tail length at the
  walls; switch off to hold the tail length fixed (needed for tail sweeps).
"""
from __future__ import annotations

from dataclasses import dataclass, field, asdict

import numpy as np

from .binding import binding_fun, boundarys_receptor, overlap
from .geometry import distance_to_axis, axis_point
from .phage import PhageState, init_phage
from .receptors import init_receptors
from .surface_walk import bulk_jumps, stokes_einstein_D, surface_step


@dataclass
class ReceptorSimParams:
    # time
    nsteps: int = 2400
    dt: float = 1e-2
    # box (um)
    LX: float = 15.0
    # phage
    P: int = 100
    Rp: float = 0.1          # capsid radius (um)
    Ltail: float = 0.05      # tail length (um)
    T: float = 300.0
    eta: float = 1e-3        # water viscosity (Pa s)
    nu: float = 1.0          # CTRW exponent (0.82 in some runs)
    tmin: float = 0.001      # CTRW minimum waiting time
    shell: float = 0.05      # initial shell thickness around the cylinder (init_phage.m)
    # receptors
    B: int = 1000
    m: float = 0.05 * np.pi  # helical pitch angle
    delta: float = 0.05      # receptor layer half-thickness
    Rr: float = 1e-4         # receptor capture radius
    k_bind: float = 1.0
    # bacterium (um)
    Rb: float = 0.8
    bac_length: float = 1.9  # half-length of cylindrical body
    # switches
    tether_force_scale: float = 0.0   # MATLAB: 0*F_proj
    binding_rule: str = "matlab"      # or "slides"
    mutate_tails: bool = True         # MATLAB boundarys_receptor mutates tails at the walls
    seed: int | None = None

    @property
    def D0(self) -> float:
        return stokes_einstein_D(self.Rp, self.T, self.eta)


@dataclass
class ReceptorSimResult:
    params: dict
    receptors: np.ndarray
    phage: PhageState
    kernel: np.ndarray            # (Nenc/Nact) * V/dt per step
    n_attached: np.ndarray        # MATLAB "mean_head_data" = sum(MatP(:,9))
    head_std: np.ndarray
    tail_median: np.ndarray
    tail_std: np.ndarray
    walk: np.ndarray              # trajectory of phage 0 (nsteps x 3)
    n_bound: np.ndarray
    t_attach: np.ndarray = None   # first attachment time per phage (NaN = never)
    t_bind: np.ndarray = None     # first receptor binding time per phage (NaN = never)
    frames: list = field(default_factory=list)

    @property
    def time(self) -> np.ndarray:
        return self.params["dt"] * np.arange(1, len(self.kernel) + 1)


def _push_out(ph: PhageState, Rb: float, L: float):
    """Move free phage that flow carried inside the exclusion shell back onto it."""
    free = ~ph.attached
    Rmin = Rb + ph.head_radius
    d = distance_to_axis(ph.pos, L)
    bad = np.flatnonzero(free & (d < Rmin))
    if bad.size:
        c = axis_point(ph.pos[bad], L)
        v = ph.pos[bad] - c
        n = np.linalg.norm(v, axis=1, keepdims=True)
        n[n == 0] = 1.0
        ph.pos[bad] = c + (Rmin[bad] * 1.0001)[:, None] * v / n
    return bad


def run_receptor_simulation(params: ReceptorSimParams | None = None, progress: bool = False,
                            flow=None, record_every: int = 0, warn_overlap: bool = True) -> ReceptorSimResult:
    p = params or ReceptorSimParams()
    rng = np.random.default_rng(p.seed)
    half = p.LX / 2
    V = p.LX ** 3
    D0 = p.D0

    MatB = init_receptors(p.B, p.delta, p.m, p.Rb, p.bac_length, p.Rr, rng)
    ph = init_phage(p.P, p.Rb, p.Rp, p.LX, p.Ltail, shell=p.shell, rng=rng)

    kernel = np.full(p.nsteps, np.nan)
    n_att = np.zeros(p.nsteps)
    head_std = np.zeros(p.nsteps)
    tail_med = np.zeros(p.nsteps)
    tail_std = np.zeros(p.nsteps)
    n_bound = np.zeros(p.nsteps)
    walk = np.zeros((p.nsteps, 3))
    t_attach = np.full(p.P, np.nan)
    t_bind = np.full(p.P, np.nan)
    frames = []

    time = 0.0
    for n in range(p.nsteps):
        if progress and n % max(p.nsteps // 10, 1) == 0:
            print(f"{100 * n / p.nsteps:5.1f}%")
        surface_step(ph, MatB, p.Rb, D0, p.k_bind, p.bac_length, p.dt, rng,
                     force_scale=p.tether_force_scale)
        overlap(ph, p.Rb, p.bac_length, warn=warn_overlap)
        if flow is not None:
            free = ~ph.attached
            if free.any():
                ph.pos[free] += flow.velocity(ph.pos[free], time) * p.dt
                _push_out(ph, p.Rb, p.bac_length)
        bulk_jumps(ph, D0, p.dt, p.tmin, p.nu, time, p.Rb, p.bac_length, rng)
        overlap(ph, p.Rb, p.bac_length, warn=warn_overlap)
        Nenc, Nact, ph = binding_fun(ph, MatB, p.Rr, rng, rule=p.binding_rule)

        now = time + p.dt
        t_attach[np.isnan(t_attach) & ph.attached] = now
        t_bind[np.isnan(t_bind) & ph.bound] = now

        kernel[n] = (Nenc / Nact) * (V / p.dt) if Nact else np.nan
        n_att[n] = ph.attached.sum()
        head_std[n] = ph.head_radius.std()
        tail_med[n] = np.median(ph.tail_length)
        tail_std[n] = ph.tail_length.std()
        n_bound[n] = ph.bound.sum()
        walk[n] = ph.pos[0]

        if p.mutate_tails:
            boundarys_receptor(ph, half, half, half, p.Rp, rng)
        else:
            boundarys_receptor(ph, half, half, half, p.Rp, rng, tail_mutation=0.0)
        overlap(ph, p.Rb, p.bac_length, warn=warn_overlap)
        if record_every and n % record_every == 0:
            frames.append(dict(t=now, pos=ph.pos.copy(), attached=ph.attached.copy(), bound=ph.bound.copy()))
        time = now

    return ReceptorSimResult(asdict(p), MatB, ph, kernel, n_att, head_std, tail_med,
                             tail_std, walk, n_bound, t_attach, t_bind, frames)
