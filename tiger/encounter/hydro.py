"""Encounter rates with hydrodynamics: phage near a swimming bacterium in its own frame.

Ports encounter_rates_hydro.m and its helpers. The same helper files appear in
the MATLAB surface_walks folder and in Phage_Bacteria_1/Encounter Rates/
encounter_rates_hydro; they are identical (originals in
archive/matlab_reference/encounter_rates_hydro/).

  MATLAB                    Python
  Langevin_fun.m            hydro_advect
  noise_gen.m               hydro_noise   (noise_gen.asv is an older autosave of it)
  encounter_fun_hydro.m     hydro_encounter (+ deposited fraction N_dep in the driver)
  boundarys_hydro.m         hydro_boundaries
  overlap.m                 hydro_overlap (+ plot_overlap for the debug figure)
  sim_plotter.m             plot_concentration, rotate_to_xy_plane, plot_cloud

The bacterium (sphere radius Rb) sits at the origin swimming with velocity
U = (-v, 0, 0). In its frame each phage is advected by the Stokes flow
U - u_s(r), plus a rotational term (a rotlet about the swimming axis)
    dy = -U_x * phi_3 dt,  dz = +U_x * phi_2 dt,  phi_i = -Rb^3 r_i / r^3,
i.e. u_rot = Rb^3 (omega x r) / r^3 with |omega| = |U_x| / (1 um). The MATLAB
uses the swim speed as the angular speed, so the swirl at the cell surface equals
the swim speed; ``swirl`` scales it (0 switches it off).
then makes a CTRW jump (redrawn while it would land inside Rb + Rp).
Contacts within Rcrit = Rb + Rp are counted, reflected specularly off the
capture sphere (event-driven correction), and the phage is deactivated until it
leaves through the upstream boundary.

``dims`` = (D_para_y, D_para_z) switches the y and z directions on/off (1/0) to
run 1D/2D/3D variants, as in the MATLAB.
"""
from __future__ import annotations

import numpy as np

from .noise import pareto_wait, uniform_kick
from .stokes import stokes_disturbance


def hydro_advect(pos, U, Rb, dt, mask, swirl=1.0):
    """Langevin_fun.m (hydro): Stokes advection + rotational flow, returns new positions.

    ``swirl`` multiplies the rotational term (1 = MATLAB, 0 = translating sphere only).
    """
    r = pos * mask
    us = stokes_disturbance(r, U, Rb) * mask
    rr = np.linalg.norm(r, axis=1)
    rr[rr == 0] = np.inf
    phi2 = -Rb ** 3 * r[:, 1] / rr ** 3
    phi3 = -Rb ** 3 * r[:, 2] / rr ** 3
    new = pos - us * dt + U * dt
    new[:, 1] += -swirl * U[0] * phi3 * dt
    new[:, 2] += swirl * U[0] * phi2 * dt
    return new * mask


def hydro_noise(pos, t_wait, time, D0, dt, tmin, nu, Rexcl, mask, rng):
    """noise_gen.m (hydro): CTRW jumps for phage whose waiting time elapsed.

    Fix: the MATLAB loop updated rows ``p = 1..len(kk)`` instead of ``kk(p)``,
    so the wrong phage jumped; here the ready phage jump. Jumps landing within
    ``Rexcl`` of the centre are redrawn.
    """
    kk = np.flatnonzero(time - t_wait >= 0)
    t_wait[kk] = time + pareto_wait(rng, kk.size, tmin, nu)
    pending = kk.copy()
    while pending.size:
        cand = (pos[pending] + uniform_kick(rng, (pending.size, 3), D0, dt)) * mask
        ok = np.linalg.norm(cand, axis=1) >= Rexcl
        pos[pending[ok]] = cand[ok]
        pending = pending[~ok]
    return pos, t_wait


def hydro_encounter(pos, prev, active, res_time, Rcrit, Rp, dt):
    """encounter_fun_hydro.m: count contacts, track residence time, reflect at Rcrit.

    Returns (Nenc, Nact, deposition_velocities, residence_times_of_hits).
    """
    Nact = int(active.sum())
    dist = np.linalg.norm(pos, axis=1)
    near = dist <= Rcrit + Rp
    res_time[near] += dt * active[near]
    hit = np.flatnonzero(dist <= Rcrit)
    Nenc = int(active[hit].sum())
    dep_v = np.empty(0)
    res_hits = res_time[hit[active[hit]]].copy()
    if hit.size:
        a = hit[active[hit]]
        v = (pos[a] - prev[a]) / dt
        pn = np.linalg.norm(prev[a], axis=1)
        pn[pn == 0] = 1
        dep_v = -np.sum(v * prev[a], 1) / pn

        # event-driven reflection off the capture sphere
        d = pos[hit] - prev[hit]
        r0 = prev[hit]
        vv = np.sum(d * d, 1)
        vv[vv == 0] = np.inf
        b = 2 * np.sum(r0 * d, 1) / vv
        c = (np.sum(r0 * r0, 1) - Rcrit ** 2) / vv
        disc = b * b - 4 * c
        ok = disc >= 0
        tau = np.where(ok, -(b + np.sqrt(np.where(ok, disc, 0))) / 2, 0.0)
        rc = d * tau[:, None] + r0
        q = pos[hit] - rc
        rcn = np.sum(rc * rc, 1)
        rcn[rcn == 0] = 1
        perp = -rc * (np.sum(rc * q, 1) / rcn)[:, None]
        newp = rc + q + 2 * perp
        pos[hit[ok]] = newp[ok]
    active[hit] = False
    res_time[hit] = 0.0
    res_time[dist > Rcrit + Rp] = 0.0
    return Nenc, Nact, dep_v, res_hits


def hydro_boundaries(pos, prev, active, half, Rp, mask, rng):
    """boundarys_hydro.m: upstream re-injection in x, periodic wrap in y and z."""
    jj = np.flatnonzero(np.abs(pos[:, 0]) >= half)
    pos[jj, 0] = half - Rp
    pos[jj, 1] = mask[1] * (2 * half * rng.random(jj.size) - half)
    pos[jj, 2] = mask[2] * (2 * half * rng.random(jj.size) - half)
    prev[jj] = pos[jj]; active[jj] = True
    for ax in (1, 2):
        jj = np.flatnonzero(np.abs(pos[:, ax]) >= half)
        pos[jj, ax] = mask[ax] * (pos[jj, ax] - 2 * half * np.sign(pos[jj, ax]))
        prev[jj] = pos[jj]; active[jj] = True
    return pos


def _tumble_rotation(pos, rng, angle_deg, dims):
    """Random small reorientation of the frame (the disabled tumble block of encounter_rates_hydro.m)."""
    a = np.deg2rad(angle_deg)
    vx = (np.cos(a) - 1) * rng.random() + 1
    s = np.sqrt(max(1 - vx * vx, 0))
    vz = dims[1] * (2 * s * rng.random() - s)
    vy = np.sign(2 * rng.random() - 1) * np.sqrt(max(1 - vx * vx - vz * vz, 0))
    th, al = -np.arctan2(vy, vx), np.arcsin(np.clip(vz, -1, 1))
    x = pos[:, 0] * np.cos(th) - pos[:, 1] * np.sin(th)
    pos[:, 1] = pos[:, 0] * np.sin(th) + pos[:, 1] * np.cos(th)
    pos[:, 0] = x * np.cos(al) + pos[:, 2] * np.sin(al)
    pos[:, 2] = -x * np.sin(al) + pos[:, 2] * np.cos(al)
    return pos


def encounter_rates_hydro(P=100_000, nsteps=1000, v=0.0, Rb=1.0, Rp=0.1, LX=10.0, dt=0.01,
                          D0=2.54, nu=1.0, tmin=0.001, omega=1.0, dims=(1, 1),
                          tumble_angle_deg=None, avg_window=500, exclusion_radius=None, swirl=1.0,
                          seed=None):
    """Run the hydro encounter simulation. Returns kernel, running average and final positions.

    In the MATLAB the swim speed came from ``v0 = linspace(0, 30, trials)`` with
    trials = 1, i.e. v = 0. Pass ``v`` to set it.

    ``exclusion_radius`` (default Rb + Rp, as in the MATLAB) is the radius inside
    which CTRW jumps are redrawn. Because it equals Rcrit, diffusion alone can
    never produce a contact and the kernel comes only from advection, so it is
    close to zero. Set ``exclusion_radius=Rb`` to let diffusion reach the capture
    shell. ``swirl`` scales the rotational flow (see ``hydro_advect``). Tumbling was disabled
    (``0*Time >= run_time``); pass ``tumble_angle_deg`` to enable a reorientation
    of up to that angle at exponential times with rate ``omega``.
    """
    rng = np.random.default_rng(seed)
    half = LX / 2
    V = LX ** 3
    Rcrit = Rb + Rp
    mask = np.array([1.0, dims[0], dims[1]])
    U = np.array([-v, 0.0, 0.0])

    lo, hi = -half + Rp, half - Rp
    pos = np.empty((P, 3)); k = 0
    while k < P:
        c = rng.uniform(lo, hi, (max(2 * (P - k), 1000), 3)) * mask
        c = c[np.linalg.norm(c, axis=1) >= Rcrit][: P - k]
        pos[k:k + len(c)] = c; k += len(c)
    initial = pos.copy()
    prev = pos.copy()
    t_wait = np.zeros(P)
    active = np.ones(P, bool)
    res_time = np.zeros(P)
    kernel = np.zeros(nsteps - 1)
    kernel_avg = np.full(nsteps - 1, np.nan)
    dep_vels, res_times = [], []
    deposited = np.zeros(nsteps - 1)
    run_time = 0.0
    time = 0.0
    for n in range(nsteps - 1):
        prev[:] = pos
        pos = hydro_advect(pos, U, Rb, dt, mask, swirl)
        pos, t_wait = hydro_noise(pos, t_wait, time, D0, dt, tmin, nu,
                                  Rcrit if exclusion_radius is None else exclusion_radius, mask, rng)
        # N_dep of encounter_fun_hydro.m: share of phage within Rcrit + Rp (the MATLAB wrote
        # 100*n/10000, i.e. a percentage for P = 10000)
        deposited[n] = np.mean(np.linalg.norm(pos, axis=1) <= Rcrit + Rp)
        Nenc, Nact, dv, rt = hydro_encounter(pos, prev, active, res_time, Rcrit, Rp, dt)
        dep_vels.append(dv); res_times.append(rt)
        kernel[n] = Nenc / Nact * V / dt if Nact else np.nan
        if n >= avg_window:
            kernel_avg[n] = np.nanmean(kernel[n - avg_window:n + 1])
        if tumble_angle_deg is not None and time >= run_time:
            pos = _tumble_rotation(pos, rng, tumble_angle_deg, dims)
            run_time = time - np.log(rng.random()) / omega
        hydro_boundaries(pos, prev, active, half, Rp, mask, rng)
        time += dt
    return dict(kernel=kernel, kernel_avg=kernel_avg, pos=pos, initial=initial, dt=dt, Rcrit=Rcrit,
                v=v, D0=D0, nu=nu, omega=omega, deposited_fraction=deposited,
                deposition_velocities=np.concatenate(dep_vels) if dep_vels else np.empty(0),
                residence_times=np.concatenate(res_times) if res_times else np.empty(0))


def plot_concentration(res, rmax=2.5, bins=40, ax=None, **kw):
    """Radial phage concentration around the bacterium (sim_plotter.m)."""
    import matplotlib.pyplot as plt
    from ..model.plotting import radial_concentration

    ax = ax or plt.gca()
    r0, c0 = radial_concentration(res["initial"], rmax, bins)
    r1, c1 = radial_concentration(res["pos"], rmax, bins)
    ax.plot(r0, c0, "k--", lw=2, label="initial")
    ax.plot(r1, c1, "r", lw=2, label="final")
    ax.set(xlabel="radius", ylabel="phage concentration")
    ax.legend()
    return ax


def hydro_overlap(pos, Rb, warn=True):
    """overlap.m (hydro): indices of phage whose centre lies inside the bacterium (r < Rb)."""
    bad = np.flatnonzero(np.linalg.norm(pos, axis=1) < Rb)
    if warn and bad.size:
        import warnings
        warnings.warn(f"{bad.size} phage inside the bacterium", RuntimeWarning, stacklevel=2)
    return bad


def plot_overlap(pos, prev, Rb, ax=None, lim=2.0):
    """The debug figure of overlap.m: current (red) and previous (green) x-z positions of overlapping phage."""
    import matplotlib.pyplot as plt
    ax = ax or plt.gca()
    bad = hydro_overlap(pos, Rb, warn=False)
    ax.plot(pos[bad, 0], pos[bad, 2], "ro", ms=12, mfc="none", label="current")
    ax.plot(prev[bad, 0], prev[bad, 2], "go", ms=12, mfc="none", label="previous")
    ax.set(xlim=(-lim, lim), ylim=(-lim, lim), xlabel="x", ylabel="z", aspect="equal")
    return ax


def rotate_to_xy_plane(pos):
    """sim_plotter.m: rotate each phage about the x axis onto the x-y plane.

    The flow is axisymmetric about the swimming (x) axis, so this folds the 3D cloud
    onto (x, +/- rho) with rho = sqrt(y^2 + z^2). As in the MATLAB (theta = -atan(z/y)),
    the sign of y is kept.
    """
    pos = np.asarray(pos, float)
    rho = np.hypot(pos[:, 1], pos[:, 2])
    sgn = np.where(pos[:, 1] < 0, -1.0, 1.0)
    return np.column_stack([pos[:, 0], sgn * rho, np.zeros(len(pos))])


def plot_cloud(res, Rb=1.0, LX=None, ax=None):
    """sim_plotter.m figure 2: the phage cloud folded onto the x-y plane, with the cell outline."""
    import matplotlib.pyplot as plt
    ax = ax or plt.gca()
    xy = rotate_to_xy_plane(res["pos"])
    ax.plot(xy[:, 0], xy[:, 1], "k.", ms=2)
    t = np.linspace(0, 2 * np.pi, 400)
    ax.plot(Rb * np.cos(t), Rb * np.sin(t), "m", lw=2)
    if LX:
        ax.set(xlim=(-LX / 2, LX / 2), ylim=(-LX / 2, LX / 2))
    ax.set(xlabel="x", ylabel="y (folded)", aspect="equal")
    return ax
