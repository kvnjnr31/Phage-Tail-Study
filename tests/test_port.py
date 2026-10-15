"""Sanity tests for the MATLAB -> Python port. Run:  python -m pytest tests -q"""
import numpy as np
import pytest

from tiger.model import (init_receptors, distance_to_axis, level_set, outward_normal, rattle_step,
                         run_receptor_simulation, ReceptorSimParams, project_to_surface)
from tiger.encounter import (encounter_rates_ballistic, ballistic_kernel, diffusion_times,
                             run_tumble_tracks, ensemble_msd, run_tumble_D, ctrw_tracks, fit_power_law,
                             stokes_disturbance)
from tiger.analysis import first_passage, ols


def test_receptors_on_surface():
    B = init_receptors(300, 0.05, 0.05 * np.pi, 0.8, 1.9, 1e-4, np.random.default_rng(0))
    assert np.allclose(distance_to_axis(B, 1.9), 0.8, atol=1e-9)
    d = np.linalg.norm(B[:, None] - B[None], axis=2) + np.eye(len(B))
    assert d.min() >= 2e-4


def test_projection():
    p = np.random.default_rng(1).normal(size=(100, 3)) * 3
    q = project_to_surface(p, 0.8, 1.9)
    assert np.allclose(distance_to_axis(q, 1.9), 0.8)


def test_rattle_keeps_constraint():
    B = init_receptors(200, 0.05, 0.05 * np.pi, 0.8, 1.9, 1e-4, np.random.default_rng(0))
    x, p = np.array([0.8, 0.3, 0.0]), np.array([0.0, 1e-16, 2e-16])
    for _ in range(200):
        x, p = rattle_step(x, p, B, 0.8, 1.9, 1e-4)
    assert abs(level_set(x, 0.8, 1.9)[0]) < 1e-10
    assert abs(p @ outward_normal(x, 1.9)) < 1e-25


def test_receptor_sim_runs_and_attached_stay_on_shell():
    r = run_receptor_simulation(ReceptorSimParams(nsteps=100, P=50, B=300, seed=2))
    ph = r.phage
    a = ph.attached
    R = 0.8 + (ph.tail_length + ph.head_radius) / 2
    assert np.allclose(distance_to_axis(ph.pos[a], 1.9), R[a])
    assert (distance_to_axis(ph.pos, 1.9) >= 0.8).all()


def test_ballistic_kernel_matches_theory():
    r = encounter_rates_ballistic(P=50_000, nsteps=60, seed=0)
    assert abs(np.mean(r["kernel"][5:]) / r["theory"] - 1) < 0.15


def test_capture_time_theory():
    r = diffusion_times(P=20_000, nsteps=500, seed=0)
    assert abs(np.nanmean(r["capture_times"]) / r["theory"] - 1) < 0.1


def test_run_tumble_diffusivity():
    tr = run_tumble_tracks(rate=10, v=1, dt=0.01, steps=3000, N=2000, seed=0)
    tau, m = ensemble_msd(tr, 0.01)
    slope = np.polyfit(tau[1000:], m[1000:], 1)[0]
    assert abs(slope / (4 * run_tumble_D(1, 10)) - 1) < 0.1


def test_ctrw_exponent():
    tau, m = ensemble_msd(ctrw_tracks(0.6, 1.0, 0.01, 2000, 500, seed=0), 0.01, origin=np.zeros(2))
    a, _ = fit_power_law(tau, m, slice(4, None, 10))
    assert abs(a - 0.6) < 0.08


def test_stokes_no_slip():
    U = np.array([-30.0, 0, 0])
    pts = np.array([[1.0, 0, 0], [0, 1.0, 0], [0, 0, -1.0], [0.6, 0.8, 0]])
    assert np.allclose(stokes_disturbance(pts, U, 1.0), U)


def test_first_passage_positive_and_monotone():
    t1, t2 = first_passage(f_rup_pN=-60), first_passage(f_rup_pN=-100)
    assert 0 < t1 < t2


def test_ols():
    x = np.arange(10.0)
    f = ols(2 + 3 * x, x)
    assert np.allclose(f.beta, [2, 3]) and f.r2 == pytest.approx(1.0)
