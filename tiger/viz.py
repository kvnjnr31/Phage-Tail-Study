"""TIGER x receptor-model movie renderer.

Folded in from 04_Phage_Adsorption_Tail_Length/tiger_encounter_coupled.py
(July 2026). That script carried its own copy of the receptor physics; here
the physics comes from :func:`tiger.model.run_receptor_simulation` (run with
``record_every``) and only the rendering is kept:

* the capsule body, its receptors and attached/bound phage rotate rigidly about
  the y-axis (the capsule's symmetry axis, so rotation is an exact symmetry of
  the receptor geometry and does not change the binding statistics);
* free phage are drawn in the lab frame;
* a TIGER.m-style travelling-wave flagellum is attached at the tail pole and
  counter-rotates.

Rotation is a rendering device only: it is not coupled to the phage motion. Use
``tiger.flow`` to couple swimming hydrodynamics to the physics.
"""
from __future__ import annotations

import shutil

import numpy as np


def rot_y(points: np.ndarray, theta: float) -> np.ndarray:
    c, s = np.cos(theta), np.sin(theta)
    out = np.array(points, float, copy=True)
    out[:, 0] = points[:, 0] * c - points[:, 2] * s
    out[:, 2] = points[:, 0] * s + points[:, 2] * c
    return out


def capsule_mesh(Rb: float, L: float, n_theta=24, n_y=5, n_cap=6):
    th = np.linspace(0, 2 * np.pi, n_theta)
    Th, Y = np.meshgrid(th, np.linspace(-L, L, n_y))
    cyl = (Rb * np.cos(Th), Y, Rb * np.sin(Th))

    def hemi(sign):
        Th2, Ph = np.meshgrid(th, np.linspace(0, np.pi / 2, n_cap))
        return (Rb * np.cos(Th2) * np.sin(Ph), sign * (L + Rb * np.cos(Ph)), Rb * np.sin(Th2) * np.sin(Ph))

    return cyl, hemi(+1), hemi(-1)


def flagellum_points(t, Rb, L, length=3.5, h=0.35, kE=1.3, k=1.3, w=-18.0, n=60):
    """TIGER.m travelling-wave flagellum, attached at the tail pole and extended along -y."""
    s = np.linspace(0, length, n)
    E = 1 - np.exp(-(kE ** 2) * s ** 2)
    return np.stack([h * E * np.sin(k * s - w * t), -(L + Rb) - s, h * E * np.cos(k * s - w * t)], axis=1)


def render_tiger_movie(res, out_path="tiger_encounter.mp4", body_revolutions=1.5, flagella_revolutions=-9.0,
                       fps=30, dpi=130, lim=(3.0, 4.5, 3.0), flag_kw=None):
    """Render the frames recorded in ``res`` (run_receptor_simulation(..., record_every=k)).

    Writes MP4 if ffmpeg is available, otherwise a GIF next to ``out_path``.
    Returns the path written.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import animation

    frames = res.frames
    if not frames:
        raise ValueError("no frames: run run_receptor_simulation(..., record_every=k)")
    Rb, L, P = res.params["Rb"], res.params["bac_length"], res.params["P"]
    tail = res.params["Ltail"]
    meshes = capsule_mesh(Rb, L)
    flag_kw = flag_kw or {}
    total_t = frames[-1]["t"]
    w_body = 2 * np.pi * body_revolutions / max(total_t, 1e-9)
    w_flag = 2 * np.pi * flagella_revolutions / max(total_t, 1e-9)

    fig = plt.figure(figsize=(9, 8))
    ax = fig.add_subplot(111, projection="3d")

    def draw(i):
        ax.clear()
        f = frames[i]
        t = f["t"]
        th = w_body * t
        for X, Y, Z in meshes:
            pts = rot_y(np.stack([X.ravel(), Y.ravel(), Z.ravel()], 1), th)
            ax.plot_surface(*(pts[:, j].reshape(X.shape) for j in range(3)), color="steelblue",
                            alpha=0.25, linewidth=0, shade=False)
        rr = rot_y(res.receptors, th)
        ax.scatter(*rr.T, c="k", s=3, alpha=0.5)
        fl = rot_y(flagellum_points(t, Rb, L, **flag_kw), th + w_flag * t)
        ax.plot(*fl.T, c="dimgray", lw=2.5)
        pos, att, bnd = f["pos"], f["attached"], f["bound"]
        for mask, rotate, kw in ((~att & ~bnd, False, dict(c="orange", s=18, label="free")),
                                 (~att & bnd, False, dict(c="gold", s=22, marker="^", label="bound (bulk)")),
                                 (att & ~bnd, True, dict(c="crimson", s=22, label="attached (searching)")),
                                 (att & bnd, True, dict(c="lime", s=40, edgecolor="k", label="bound"))):
            if mask.any():
                pp = rot_y(pos[mask], th) if rotate else pos[mask]
                ax.scatter(*pp.T, **kw)
        ax.set(xlim=(-lim[0], lim[0]), ylim=(-lim[1], lim[1]), zlim=(-lim[2], lim[2]),
               xlabel="x", ylabel="y (body axis)", zlabel="z")
        ax.set_box_aspect((1, lim[1] / lim[0], 1))
        ax.view_init(elev=18, azim=-60 + 25 * np.sin(0.3 * t))
        ax.set_title(f"TIGER receptor model | tail = {tail * 1000:.0f} nm | t = {t:5.2f} s | "
                     f"bound = {int(bnd.sum())}/{P}")
        if i == 0:
            ax.legend(loc="upper left", fontsize=8)
        return []

    anim = animation.FuncAnimation(fig, draw, frames=len(frames), interval=1000 / fps, blit=False)
    if shutil.which("ffmpeg") and str(out_path).endswith(".mp4"):
        anim.save(out_path, writer=animation.FFMpegWriter(fps=fps, bitrate=3000), dpi=dpi)
        written = str(out_path)
    else:
        written = str(out_path).rsplit(".", 1)[0] + ".gif"
        anim.save(written, writer=animation.PillowWriter(fps=min(fps, 20)), dpi=min(dpi, 90))
    plt.close(fig)
    return written
