"""TIGER: kinematic model of a swimming E. coli (rotating body + helical flagellum).

Port of Phage_Bacteria_1/bacteria_model/TIGER.m (same file as CODE/MATLAB/TIGER.m).

Body: 7 meridians of 10 traction points on a prolate half-ellipsoid of
semi-axis 3 Rb (x from -6 Rb to 0) and radius Rb, spaced by 60 degrees. Each
frame the body is rotated by -pi/24 about the x-axis (counter-rotation).
Flagellum: helix along +x with growing amplitude
  y = h E(x) sin(k x - w t),  z = h E(x) cos(k x - w t),  E(x) = 1 - exp(-kE^2 x^2).
"""
from __future__ import annotations

import numpy as np


def body_points(Rb=1.0, n_x=10, n_meridians=7):
    """Traction points on the body (x, y, z): 3 x (n_x * n_meridians)."""
    x0 = np.linspace(-2 * (3 * Rb), 0, n_x)
    hz = Rb * np.sqrt(np.clip(1 - ((x0 + 3 * Rb) / (3 * Rb)) ** 2, 0, None))
    yz0 = np.vstack([0 * x0, hz])
    th = np.pi / 3
    R = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
    xs, yzs = [x0], [yz0]
    for _ in range(n_meridians - 1):
        yz0 = R @ yz0
        xs.append(x0); yzs.append(yz0)
    yz = np.hstack(yzs)
    return np.vstack([np.concatenate(xs), yz])


def flagellum(t, h=0.77, kE=1.3, k=1.3, w=-100.0, length=7.0, n=25):
    """Helical flagellum points (3 x n) at time t."""
    x = np.linspace(0, length, n)
    E = 1 - np.exp(-(kE ** 2) * x ** 2)
    return np.vstack([x, h * E * np.sin(k * x - w * t), h * E * np.cos(k * x - w * t)])


def tiger_frames(n_frames=100, dt=0.05, Rb=1.0, body_rot=-np.pi / 24, **flag_kw):
    """Yield (t, body_xyz, flagellum_xyz) for each frame, as in TIGER.m's loop."""
    body = body_points(Rb)
    R = np.array([[np.cos(body_rot), -np.sin(body_rot)], [np.sin(body_rot), np.cos(body_rot)]])
    t = 0.0
    for _ in range(n_frames):
        body = body.copy()
        body[1:] = R @ body[1:]
        f = flagellum(t, **flag_kw)
        t += dt
        yield t, body, f


def animate_tiger(path="TIGER.mp4", n_frames=100, dt=0.05, Rb=1.0, fps=20, **flag_kw):
    """Render the TIGER animation (replaces movie2avi). Use .gif if ffmpeg is missing."""
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation

    frames = list(tiger_frames(n_frames, dt, Rb, **flag_kw))
    fig = plt.figure(figsize=(7, 6))
    ax = fig.add_subplot(projection="3d")

    def draw(i):
        ax.cla()
        t, b, f = frames[i]
        ax.plot(*b, "ko--", ms=4, lw=1)
        ax.plot(*f, "ko--", ms=4, lw=1)
        ax.set(xlim=(-6 * Rb, 7), ylim=(-2, 2), zlim=(-2, 2), xlabel="L_x", ylabel="L_y", zlabel="L_z",
               title=f"TIGER time: {t:.2f} s")
        ax.view_init(-22, 10)

    anim = FuncAnimation(fig, draw, frames=len(frames))
    anim.save(path, fps=fps)
    plt.close(fig)
    return path
