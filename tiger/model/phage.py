"""Phage population state.

The MATLAB code stores phage in an N x 11 matrix ``MatP``. Here the columns
are named fields of :class:`PhageState`:

=======  =============================  ==============================
MatP col PhageState field               meaning
=======  =============================  ==============================
1-3      pos[:, 0:3]                    current position (x, y, z)
4-6      prev[:, 0:3]                   previous / initial position
                                        (MATLAB also called these the
                                        "drift velocity" columns)
7        t_wait                         CTRW time of next allowed jump
8        receptor                       index of bound receptor, -1 if
                                        unbound (MATLAB: 0 = unbound,
                                        k = 1-based receptor index)
9        attached                       True once on the bacterium surface
10       tail_length                    tail length (um)
11       head_radius                    capsid radius Rp (um)
=======  =============================  ==============================
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class PhageState:
    pos: np.ndarray
    prev: np.ndarray
    t_wait: np.ndarray
    receptor: np.ndarray
    attached: np.ndarray
    tail_length: np.ndarray
    head_radius: np.ndarray

    @property
    def n(self) -> int:
        return self.pos.shape[0]

    @property
    def bound(self) -> np.ndarray:
        return self.receptor >= 0

    def copy(self) -> "PhageState":
        return PhageState(*(getattr(self, f).copy() for f in
                            ("pos", "prev", "t_wait", "receptor", "attached",
                             "tail_length", "head_radius")))

    def to_matp(self) -> np.ndarray:
        """Return the MATLAB-layout N x 11 matrix (1-based receptor index, 0 = unbound)."""
        M = np.zeros((self.n, 11))
        M[:, 0:3] = self.pos
        M[:, 3:6] = self.prev
        M[:, 6] = self.t_wait
        M[:, 7] = self.receptor + 1
        M[:, 8] = self.attached
        M[:, 9] = self.tail_length
        M[:, 10] = self.head_radius
        return M


def init_phage(P: int, Rb: float, Rp: float, box: float, Ltail: float,
               shell: float = 0.05, rng: np.random.Generator | None = None) -> PhageState:
    """Port of init_phage.m.

    Phage are drawn uniformly in the cube [-box/2, box/2]^3 (shrunk by Rp) and
    kept only if their distance from the bacterium's long (y) axis lies in
    [Rb + Rp, Rb + Rp + shell]; i.e. they start in a thin shell around the
    (infinite) cylinder, as in the MATLAB code. All start detached and unbound.
    """
    rng = np.random.default_rng() if rng is None else rng
    lo, hi = -box / 2 + Rp, box / 2 - Rp
    pos = np.empty((P, 3))
    filled = 0
    while filled < P:
        cand = rng.uniform(lo, hi, size=(max(4 * (P - filled), 1000), 3))
        r = np.hypot(cand[:, 0], cand[:, 2])
        ok = cand[(r >= Rb + Rp) & (r <= Rb + Rp + shell)]
        take = ok[: P - filled]
        pos[filled: filled + len(take)] = take
        filled += len(take)
    return PhageState(
        pos=pos,
        prev=pos.copy(),
        t_wait=np.zeros(P),
        receptor=-np.ones(P, dtype=int),
        attached=np.zeros(P, dtype=bool),
        tail_length=np.full(P, float(Ltail)),
        head_radius=np.full(P, float(Rp)),
    )
