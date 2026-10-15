import numpy as np

# === Minimal LBM Scaffold for TIGER Simulation ===

# D3Q19 lattice
c = np.array([[ 0, 0, 0], [1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0],
              [0, 0, 1], [0, 0, -1], [1, 1, 0], [-1, -1, 0], [1, -1, 0],
              [-1, 1, 0], [1, 0, 1], [-1, 0, -1], [1, 0, -1], [-1, 0, 1],
              [0, 1, 1], [0, -1, -1], [0, 1, -1], [0, -1, 1]])

w = np.array([1/3] + [1/18]*6 + [1/36]*12)


def initialize_lbm(nx, ny, nz, tau):
    f = np.ones((nx, ny, nz, 19)) * (1 / 19)
    rho = np.ones((nx, ny, nz))
    u = np.zeros((nx, ny, nz, 3))
    return f, rho, u, c, w


def collision(f, rho, u, c, w, tau):
    nx, ny, nz = rho.shape
    feq = np.zeros_like(f)
    for i in range(19):
        cu = np.einsum('ijk,->ijk', u, c[i])
        feq[..., i] = w[i] * rho * (1 + 3*cu + 9/2*cu**2 - 3/2*np.sum(u**2, axis=-1))
    f += -(1 / tau) * (f - feq)
    return f


def streaming(f, c):
    f_new = np.zeros_like(f)
    for i in range(19):
        f_new[..., i] = np.roll(f[..., i], shift=c[i], axis=(0, 1, 2))
    return f_new


def macroscopic(f, c):
    rho = np.sum(f, axis=-1)
    u = np.einsum('...i,ij->...j', f, c)
    u /= rho[..., None] + 1e-12  # avoid divide by zero
    return rho, u
