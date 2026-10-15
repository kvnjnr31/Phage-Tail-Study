# hybrid_simulation.py
# TIGER–LBM with Guo forcing + smoothed deposition
# Physically calibrated to an E. coli–like swimmer (Watari & Larson / Berg-style).
#
# Run:
#   python simulation/hybrid_simulation.py
#
# Outputs:
#   data/numerical/velocity_field.npy
#   data/numerical/velocity_field_small.npy
#   data/numerical/bacterium_path.npy
#   data/numerical/tail_coords.npy
#
# All key physical + forcing knobs are now in tiger_simulation_config.yaml
# under `physical`, `bacterium`, and `forcing`.

import os
import numpy as np
import yaml

# ---------------------------------------------------------------------
# 1. LATTICE + PHYSICAL CONFIG (from YAML)
# ---------------------------------------------------------------------

CONFIG_PATH = os.path.join(os.path.dirname(__file__), "tiger_simulation_config.yaml")
with open(CONFIG_PATH, "r") as f:
    config = yaml.safe_load(f)

# Grid / fluid / time-stepping
nx    = int(config["grid"]["nx"])
ny    = int(config["grid"]["ny"])
nz    = int(config["grid"]["nz"])

tau   = float(config["fluid"]["tau"])              # LBM relaxation time

steps = int(config["simulation"]["steps"])         # number of time steps
dt    = float(config["simulation"]["dt"])          # simulation time step (lattice units)

# Bacterium geometry in lattice units
RB_LB = float(config["bacterium"]["Rb"])           # body radius in lattice units

# Physical parameters (Watari & Larson / Berg-style)
phys_cfg = config["physical"]

RB_PHYS_UM          = float(phys_cfg.get("Rb_um", 1.0))           # body radius [µm] (~1 µm)
HELIX_RADIUS_FACTOR = float(phys_cfg.get("helix_radius_factor", 0.25))
HELIX_PITCH_FACTOR  = float(phys_cfg.get("helix_pitch_factor", 0.60))

NU_PHYS_TARGET      = float(phys_cfg.get("nu_phys_target", 1.0e-6))   # water ν [m^2/s]
RHO_PHYS            = float(phys_cfg.get("rho_phys", 1000.0))         # density [kg/m^3]

# Target physical bundle rotation (Hz) and measured lattice frequency
F_BUNDLE_TARGET_HZ  = float(phys_cfg.get("bundle_target_Hz", 100.0))  # e.g. 100–130 Hz (Berg)
F_BUNDLE_LB_MEAS    = float(phys_cfg.get("f_bundle_lb_measured", 1.6978e-2))  # [1 / tu]

# Geometric scale: Δx from matching body radius
DX_PHYS_M = (RB_PHYS_UM * 1e-6) / RB_LB    # meters per lattice cell

# ---------------------------------------------------------------------
# 2. FORCING / GEOMETRY KNOBS (dimensionless)
# ---------------------------------------------------------------------

forcing_cfg = config["forcing"]

# Source models
USE_STRESSLET            = bool(forcing_cfg.get("use_stresslet", True))
USE_ROTLET_RING_BODY     = bool(forcing_cfg.get("use_rotlet_ring_body", True))
USE_HELICAL_BUNDLE       = bool(forcing_cfg.get("use_helical_bundle", True))

USE_TAIL                 = bool(forcing_cfg.get("use_tail", True))
FORCE_FREE_TAIL          = bool(forcing_cfg.get("force_free_tail", True))
HEAD_DRAG_OFFSET         = float(forcing_cfg.get("head_drag_offset", -0.25))  # in units of Rb
AUTO_DISABLE_STRESSLET_WHEN_TAIL = bool(
    forcing_cfg.get("auto_disable_stresslet_when_tail", True)
)

# Tail forcing
TAIL_SAMPLES   = int(forcing_cfg.get("tail_samples", 200))
FT_TAIL        = float(forcing_cfg.get("ft_tail", 1.0e-1))  # per-sample thrust along +x

# Stresslet parameters
F_STRESS       = float(forcing_cfg.get("f_stress", 0.05))   # stresslet force magnitude
D_STRESS       = forcing_cfg.get("d_stress", None)          # separation; if None → 0.8*Rb

# Torque parameters (shared by bundle and body ring)
TORQUE_GAMMA   = float(forcing_cfg.get("torque_gamma", 1.0))  # total |Γ| for the bundle

SIGMA          = float(forcing_cfg.get("sigma", 2.2))        # Gaussian deposition width [lu]
STENCIL_R      = int(forcing_cfg.get("stencil_r", 3))        # deposition half-width

# Helical bundle geometry (about x̂), from Watari & Larson ratios
N_HELICES      = int(forcing_cfg.get("n_helices", 1))        # 1 or 2 helices (phase-shifted)
N_TURN_INTEGER = int(forcing_cfg.get("n_turn_integer", 3))   # integer number of turns
PTS_PER_TURN   = int(forcing_cfg.get("pts_per_turn", 120))   # resolution along each helix

# Ratios relative to Rb; converted to lattice units below
XOFF_HELIX_FACTOR = float(forcing_cfg.get("xoff_helix_factor", 0.6))  # axial offset in Rb units

# Bundle / body angular frequencies in lattice units (rad / tu)
OMEGA_BUNDLE = float(forcing_cfg.get("omega_bundle", 6.0))
OMEGA_BODY   = float(forcing_cfg.get("omega_body", -6.0))

# Body counter-ring geometry
RING_R_BODY_FACTOR = float(forcing_cfg.get("ring_r_body_factor", 0.55))
XOFF_BODY          = float(forcing_cfg.get("xoff_body", 0.0))
N_TORQUE_PTS       = int(forcing_cfg.get("n_torque_pts", 48))

# ---------------------------------------------------------------------
# 3. DERIVED LATTICE SCALES & PHYSICAL MAPPING
# ---------------------------------------------------------------------

CS2   = 1.0 / 3.0
NU_LB = CS2 * (tau - 0.5)                     # lattice viscosity [lu^2 / tu]

# seconds per lattice time unit:
#   f_phys = f_lb / Δt  →  Δt = f_lb / f_phys
DT_UNIT_S = abs(F_BUNDLE_LB_MEAS) / F_BUNDLE_TARGET_HZ   # [s / lattice time unit]

# physical seconds per simulation *step* (dt is in lattice time units)
DT_PHYS_S = dt * DT_UNIT_S                                   # [s / simulation step]

# implied physical viscosity from LBM
NU_PHYS_FROM_LBM = NU_LB * DX_PHYS_M**2 / DT_UNIT_S          # [m^2/s]

# Helical geometry in lattice units, using Watari & Larson ratios
# Watari & Larson: helix radius ≈ 0.25 µm, pitch ≈ 2.5 µm; body radius ≈ 1 µm.
R_HELIX_LB        = HELIX_RADIUS_FACTOR * RB_LB
PITCH_PER_TURN_LB = HELIX_PITCH_FACTOR * RB_LB
XOFF_HELIX_LB     = XOFF_HELIX_FACTOR * RB_LB
RING_R_BODY_LB    = RING_R_BODY_FACTOR * RB_LB

def print_calibration_summary():
    """
    Print how the current lattice parameters map to SI units (Watari & Larson / Berg-style).
    """
    print("\n=== LBM → SI calibration (E. coli run-and-tumble) ===")
    print(f"Grid size              : nx={nx}, ny={ny}, nz={nz}")
    print(f"Rb (lattice units)     : {RB_LB:.3f}")
    print(f"Rb (physical)          : {RB_PHYS_UM:.3f} µm")
    print(f"Δx                     : {DX_PHYS_M*1e6:.3f} µm / cell")
    print(f"tau (LBM)              : {tau:.4f}")
    print(f"ν_lattice              : {NU_LB:.4e} [lu^2/tu]")
    print(f"Measured f_lb          : {F_BUNDLE_LB_MEAS:.4e} [1/tu]")
    print(f"Target f_phys          : {F_BUNDLE_TARGET_HZ:.2f} [Hz]")
    print(f"Δt (phys per tu)       : {DT_UNIT_S:.4e} s/tu")
    print(f"dt (sim step)          : {dt:.3f} tu → {DT_PHYS_S:.4e} s per step")
    print(f"ν_phys (from LBM)      : {NU_PHYS_FROM_LBM:.4e} m^2/s")
    print(f"ν_phys (target)        : {NU_PHYS_TARGET:.4e} m^2/s")
    print(f"ρ_phys (assumed)       : {RHO_PHYS:.1f} kg/m^3")
    print(f"Helix radius (lu)      : {R_HELIX_LB:.3f} ≈ {HELIX_RADIUS_FACTOR:.2f} Rb")
    print(f"Helix pitch/turn (lu)  : {PITCH_PER_TURN_LB:.3f} ≈ {HELIX_PITCH_FACTOR:.2f} Rb")
    print("======================================================\n")

# ---------------------------------------------------------------------
# 4. TIGER BODY (HEAD + FLEXIBLE TAIL)
# ---------------------------------------------------------------------

class TIGERBacterium:
    """
    Simple TIGER body with a rotating helical tail.
    Geometry is defined in lattice units; physical size is set
    by RB_LB and RB_PHYS_UM at the top of the file.
    """
    def __init__(self, center, Rb, kE, k, w, h):
        self.Rb   = Rb
        self.kE   = kE
        self.k    = k
        self.w    = w
        self.h    = h

        # Tail parametric x coordinate and exponential envelope
        self.x    = np.linspace(0, 7, 25)
        self.Ex   = 1.0 - np.exp(-(kE**2) * self.x**2)

        self.center   = np.array(center, dtype=float)
        self.theta    = np.pi/3
        self.theta_step = -np.pi/24

        self.velocity = np.zeros(3)
        self._build_head_arc()

    def _build_head_arc(self):
        # Build a spheroidal head arc around the bacterium center in yz-plane.
        head_coord_x0 = np.linspace(-2*(3*self.Rb), 0, 10)
        hy = np.zeros_like(head_coord_x0)
        hz = self.Rb * np.sqrt(1.0 - ((head_coord_x0 + (3*self.Rb))/(3*self.Rb))**2)

        head_theta0 = np.vstack((hy, hz))
        head_coord_yz = head_theta0.copy()

        R_theta = np.array([
            [np.cos(self.theta), -np.sin(self.theta)],
            [np.sin(self.theta),  np.cos(self.theta)],
        ])

        for _ in range(6):
            head_theta0   = R_theta @ head_theta0
            head_coord_yz = np.hstack((head_coord_yz, head_theta0))

        self.head_coord_x  = np.tile(head_coord_x0, 7)
        self.head_coord_yz = head_coord_yz

        self.R_theta_time = np.array([
            [np.cos(self.theta_step), -np.sin(self.theta_step)],
            [np.sin(self.theta_step),  np.cos(self.theta_step)],
        ])

    def update_geometry(self, t):
        # Rotate head cross-section in time
        self.head_coord_yz = self.R_theta_time @ self.head_coord_yz

        # Tail: traveling helix with exponential envelope
        self.tail_y = self.h * self.Ex * np.sin(self.k*self.x - self.w*t)
        self.tail_z = self.h * self.Ex * np.cos(self.k*self.x - self.w*t)
        self.tail_x = self.x

    def advect(self, u_field, dt_local):
        idx = tuple(np.clip(self.center.astype(int), [0, 0, 0],
                            np.array(u_field.shape[:3]) - 1))
        fluid_u = u_field[idx]
        if not np.all(np.isfinite(fluid_u)):
            fluid_u = np.zeros_like(fluid_u)
        gamma = 0.5
        self.velocity += (-gamma * (self.velocity - fluid_u)) * dt_local
        self.center   += self.velocity * dt_local

# ---------------------------------------------------------------------
# 5. LBM CORE (D3Q19)
# ---------------------------------------------------------------------

c = np.array([
    [ 0, 0, 0],
    [ 1, 0, 0], [-1, 0, 0],
    [ 0, 1, 0], [ 0,-1, 0],
    [ 0, 0, 1], [ 0, 0,-1],
    [ 1, 1, 0], [-1,-1, 0],
    [ 1,-1, 0], [-1, 1, 0],
    [ 1, 0, 1], [-1, 0,-1],
    [ 1, 0,-1], [-1, 0, 1],
    [ 0, 1, 1], [ 0,-1,-1],
    [ 0, 1,-1], [ 0,-1, 1]
], dtype=float)

w = np.array([1/3] + [1/18]*6 + [1/36]*12, dtype=float)

def initialize_lbm(nx_local, ny_local, nz_local, tau_local):
    f = np.zeros((nx_local, ny_local, nz_local, 19), dtype=float)
    rho = np.ones((nx_local, ny_local, nz_local), dtype=float)
    u   = np.zeros((nx_local, ny_local, nz_local, 3), dtype=float)
    for i in range(19):
        f[..., i] = w[i] * rho
    return f, rho, u

def streaming(f):
    f_new = np.zeros_like(f)
    for i in range(19):
        f_new[..., i] = np.roll(f[..., i], shift=c[i].astype(int), axis=(0, 1, 2))
    return f_new

def macroscopic(f):
    rho = np.sum(f, axis=-1)
    mom = np.einsum("...i,ij->...j", f, c)
    u   = mom / (rho[..., None] + 1e-12)
    return rho, u

# ---------------------------------------------------------------------
# 6. FORCING UTILITIES (Guo forcing)
# ---------------------------------------------------------------------

def deposit_force(F, pos, vec, sigma=1.0, stencil_r=1):
    """
    Gaussian-smoothed deposition to a (2r+1)^3 stencil around pos (grid coords).
    """
    nx_local, ny_local, nz_local, _ = F.shape
    x0, y0, z0 = map(float, pos)
    # Guard against NaN / inf positions to avoid integer conversion errors
    if (not np.isfinite(x0)) or (not np.isfinite(y0)) or (not np.isfinite(z0)):
        return
    ix, iy, iz = int(round(x0)), int(round(y0)), int(round(z0))

    for dx in range(-stencil_r, stencil_r+1):
        for dy in range(-stencil_r, stencil_r+1):
            for dz in range(-stencil_r, stencil_r+1):
                x = ix + dx
                y = iy + dy
                z = iz + dz
                if 0 <= x < nx_local and 0 <= y < ny_local and 0 <= z < nz_local:
                    r2 = (x - x0)**2 + (y - y0)**2 + (z - z0)**2
                    wgt = np.exp(-0.5 * r2 / (sigma * sigma))
                    F[x, y, z, :] += wgt * vec

def resample_curve(P, n):
    """
    Resample a 3D polyline P (N×3) to n points (uniform in index-arc).
    """
    N = len(P)
    if N == 0:
        return P
    t = np.linspace(0, N-1, N)
    ti = np.linspace(0, N-1, n)
    Pi = np.column_stack([
        np.interp(ti, t, P[:, 0]),
        np.interp(ti, t, P[:, 1]),
        np.interp(ti, t, P[:, 2]),
    ])
    return Pi

def deposit_torque_ring(F, center, R, xoff, torque_gamma, npts, omega, t, sigma, stencil_r):
    """
    Rotlet about x̂ as N tangential forces on a ring at x = x_c + xoff.
    Each point uses F0 = Γ / (N * R).  (r × f)_x = R * F0  → Σ = Γ.
    """
    if R is None or R <= 0.0 or npts < 2:
        return

    F0 = torque_gamma / (npts * R)
    xc = center[0] + xoff
    yc = center[1]
    zc = center[2]

    for j in range(npts):
        phi = 2.0*np.pi*j/npts + omega*t
        y = yc + R*np.cos(phi)
        z = zc + R*np.sin(phi)
        ty = -np.sin(phi)
        tz =  np.cos(phi)   # t̂_φ in YZ
        vec = F0 * np.array([0.0, ty, tz])
        deposit_force(F, (xc, y, z), vec, sigma=sigma, stencil_r=stencil_r)

def deposit_torque_helix(F, center, R, pitch_per_turn, n_turns, nsamp,
                         xoff, torque_gamma, omega, t,
                         n_helices=1, phase0=0.0, sigma=1.0, stencil_r=2):
    """
    Non-axisymmetric rotlet: forces placed on one or more helices about x̂.

    Geometry:
      x(s) = x_c + xoff + s,   φ(s) = 2π s / pitch + ω t + φ0,
      y    = y_c + R cosφ,     z    = z_c + R sinφ

    Force at each sample is tangential around x̂: (0, -sinφ, cosφ).
    Choose F0 = Γ / (Nsamples_total * R) so Σ (r × f)_x ≈ Γ.
    Use integer n_turns to keep net force ~ 0.
    """
    if R is None or R <= 0.0 or nsamp < 2 or n_turns <= 0:
        return

    L = n_turns * pitch_per_turn
    nsamp_total = nsamp * max(1, n_helices)
    F0 = torque_gamma / (nsamp_total * R)

    xc0, yc, zc = center[0] + xoff, center[1], center[2]

    for h in range(n_helices):
        phi0 = phase0 + h * (2.0*np.pi / max(1, n_helices))
        for k in range(nsamp):
            s = (k + 0.5) * (L / nsamp)  # 0..L
            phi = 2.0*np.pi * (s / pitch_per_turn) + omega*t + phi0
            x = xc0 + s - 0.5*L
            y = yc + R*np.cos(phi)
            z = zc + R*np.sin(phi)
            ty = -np.sin(phi)
            tz =  np.cos(phi)
            vec = F0 * np.array([0.0, ty, tz])  # pure azimuthal around x̂
            deposit_force(F, (x, y, z), vec, sigma=sigma, stencil_r=stencil_r)

# ---------------------------------------------------------------------
# 7. FORCING ASSEMBLY
# ---------------------------------------------------------------------

def build_forces(bact, u_field_shape, t):
    """
    Assemble body-force field F for this step.

    Includes:
      - Helical bundle torque (Watari & Larson-style)
      - Body counter-rotlet ring (to keep total torque ≈ 0)
      - Optional tail thrust (force-free option)
      - Optional stresslet dipole
    """
    F = np.zeros(u_field_shape, dtype=float)
    Rb = bact.Rb

    # Decide if we keep a separate stresslet when a force-free tail is active
    apply_stresslet = USE_STRESSLET
    if USE_TAIL and FORCE_FREE_TAIL and AUTO_DISABLE_STRESSLET_WHEN_TAIL:
        apply_stresslet = False

    # --- Helical bundle torque (non-axisymmetric) ---
    if USE_HELICAL_BUNDLE:
        deposit_torque_helix(
            F,
            bact.center,
            R_HELIX_LB,                 # from Watari & Larson ratios
            PITCH_PER_TURN_LB,
            N_TURN_INTEGER,
            PTS_PER_TURN,
            XOFF_HELIX_LB,
            TORQUE_GAMMA,
            OMEGA_BUNDLE,
            t,
            n_helices=N_HELICES,
            phase0=0.0,
            sigma=SIGMA,
            stencil_r=STENCIL_R,
        )

    # --- Body counter-ring (−Γ) to ensure total torque = 0 ---
    if USE_ROTLET_RING_BODY:
        deposit_torque_ring(
            F,
            bact.center,
            RING_R_BODY_LB,
            XOFF_BODY,
            -TORQUE_GAMMA,
            N_TORQUE_PTS,
            OMEGA_BODY,
            t,
            sigma=SIGMA,
            stencil_r=STENCIL_R,
        )

    # --- Tail thrust along +x (optionally force-free) ---
    if USE_TAIL:
        tail_curve = np.vstack((
            bact.center[0] + bact.tail_x,
            bact.center[1] + bact.tail_y,
            bact.center[2] + bact.tail_z
        )).T
        tail_curve = resample_curve(tail_curve, TAIL_SAMPLES)

        F_tail_total = np.zeros(3)
        for px, py, pz in tail_curve:
            fvec = np.array([FT_TAIL, 0.0, 0.0])   # thrust along +x
            deposit_force(F, (px, py, pz), fvec, sigma=SIGMA, stencil_r=STENCIL_R)
            F_tail_total += fvec

        if FORCE_FREE_TAIL:
            # equal-and-opposite drag on the head region so net force ≈ 0
            drag_pos = tuple(bact.center + np.array([HEAD_DRAG_OFFSET*Rb, 0.0, 0.0]))
            deposit_force(F, drag_pos, -F_tail_total, sigma=SIGMA, stencil_r=STENCIL_R)

    # --- Stresslet (dipole) along x̂ (optional; avoid when tail provides dipole) ---
    if apply_stresslet:
        d = (0.8*Rb) if (D_STRESS is None) else float(D_STRESS)
        deposit_force(
            F,
            tuple(bact.center + np.array([+d, 0.0, 0.0])),
            np.array([+F_STRESS, 0.0, 0.0]),
            sigma=SIGMA,
            stencil_r=STENCIL_R,
        )
        deposit_force(
            F,
            tuple(bact.center + np.array([-d, 0.0, 0.0])),
            np.array([-F_STRESS, 0.0, 0.0]),
            sigma=SIGMA,
            stencil_r=STENCIL_R,
        )

    return F

# ---------------------------------------------------------------------
# 8. COLLISION (BGK + Guo forcing)
# ---------------------------------------------------------------------

def collision_with_guo(f, rho, u, tau_local, F, dt_local):
    """
    BGK collision step with Guo forcing.
    """
    feq = np.zeros_like(f)
    uu = np.sum(u*u, axis=-1)

    for i in range(19):
        cu = u[..., 0]*c[i, 0] + u[..., 1]*c[i, 1] + u[..., 2]*c[i, 2]
        feq[..., i] = w[i] * rho * (
            1.0
            + (cu/CS2)
            + 0.5*(cu*cu)/(CS2*CS2)
            - 0.5*(uu/CS2)
        )

    f += -(1.0/tau_local) * (f - feq)

    a = F / (rho[..., None] + 1e-12)
    pref = (1.0 - 1.0/(2.0*tau_local)) * dt_local

    for i in range(19):
        ci = c[i]
        ci_dot_u = u[..., 0]*ci[0] + u[..., 1]*ci[1] + u[..., 2]*ci[2]
        ci_dot_a = a[..., 0]*ci[0] + a[..., 1]*ci[1] + a[..., 2]*ci[2]

        term1 = (
            a[..., 0]*(ci[0] - u[..., 0])
            + a[..., 1]*(ci[1] - u[..., 1])
            + a[..., 2]*(ci[2] - u[..., 2])
        ) / CS2
        term2 = (ci_dot_u * ci_dot_a) / (CS2*CS2)

        f[..., i] += pref * w[i] * (term1 + term2)

    return f

# ---------------------------------------------------------------------
# 9. MAIN SIMULATION DRIVER
# ---------------------------------------------------------------------

def simulate_hybrid_lbm_tiger():
    print_calibration_summary()

    # Initialize LBM fields
    f, rho, u = initialize_lbm(nx, ny, nz, tau)

    # Instantiate TIGER bacterium from config
    bact = TIGERBacterium(
        center=np.array(config["bacterium"]["center"], dtype=float),
        Rb=float(config["bacterium"]["Rb"]),
        kE=float(config["bacterium"]["kE"]),
        k=float(config["bacterium"]["k"]),
        w=float(config["bacterium"]["w"]),
        h=float(config["bacterium"]["h"]),
    )

    velocity_frames = []
    bacterium_path  = []
    tail_coords     = []

    if abs(OMEGA_BUNDLE) > 0.0:
        steps_per_rot = int(round((2.0*np.pi) / (abs(OMEGA_BUNDLE) * dt)))
        print(f"~{steps_per_rot} steps per bundle rotation (dt={dt})")

    for step in range(steps):
        rho, u = macroscopic(f)
        t = step * dt

        bact.update_geometry(t)

        # Build forces (helical bundle + body ring + optional tail + optional stresslet)
        F = build_forces(bact, u.shape, t)

        # LBM collide/stream with Guo forcing
        f = collision_with_guo(f, rho, u, tau, F, dt)
        f = streaming(f)

        # Bacterium kinematics
        bact.advect(u, dt)

        if step % 50 == 0:
            umax = float(np.max(np.linalg.norm(u, axis=-1)))
            print(f"[{step:05d}] |u|max={umax:.3e}  center={bact.center}")

        velocity_frames.append(u.copy())
        bacterium_path.append(bact.center.copy())
        tail_coords.append(np.vstack([bact.tail_x, bact.tail_y, bact.tail_z]).T)

    out_dir = os.path.abspath(
        os.path.join(os.path.dirname(__file__), "..", config["output"]["path"])
    )
    os.makedirs(out_dir, exist_ok=True)

    np.save(os.path.join(out_dir, config["output"]["velocity_field"]),
            np.array(velocity_frames))
    np.save(os.path.join(out_dir, config["output"]["velocity_field_small"]),
            np.array(velocity_frames)[::10])
    np.save(os.path.join(out_dir, config["output"]["bacterium_path"]),
            np.array(bacterium_path))
    np.save(os.path.join(out_dir, config["output"]["tail_coords"]),
            np.array(tail_coords))

    meta = np.array([OMEGA_BUNDLE, OMEGA_BODY, dt], dtype=float)
    np.save(os.path.join(out_dir, "rotation_meta.npy"), meta)

    print(f"Saved: {out_dir}")

if __name__ == "__main__":
    simulate_hybrid_lbm_tiger()
