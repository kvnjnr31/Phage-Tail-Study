def compute_vorticity_z(u):
    """
    Compute the Z-component of the vorticity vector field from velocity field u.
    u should be a (nx, ny, nz, 3) array representing (ux, uy, uz).
    """
    ux = u[..., 0]
    uy = u[..., 1]

    dudy = np.gradient(uy, axis=1)
    dudx = np.gradient(ux, axis=0)

    omega_z = dudy - dudx  # ∂uy/∂x - ∂ux/∂y
    return omega_z


def plot_vorticity_frame(u, data/image/frame_idx=0, output_path="vorticity_framedata/image/%s.png", bacterium_center=None):
    """
    Visualize the Z-component of the vorticity field in the X-Z plane (middle Y slice).
    """
    omega_z = compute_vorticity_z(u)

    nx, ny, nz = omega_z.shape
    mid_y = ny // 2
    vorticity_slice = omega_z[:, mid_y, :]

    plt.figure(figsize=(8, 6))
    plt.imshow(vorticity_slice.T, origin="lower", cmap="seismic", aspect="auto")
    plt.colorbar(label="ω_z (Vorticity)")
    plt.title(f"Z-Vorticity in XZ Plane @ y={mid_y} (Frame {data/image/frame_idx})")
    if bacterium_center is not None:
        plt.plot(bacterium_center[0], bacterium_center[2], 'ko', markersize=5)
    plt.xlabel("X")
    plt.ylabel("Z")
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
