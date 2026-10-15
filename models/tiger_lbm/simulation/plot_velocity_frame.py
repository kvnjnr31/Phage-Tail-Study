import numpy as np
import matplotlib.pyplot as plt
import os

def plot_velocity_frame(u, bacterium, t, output_dir="frames", dpi=300,
                        bacterium_velocity=None, tail_points=None, debug=True):
    """
    Web-optimized 2x2 velocity panel from TIGER simulation.
    Includes: large fonts, timestamp, darker theme, and clearer visual cues.
    """
    os.makedirs(output_dir, exist_ok=True)

    if bacterium_velocity is not None:
        u = u - np.array(bacterium_velocity)

    fig, axs = plt.subplots(2, 2, figsize=(12, 9))  # larger figure
    fig.suptitle("TIGER Simulation – Velocity Field Cross-Sections", fontsize=22, weight='bold')

    views = ["Left", "Right", "Front", "Back"]
    slice_offsets = [-5, 5, -5, 5]
    bacterium = np.array(bacterium).astype(int)
    pos_map = [(0, 0), (0, 1), (1, 0), (1, 1)]

    for idx, (view, offset, pos) in enumerate(zip(views, slice_offsets, pos_map)):
        ax = axs[pos]

        if view in ["Front", "Back"]:
            slice_x = np.clip(bacterium[0] + offset, 0, u.shape[0] - 1)
            uy = u[slice_x, :, :, 1]
            uz = u[slice_x, :, :, 2]
            X, Y = np.meshgrid(np.arange(uy.shape[0]), np.arange(uy.shape[1]), indexing="xy")
            speed = np.sqrt(uy**2 + uz**2)
            im = ax.imshow(speed.T, origin="lower", cmap="plasma", alpha=0.85)
            ax.quiver(X, Y, uy.T, uz.T, scale=30, color='white', linewidth=0.4)
            ax.set_title(f"{view} View @ x={slice_x}", fontsize=16)

            if tail_points is not None:
                tail = np.array(tail_points).astype(int)
                mask = tail[:, 0] == slice_x
                ax.plot(tail[mask][:, 1], tail[mask][:, 2], 'go', markersize=3)

            if debug:
                ax.plot(bacterium[1], bacterium[2], 'rx', markersize=6)

        elif view in ["Left", "Right"]:
            slice_y = np.clip(bacterium[1] + offset, 0, u.shape[1] - 1)
            ux = u[:, slice_y, :, 0]
            uz = u[:, slice_y, :, 2]
            X, Y = np.meshgrid(np.arange(ux.shape[0]), np.arange(ux.shape[1]), indexing="xy")
            speed = np.sqrt(ux**2 + uz**2)
            im = ax.imshow(speed.T, origin="lower", cmap="plasma", alpha=0.85)
            ax.quiver(X, Y, ux.T, uz.T, scale=30, color='white', linewidth=0.4)
            ax.set_title(f"{view} View @ y={slice_y}", fontsize=16)

            if tail_points is not None:
                tail = np.array(tail_points).astype(int)
                mask = tail[:, 1] == slice_y
                ax.plot(tail[mask][:, 0], tail[mask][:, 2], 'go', markersize=3)

            if debug:
                ax.plot(bacterium[0], bacterium[2], 'rx', markersize=6)

        # === Styling ===
        ax.set_facecolor('black')
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_aspect('equal')

    # === Timestamp in corner ===
    fig.text(0.82, 0.92, f"t = {t * 0.05:.2f}s", fontsize=14, color='gray')

    # === Optional colorbar ===
    cbar_ax = fig.add_axes([0.92, 0.25, 0.015, 0.5])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.set_label("Velocity Magnitude", fontsize=12)
    cbar.ax.tick_params(labelsize=10)

    fig.patch.set_facecolor('black')
    plt.tight_layout(rect=[0, 0, 0.9, 0.95])
    save_path = os.path.join(output_dir, f"frame_{t:04d}.png")
    plt.savefig(save_path, dpi=dpi)
    plt.close()
