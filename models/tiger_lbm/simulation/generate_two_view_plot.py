import numpy as np
import matplotlib.pyplot as plt
import os

# === CONFIGURATION ===
data_path = "../data/numerical/velocity_field_small.npy"
output_path = "../data/image/quiver_velocity_views.png"

# Choose slicing planes
y_index = 50  # For X-Z slice (Right View)
x_index = 30  # For Y-Z slice (Back View)

# Subsampling for clarity in quiver plot
step = 2  # increase for fewer arrows

# === LOAD FIELD ===
velocity = np.load(data_path)
velocity = np.squeeze(velocity)  # Remove leading singleton dimension if present
assert velocity.ndim == 5 and velocity.shape[-1] == 3, \
    f"Expected shape (nx, ny, nz, 3), got {velocity.shape}"

# === EXTRACT VIEWS ===

# RIGHT VIEW: X-Z Plane @ y_index
u_r = velocity[::step, y_index, ::step, 0]  # x-component
w_r = velocity[::step, y_index, ::step, 2]  # z-component
X_r, Z_r = np.meshgrid(
    np.arange(0, velocity.shape[0], step),
    np.arange(0, velocity.shape[2], step),
    indexing='ij'
)

# BACK VIEW: Y-Z Plane @ x_index
v_b = velocity[x_index, ::step, ::step, 1]  # y-component
w_b = velocity[x_index, ::step, ::step, 2]  # z-component
Y_b, Z_b = np.meshgrid(
    np.arange(0, velocity.shape[1], step),
    np.arange(0, velocity.shape[2], step),
    indexing='ij'
)

# === PLOT ===
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
fig.suptitle("TIGER Simulation – Velocity Field Quiver Views", fontsize=16)

# RIGHT VIEW
ax0 = axes[0]
ax0.quiver(X_r, Z_r, u_r, w_r, scale=5, color='white')
ax0.set_title(f"Right View (X–Z plane) @ y={y_index}")
ax0.set_xlabel("X")
ax0.set_ylabel("Z")
ax0.set_aspect('equal')
ax0.set_facecolor('black')

# BACK VIEW
ax1 = axes[1]
ax1.quiver(Y_b, Z_b, v_b, w_b, scale=5, color='white')
ax1.set_title(f"Back View (Y–Z plane) @ x={x_index}")
ax1.set_xlabel("Y")
ax1.set_ylabel("Z")
ax1.set_aspect('equal')
ax1.set_facecolor('black')

plt.tight_layout(rect=[0, 0.03, 1, 0.95])
os.makedirs(os.path.dirname(output_path), exist_ok=True)
plt.savefig(output_path, dpi=300)
plt.close()

print(f"✅ Quiver plots saved to: {output_path}")
