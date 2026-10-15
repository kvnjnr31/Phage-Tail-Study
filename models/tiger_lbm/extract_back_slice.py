#!/usr/bin/env python3
import os
import numpy as np

# Paths to your data (match your screenshot)
VEL_PATH = os.path.join("data", "numerical", "velocity_field_small.npy")
OUT_PATH = os.path.join("data", "numerical", "back_slice.npy")

# For this run, the bacterium starts at [48, 32, 32] with Rb = 10
# Helical bundle is centered around x ≈ center_x + 0.6 * Rb
center_x = 48.0
Rb = 10.0
x_index = int(round(center_x + 0.6 * Rb))  # = 54

print(f"Using x_index = {x_index}")

# Load the 3D velocity field
vel = np.load(VEL_PATH)          # shape (T, nx, ny, nz, 3)
T, nx, ny, nz, _ = vel.shape
print(f"Loaded velocity_field_small: {vel.shape}")

if not (0 <= x_index < nx):
    raise ValueError(f"x_index {x_index} is out of bounds for nx={nx}")

# Extract rear-view Y–Z slice for all times:
# result shape: (T, ny, nz, 3)
back_slice = vel[:, x_index, :, :, :]

# Save it
os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
np.save(OUT_PATH, back_slice)
print(f"Saved back-view slice to {OUT_PATH} with shape {back_slice.shape}")
