def update_bacteria_direction(self):
    """
    Update phage direction by rotating their positions in the bacterium-centered frame.
    The bacterium remains fixed; this simulates a change in bacterial swimming direction.
    """
    MatP = self.MatP
    omega = self.params["omega"]
    D_para_y = self.params["D_para_y"]
    D_para_z = self.params["D_para_z"]

    # --- Generate new random direction vector in spherical shell ---
    angle_deg = 0 * 3.5  # Optional small deflection angle
    angle_rad = angle_deg * (2 * np.pi / 360)

    vx = (np.cos(angle_rad) - 1) * np.random.rand() + 1
    vz = D_para_z * (2 * np.sqrt(1 - vx**2) * np.random.rand() - np.sqrt(1 - vx**2))
    b = 2 * np.random.rand() - 1
    vy = np.sign(b) * np.sqrt(1 - vx**2 - vz**2)

    # --- Update next tumble time ---
    tau = -np.log(np.random.rand()) / omega
    self.run_time = self.Time + tau

    # --- Compute rotation angles ---
    theta = -np.arctan2(vy, vx)  # rotate in x-y plane (z-axis rotation)
    alpha = np.arcsin(vz)        # rotate in x-z plane (y-axis rotation)

    # --- Rotation: about z-axis ---
    x = MatP[:, 0]
    y = MatP[:, 1]
    z = MatP[:, 2]

    x_new = x * np.cos(theta) - y * np.sin(theta)
    y_new = x * np.sin(theta) + y * np.cos(theta)
    z_new = z

    # --- Rotation: about y-axis (applied after z-rotation) ---
    x_rot = x_new * np.cos(alpha) + z_new * np.sin(alpha)
    y_rot = y_new
    z_rot = -x_new * np.sin(alpha) + z_new * np.cos(alpha)

    # --- Update phage positions (new direction) ---
    MatP[:, 0] = x_rot
    MatP[:, 1] = y_rot
    MatP[:, 2] = z_rot
