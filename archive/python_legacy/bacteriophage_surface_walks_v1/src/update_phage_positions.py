def update_phage_positions(self, dt):
    """
    Langevin_fun: Update phage positions using Stokes flow and bacteria velocity field.
    """
    MatP = self.MatP
    MatB = self.MatB
    Rb = self.params["Rb"]
    D_para_y = self.params["D_para_y"]
    D_para_z = self.params["D_para_z"]

    # Compute radial distance and higher powers
    r = np.sqrt(MatP[:, 0]**2 + D_para_y * MatP[:, 1]**2 + D_para_z * MatP[:, 2]**2)
    r2 = r**2
    r3 = r**3
    r4 = r**4

    # Dot product of bacteria velocity with phage position
    dot = (MatB[0, 3] * MatP[:, 0] +
           D_para_y * MatB[0, 4] * MatP[:, 1] +
           D_para_z * MatB[0, 5] * MatP[:, 2])

    # Velocity components from Stokes flow
    v_comp_x = (3 * Rb / (4 * r)) * (1 + (Rb**2) / (3 * r2)) * MatB[0, 3]
    v_comp_y = D_para_y * (3 * Rb / (4 * r)) * (1 + (Rb**2) / (3 * r2)) * MatB[0, 4]
    v_comp_z = D_para_z * (3 * Rb / (4 * r)) * (1 + (Rb**2) / (3 * r2)) * MatB[0, 5]

    # Non-radial correction (proportional to dot product)
    v_corr = (3 * Rb**3 / (4 * r4)) * dot
    v_corr_x = v_corr * MatP[:, 0]
    v_corr_y = v_corr * D_para_y * MatP[:, 1]
    v_corr_z = v_corr * D_para_z * MatP[:, 2]

    # Combine fields
    vx = v_comp_x + v_corr_x
    vy = v_comp_y + v_corr_y
    vz = v_comp_z + v_corr_z

    # Update phage positions
    self.MatP[:, 0] += dt * vx
    self.MatP[:, 1] += dt * vy
    self.MatP[:, 2] += dt * vz
