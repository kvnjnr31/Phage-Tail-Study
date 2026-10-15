def apply_boundary_conditions(self):
    """
    Apply reflective boundary conditions to phage at edges of the simulation domain.
    Implements logic from boundarys_hydro.m.
    """
    MatP = self.MatP
    Rp = self.params["Rp"]
    x_right = self.params["x_right"]
    y_right = self.params["y_right"]
    z_right = self.params["z_right"]
    D_para_y = self.params["D_para_y"]
    D_para_z = self.params["D_para_z"]

    # X-boundary
    jj = np.where(np.abs(MatP[:, 0]) >= x_right)[0]
    MatP[jj, 0] = x_right - Rp
    MatP[jj, 1] = D_para_y * (2 * y_right * np.random.rand(len(jj)) - y_right)
    MatP[jj, 2] = D_para_z * (2 * z_right * np.random.rand(len(jj)) - z_right)

    # Record old and initial positions
    MatP[jj, 7:10] = MatP[jj, 0:3]
    MatP[jj, 10:13] = MatP[jj, 0:3]
    MatP[jj, 5] = 1  # activate

    # Y-boundary
    jj = np.where(np.abs(MatP[:, 1]) >= y_right)[0]
    MatP[jj, 1] = D_para_y * (MatP[jj, 1] - 2 * y_right * MatP[jj, 1] / np.abs(MatP[jj, 1]))
    MatP[jj, 2] = D_para_z * MatP[jj, 2]

    MatP[jj, 7:10] = MatP[jj, 0:3]
    MatP[jj, 10:13] = MatP[jj, 0:3]
    MatP[jj, 5] = 1  # activate

    # Z-boundary
    jj = np.where(np.abs(MatP[:, 2]) >= z_right)[0]
    MatP[jj, 2] = D_para_z * (MatP[jj, 2] - 2 * z_right * MatP[jj, 2] / np.abs(MatP[jj, 2]))

    MatP[jj, 7:10] = MatP[jj, 0:3]
    MatP[jj, 10:13] = MatP[jj, 0:3]
    MatP[jj, 5] = 1  # activate
