def generate_noise(self, dt, tmin, nu, time):
    """
    Apply Brownian motion noise to phages that are ready to jump.

    This method replicates noise_gen.m by finding phages whose transition timers have expired
    and assigning them new stochastic positions while respecting minimum proximity to the central bacterium.
    """
    MatP = self.MatP
    MatB = self.MatB
    D0 = self.params["D0"]
    Rb = self.params["Rb"]
    Rp = self.params["Rp"]
    D_para_y = self.params["D_para_y"]
    D_para_z = self.params["D_para_z"]

    # Find phages whose waiting time has elapsed
    ready_idx = np.where(time - MatP[:, 3] >= 0)[0]

    # Generate new waiting times (power-law distributed)
    U1 = np.random.rand(len(ready_idx))
    tau = tmin * (1.0 / U1)**(1.0 / nu)
    MatP[ready_idx, 3] = time + tau  # Update transition timers

    # Apply Brownian motion with spatial rejection to avoid overlap
    for i, p in enumerate(ready_idx):
        r = 0
        while r < (Rb + Rp):
            # Generate Brownian displacements
            wx = np.sqrt(24 * D0 * dt) * (np.random.rand() - 0.5)
            wy = np.sqrt(24 * D0 * dt) * (np.random.rand() - 0.5)
            wz = np.sqrt(24 * D0 * dt) * (np.random.rand() - 0.5)

            # Proposed new positions
            x_new = MatP[p, 0] + wx
            y_new = D_para_y * (MatP[p, 1] + wy)
            z_new = D_para_z * (MatP[p, 2] + wz)

            # Distance to central bacterium
            r = np.sqrt((x_new - MatB[0, 0])**2 +
                        D_para_y * (y_new - MatB[0, 1])**2 +
                        D_para_z * (z_new - MatB[0, 2])**2)

        # Accept displacement if outside exclusion zone
        MatP[p, 0] = x_new
        MatP[p, 1] = y_new
        MatP[p, 2] = z_new
