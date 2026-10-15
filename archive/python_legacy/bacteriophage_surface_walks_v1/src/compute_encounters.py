def compute_encounters(self, res_times, dep_vels, time):
    """
    Detects phage-bacteria encounters and updates particle state accordingly.
    Implements logic from encounter_fun_hydro.m.
    """
    MatP = self.MatP
    MatB = self.MatB
    dt = self.params["dt"]
    Rp = self.params["Rp"]
    Rcrit = self.params["Rb"] + Rp
    D_para_y = self.params["D_para_y"]
    D_para_z = self.params["D_para_z"]
    LX = self.params["LX"]

    # Count active phages
    Nact = np.sum(MatP[:, 5])

    # Distance from bacterium
    dist = np.sqrt(
        MatP[:, 0] ** 2 +
        D_para_y * MatP[:, 1] ** 2 +
        D_para_z * MatP[:, 2] ** 2
    )

    # Update residence time for phages within Rcrit + Rp
    near_idx = np.where(dist <= Rcrit + Rp)[0]
    MatP[near_idx, 5] *= 1  # make sure activity is boolean
    MatP[near_idx, 6] += dt * MatP[near_idx, 5]
    self.N_dep = 100 * len(near_idx) / 10000

    # Count conditional encounters (active particles within Rcrit)
    enc_idx = np.where(dist <= Rcrit)[0]
    Nenc = int(np.sum(MatP[enc_idx, 5]))

    if len(enc_idx) > 0 and time >= 0:
        active_idx = enc_idx[MatP[enc_idx, 5] == 1]

        if len(active_idx) > 0:
            # Compute deposition velocity
            dist_v = np.column_stack([
                MatP[active_idx, 7],
                D_para_y * MatP[active_idx, 8],
                D_para_z * MatP[active_idx, 9]
            ])
            vel_v = np.column_stack([
                (MatP[active_idx, 0] - MatP[active_idx, 7]) / dt,
                D_para_y * (MatP[active_idx, 1] - MatP[active_idx, 8]) / dt,
                D_para_z * (MatP[active_idx, 2] - MatP[active_idx, 9]) / dt
            ])

            dot = np.sum(vel_v * dist_v, axis=1)
            proj = dot / np.sqrt(np.sum(dist_v ** 2, axis=1))
            # dep_vels.extend(-proj.tolist())  # Optional storage of velocity

            # Event-driven collision correction
            vel = MatP[enc_idx, 0:3] - MatP[enc_idx, 7:10]
            r_old = MatP[enc_idx, 7:10]
            v_mag_sq = np.sum(vel ** 2, axis=1)

            b = 2 * np.sum(r_old * vel, axis=1) / v_mag_sq
            c = (np.sum(r_old ** 2, axis=1) - Rcrit ** 2) / v_mag_sq
            tau = -(b + np.sqrt(b ** 2 - 4 * c)) / 2

            rc = vel * tau[:, None] + r_old
            q = MatP[enc_idx, 0:3] - rc
            perp = -rc * np.sum(rc * q, axis=1)[:, None] / np.sum(rc ** 2, axis=1)[:, None]
            MatP[enc_idx, 0:3] = rc + q + 2 * perp

    # Reset state for deactivated or distant particles
    MatP[enc_idx, 5] = 0  # deactivate
    MatP[enc_idx, 6] = 0  # reset res. times
    MatP[dist > Rcrit + Rp, 6] = 0

    return Nenc, Nact, res_times, dep_vels, self.N_dep
