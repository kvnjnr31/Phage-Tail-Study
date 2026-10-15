
import numpy as np

class PhageSimulator:
    def __init__(self, B, P, params):
        self.B = B
        self.P = P
        self.params = params
        self.Time = 0
        self.n = 1
        self.N_dep = 0
        self.kernel_new = np.zeros((1, params["nsteps"]))
        self.kernel_avg = np.zeros((1, params["nsteps"]))
        self._initialize_system()

    def _initialize_system(self):
        # Initialize MatB and MatP
        pass

    def update_phage_positions(self, dt):
        """ Langevin_fun equivalent """
        # TODO: Translate Langevin_fun logic here
        pass

    
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

    
def resolve_overlap(self):
    """
    Check for and optionally respond to phage-bacterium overlaps.

    This method translates overlap.m. In its original form, it visualized overlapping
    phages. Here, we use it to optionally log or flag overlaps, or extend it in the future
    to apply collision logic.
    """
    MatP = self.MatP
    MatB = self.MatB
    Rb = self.params["Rb"]
    D_para_y = self.params["D_para_y"]
    D_para_z = self.params["D_para_z"]

    for p in range(MatP.shape[0]):
        r = np.sqrt(
            (MatP[p, 0] - MatB[0, 0])**2 +
            D_para_y * (MatP[p, 1] - MatB[0, 1])**2 +
            D_para_z * (MatP[p, 2] - MatB[0, 2])**2
        )
        if r < Rb:
            # Placeholder: add logging, counter, or corrective action here
            # For now, we simply acknowledge overlap was detected
            pass  # TODO: implement actual handling logic if desired

    def compute_encounters(self, res_times, dep_vels, time):
        """ encounter_fun_hydro equivalent """
        # TODO: Translate encounter_fun_hydro logic here
        Nenc, Nact = 0, 1  # placeholders
        return Nenc, Nact, res_times, dep_vels, self.N_dep

    def update_bacteria_direction(self):
        # Generate new random run direction and rotate phage positions
        pass

    def apply_boundary_conditions(self):
        # Apply hydro-style boundary conditions
        pass

    def plot_simulation_state(self):
        # Final plot logic goes here
        pass



def update_bacteria_direction(self):
    """
    Select a new bacteria run direction and rotate the phage positions accordingly.
    This models a tumbling event and changes the phage swarm direction.
    """
    MatP = self.MatP
    omega = self.params["omega"]
    D_para_y = self.params["D_para_y"]
    D_para_z = self.params["D_para_z"]

    # Choose a small angle of change (currently fixed at 0 degrees)
    angle_deg = 0 * 3.5
    angle_rad = angle_deg * (2 * np.pi / 360)

    # Generate constrained random unit vector
    vx = (np.cos(angle_rad) - 1) * np.random.rand() + 1
    vz = D_para_z * (2 * np.sqrt(1 - vx**2) * np.random.rand() - np.sqrt(1 - vx**2))
    b = 2 * np.random.rand() - 1
    vy = np.sign(b) * np.sqrt(1 - vx**2 - vz**2)

    # Update next run time
    tau = -np.log(np.random.rand()) / omega
    self.run_time = self.Time + tau

    # Rotation angles
    theta = -np.arctan2(vy, vx)  # rotation around z
    alpha = np.arcsin(vz)        # rotation around y

    # Rotation about z-axis
    x_new = MatP[:, 0] * np.cos(theta) - MatP[:, 1] * np.sin(theta)
    MatP[:, 1] = MatP[:, 0] * np.sin(theta) + MatP[:, 1] * np.cos(theta)

    # Rotation about y-axis
    MatP[:, 0] = x_new * np.cos(alpha) + MatP[:, 2] * np.sin(alpha)
    MatP[:, 2] = -x_new * np.sin(alpha) + MatP[:, 2] * np.cos(alpha)



def plot_simulation_state(self):
    """
    Plot the current phage concentration and spatial distribution.
    Includes:
    - Radial concentration histogram (normalized)
    - 3D scatter of rotated phage positions
    """
    import matplotlib.pyplot as plt

    MatP = self.MatP
    MatB = self.MatB
    Rb = self.params["Rb"]
    Rp = self.params["Rp"]
    LX = self.params["LX"]
    LY = self.params.get("LY", LX)
    Time = self.Time
    N_dep = self.N_dep

    # Rotate phage positions onto x-y plane for visualization
    MatR = np.zeros_like(MatP[:, 0:3])
    for p in range(MatP.shape[0]):
        v = MatP[p, 0:3]
        theta = -np.arctan2(v[2], v[1])
        R = np.array([
            [1, 0, 0],
            [0, np.cos(theta), -np.sin(theta)],
            [0, np.sin(theta),  np.cos(theta)]
        ])
        MatR[p, :] = R @ v

    r = np.linalg.norm(MatP[:, 0:3], axis=1)
    r = r[r <= 2.5]
    n, edges = np.histogram(r, bins=40)
    radius = 0.5 * (edges[1:] + edges[:-1])
    dr = (radius[1] - radius[0]) / 2
    volumes = (4/3) * np.pi * ((radius + dr)**3 - (radius - dr)**3)

    plt.figure(1)
    plt.plot(radius, n / (np.sum(n) * volumes), 'r', linewidth=2)
    plt.xlabel('Radius')
    plt.ylabel('Phage Concentration')
    plt.axis('square')
    plt.grid(True)
    plt.gca().set(title='Radial Distribution', fontsize=14, linewidth=2)

    plt.figure(2)
    ax = plt.axes(projection='3d')
    ax.scatter(MatR[:, 0], MatR[:, 1], MatR[:, 2], c='k', s=2)
    ax.set_xlim([-LX/2, LX/2])
    ax.set_ylim([-LY/2, LY/2])
    ax.set_zlim([-LY/2, LY/2])
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_title('Phage Swarm View')
    ax.view_init(45, 45)
    plt.tight_layout()
    plt.show()
