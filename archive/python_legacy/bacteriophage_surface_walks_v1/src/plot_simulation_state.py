def plot_simulation_state(self):
    """
    Generate simulation plots:
    1. Radial phage concentration distribution
    2. 3D phage scatter after planar rotation
    """
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    MatP = self.MatP
    Rb = self.params["Rb"]
    Rp = self.params["Rp"]
    LX = self.params["LX"]
    LY = self.params.get("LY", LX)
    Time = self.Time
    N_dep = self.N_dep

    # --- Rotate phage positions onto x-y plane ---
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

    # --- Histogram of phage radial distribution ---
    r = np.linalg.norm(MatP[:, :3], axis=1)
    r = r[r <= 2.5]  # Remove outliers

    n, edges = np.histogram(r, bins=40)
    radius = 0.5 * (edges[1:] + edges[:-1])
    dr = (radius[1] - radius[0]) / 2
    volumes = (4/3) * np.pi * ((radius + dr)**3 - (radius - dr)**3)

    plt.figure(1)
    plt.plot(radius, n / (np.sum(n) * volumes), 'r', linewidth=2)
    plt.xlabel('Radius')
    plt.ylabel('Phage Concentration')
    plt.title('Radial Distribution')
    plt.axis('square')
    plt.grid(True)
    plt.gca().tick_params(labelsize=14, width=2)
    plt.tight_layout()

    # --- 3D Scatter of rotated phage cloud ---
    plt.figure(2)
    ax = plt.axes(projection='3d')
    ax.scatter(MatR[:, 0], MatR[:, 1], MatR[:, 2], c='k', s=2)
    ax.set_xlim([-LX / 2, LX / 2])
    ax.set_ylim([-LY / 2, LY / 2])
    ax.set_zlim([-LY / 2, LY / 2])
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_title('Phage Swarm View')
    ax.view_init(45, 45)
    ax.set_box_aspect([1, 1, 1])
    plt.tight_layout()
    plt.show()
