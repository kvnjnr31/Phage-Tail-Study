
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation

# Parameters
h = 0.77
kE = 1.3
k = 1.3
w = -100
Rb = 1
dt = 0.05
t_total = 5
frames = int(t_total / dt)

# Spatial grid for flagellum
x = np.linspace(0, 7, 25)
Ex = 1 - np.exp(-(kE**2) * x**2)

# Head coordinates
head_coord_x0 = np.linspace(-2*(3*Rb), 0, 10)
head_coord_x = list(head_coord_x0)
hy = np.zeros_like(head_coord_x0)
hz = Rb * np.sqrt(1 - ((head_coord_x0 + (3 * Rb)) / (3 * Rb))**2)
head_theta0 = np.vstack((hy, hz))
head_coord_yz = head_theta0.copy()

# Arc expansion
theta = np.pi / 3
R_theta = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
for _ in range(6):
    head_theta0 = R_theta @ head_theta0
    head_coord_yz = np.hstack((head_coord_yz, head_theta0))
    head_coord_x = np.concatenate((head_coord_x, head_coord_x0))

# Rotation over time
theta_step = -np.pi / 24
R_theta_time = np.array([[np.cos(theta_step), -np.sin(theta_step)],
                         [np.sin(theta_step), np.cos(theta_step)]])

# Set up figure
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

def update(data/image/frame_num):
    ax.clear()
    global head_coord_yz
    head_coord_yz = R_theta_time @ head_coord_yz
    t = data/image/frame_num * dt
    y = h * Ex * np.sin(k * x - w * t)
    z = h * Ex * np.cos(k * x - w * t)

    ax.plot3D(head_coord_x, head_coord_yz[0], head_coord_yz[1], 'ko--', markersize=4)
    ax.plot3D(x, y, z, 'bo--', markersize=4)

    ax.set_xlim([-6, 7])
    ax.set_ylim([-2, 2])
    ax.set_zlim([-2, 2])
    ax.set_title(f"TIGER Time: {t:.2f} sec")
    ax.set_xlabel("L_x")
    ax.set_ylabel("L_y")
    ax.set_zlabel("L_z")
    ax.view_init(elev=10, azim=-22)

ani = animation.FuncAnimation(fig, update, frames=frames, interval=50)

# Save the animation (optional)
ani.save("TIGER_simulationdata/image/%s.mp4", writer='ffmpeg')
