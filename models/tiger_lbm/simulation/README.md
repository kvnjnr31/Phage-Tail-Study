# TIGER: Tail Integrated Grid Embedded Rigid-body Model

**TIGER** models bacteria swimming using explicit tail-fluid coupling:

- The **bacterium head** is a rigid ellipsoid with time-evolving position $ \vec{x}(t) $ and orientation angle $ \theta(t) $.
- The **flagellar tail** is a semi-flexible arc discretized into multiple segments, oscillating in a plane using:

$$
x_i = -2Rb + i \cdot \Delta s,\quad y_i = h \cdot \sqrt{1 - \left(\frac{x_i}{2Rb}\right)^2}
$$

The tail arc is rotated using a time-dependent angle $ \theta(t) = \theta_0 + \omega t $, and projected into the fluid frame using rotation matrices:

$$
R(\theta) =
\begin{bmatrix}
\cos\theta & -\sin\theta \\
\sin\theta & \cos\theta
\end{bmatrix}
$$

At each timestep, momentum is injected into nearby fluid nodes from the tail using a projection kernel, simulating swimming.

## Lattice Boltzmann Method (LBM)

We solve mesoscale fluid dynamics using the Lattice Boltzmann method on a D3Q19 lattice. The evolution of the single-particle distribution function \( f(\vec{x}, \vec{\xi}, t) \) obeys a discretized Boltzmann equation:

$$
f_\alpha(\vec{x} + \vec{c}_\alpha \Delta t, t + \Delta t) = f_\alpha(\vec{x}, t) - \frac{\Delta t}{\tau} \left[f_\alpha(\vec{x}, t) - f_\alpha^{eq}(\vec{x}, t)\right]
$$

- $ f_\alpha $: particle distribution in direction $ \alpha $
- $ \tau $: relaxation time controlling viscosity
- $ f_\alpha^{eq} $: Maxwell-Boltzmann equilibrium distribution
- $ \vec{c}_\alpha $: discrete lattice velocities (D3Q19)
- Boundary conditions: bounce-back (no-slip) and periodic

The macroscopic velocity field $ \vec{u} $ and density $ \rho $ are computed from the zeroth and first moments of $f$.



## Output Data

- `velocity_field.npy` — fluid velocity $ \vec{u}(\vec{x}, t) $
- `bacterium_path.npy` — head position over time
- `tail_coords.npy` — 3D coordinates of flagellar tail
- `fluid_bacterium_movie.mp4` — composite velocity field render

