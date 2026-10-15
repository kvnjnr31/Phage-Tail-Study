# Analytic TWIN

## Overview

This repository contains the core framework for the **Analytic TWIN** — a physics-guided digital twin architecture designed to replicate and analyze microscale bacterial-phage interactions under hydrodynamic constraints.

With the emergence of the **Analytic TWIN**, we now define a lightweight bacterial agent that acts as a surrogate for the high-resolution TIGER simulation. This twin enables fast, scalable modeling in large systems, marking the birth of our first **bacterial agent**.

## Equation of Motion

The bacterial agent evolves according to a simplified form of Newtonian motion under drag:
\[
\frac{d\mathbf{v}}{dt} = -\gamma (\mathbf{v} - \mathbf{u}_{\text{fluid}})
\]
where:
- \( \mathbf{v} \) is the velocity of the bacterium
- \( \mathbf{u}_{\text{fluid}} \) is the local fluid velocity (from simulation or model)
- \( \gamma \) is the damping coefficient (fluid resistance)

This can be further recast into an **Euler–Lagrange** formulation in the upcoming module, laying the groundwork for general variational modeling, optimal control, and physics-informed training of surrogate models.

## Folder Structure

```
.
├── data/                 # Stores output files: velocity fields, paths, and tail coordinates
├── simulation/           # Main hybrid TIGER-LBM simulation scripts
├── utils/                # Shared utilities for plotting and analysis
├── figures/              # Rendered plots and generated visuals
├── render/               # Output PNGs and movies for visualization
├── README.md             # This file
```

## Agent-Based Extension

The Analytic TWIN now forms the basis for **bacterial agents** that are:
- Lightweight and scalable
- Analytically defined (ODE-based)
- Easily embeddable in large phage-bacteria population simulations

These agents evolve based on far-field approximations of drag, rotation, and net fluid coupling — enabling mass simulation without costly LBM computation.

## Key Components

- `hybrid_simulation.py`: Core simulation loop using TIGER bacterium and LBM fluid coupling.
- `generate_frames.py`: 2D velocity slice visualizer.
- `make_phage_movie.py`: Renders 3D bacterium tail motion over time.
- `plot_velocity_field.py`: Utilities for quiver plots and field rendering.

## Dependencies

- Python 3.8+
- NumPy
- Matplotlib
- ImageIO

## Getting Started

1. Run the main simulation:
   ```bash
   python simulation/hybrid_simulation.py
   ```

2. Generate 2D velocity slices:
   ```bash
   python simulation/generate_frames.py
   ```

3. Create a 3D tail-motion movie:
   ```bash
   python simulation/make_phage_movie.py render
   ```

## Author

**Kevin L. Joiner**  
AI Low Voltage Innovations LLC  
Veteran-Owned | Research-Led

---

> _This README will be updated as the Analytic TWIN architecture expands to include deterministic-stochastic coupling, boundary confinement models, and embedded optimization routines._


## Architectural Interoperability

To ensure robust simulation and collaboration across biological agents, the folder structure of this repository is intentionally designed to mirror its counterpart in `bacteriophage_ABM/`. This parallel architecture supports:

- ✅ **Modular Communication**: Both phage and bacterium agents can share protocols, interfaces, and environmental states.
- ✅ **Plug-and-Play Components**: Shared code between agents, such as the fluid model or control systems, can be reused cleanly.
- ✅ **Scalable Integration**: Simulation engines can run thousands of phages and bacteria concurrently with matching logic trees.
- ✅ **Agent Co-Evolution**: Facilitates a shared simulation environment where each agent (bacteria or phage) can evolve interactively.

This design ensures that the `bacteria_ABM/` and `bacteriophage_ABM/` directories serve as aligned, high-performance units capable of supporting emergent behaviors, mutual interactions, and future reinforcement learning or variational optimization pipelines.
