# Analytic TWIN

## Overview

This repository contains the core framework for the **Analytic TWIN** — a physics-guided digital twin architecture designed to replicate and analyze microscale bacterial-phage interactions under hydrodynamic constraints.

With the emergence of the **Analytic TWIN**, we now define a lightweight bacterial agent that acts as a surrogate for the high-resolution TIGER simulation. This twin enables fast, scalable modeling in large systems, marking the birth of our first **bacterial agent**.

## Equation of Motion

The bacterial agent evolves according to a simplified form of Newtonian motion under drag:
$$
\frac{d\mathbf{v}}{dt} = -\gamma (\mathbf{v} - \mathbf{u}_{\text{fluid}})
$$
where:
- $\mathbf{v}$ is the velocity of the bacterium
- $\mathbf{u}_{\text{fluid}}$ is the local fluid velocity (from simulation or model)
- $\gamma$ is the damping coefficient (fluid resistance)

This can be further recast into an **Euler–Lagrange** formulation in the upcoming module, laying the groundwork for general variational modeling, optimal control, and physics-informed training of surrogate models.


## Agent-Based Extension

The Analytic TWIN forms the basis for **bacterial agents** that are:
- Lightweight and scalable
- Analytically defined (ODE-based)
- Easily embeddable in large phage-bacteria population simulations

These agents evolve based on far-field approximations of drag, rotation, and net fluid coupling — enabling mass simulation without costly LBM computation.

## Key Components

- `hybrid_simulation.py`: Core simulation loop using TIGER bacterium and LBM fluid coupling.
- `phage_bulk_surface_exchange.py`: 3D random walks of phage in bulk and on the bacterium surface.
- `phage_surface_walk.py`: 3D random walks of phage on the bacterium surface only.
- `validate_phage_surface_walk.py`: 3D random walks of phage in bulk and on the bacterium surface.
- `PhageSimpleTail.py`: a simple model that represents a bacteriophage with a round head and a straight line of beads for its tail, used to check binding and track its position in space


## Mathematical Modeling

Our modeling framework builds directly on the conceptual foundations laid out in Schwartz (1976), “Random walk models for phage motion and adsorption,” which established how phages diffuse and attach to bacterial hosts. We extend this idea with a simplified beads-and-sticks representation, where each phage is modeled as a mobile entity with a semi-rigid head and a flexible tail filament. This structural abstraction captures the essential geometry while remaining computationally efficient, allowing us to track phage–bacterium encounter rates with clarity and scalability.

The encounter dynamics combine 3D Brownian motion with optional advection from flow fields, providing flexibility to simulate both quiescent and fluid-driven environments. Adsorption is treated as a contact-based probabilistic process, with binding outcomes governed by stochastic rules that can be naturally extended into continuous-time Markov models. Together, these design choices create a framework that is simple enough for tractable analysis, yet robust enough to accommodate future enhancements in adsorption probability, receptor heterogeneity, and phage–host interaction pathways.
