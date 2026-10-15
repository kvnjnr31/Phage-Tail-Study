# Bacteriophage Encounter Dynamics (MATLAB to Python Conversion)

## Objective

This simulation models bacteriophage-particle encounters with diffusing bacteria using stochastic particle-based methods. Originally implemented in MATLAB using `encounter_rates_hydro.m` and `encounter_fun_hydro.m`, this Python version provides modular, scalable code for simulating residence time and encounter rates.

---

## Core Model Parameters

- `P = 1e5`: Number of phage particles
- `D0 = 2.54`: Diffusion constant
- `omega = 1`: Bacterial tumble rate
- `dt = 0.01`: Time step
- `Rcrit`: Critical interaction radius
- `Rp`: Particle radius
- `nu = 1`: Diffusion exponent
- `D_para_y`, `D_para_z`: Toggles for dimensionality

---

## Core Structures

### `MatP` - Phage Particle Matrix
| Column | Description           |
|--------|------------------------|
| 0      | x position             |
| 1      | y position             |
| 2      | z position             |
| 4      | is_active flag         |
| 5      | residence time tracker |

---

## Main Functions

### `encounter_fun_hydro`
This function checks particle distance to a central point and:
- Updates residence time if within `Rcrit + Rp`
- Counts conditional encounters (if distance ≤ `Rcrit`)
- Returns:
  - `Nenc`: Active encounters
  - `Nact`: Active phages
  - `N_dep`: Deposited phage %
  - Updated `MatP`, `res_times`, `dep_vels`

---

## Python Port Plan

1. Set parameters and initialize `MatP`
2. Loop over time steps:
   - Update particle positions (Brownian motion)
   - Call `encounter_fun_hydro`
3. Collect and visualize results

---

## Next Steps

- Implement `update_particle_positions`
- Extend logging of encounter/residence stats
- Visualize spatial distributions and time evolution
---

## Methods

### 🔧 Design Upgrade: Modular Class-Based Architecture

To support future extensions—including **parallel 3D simulation**, **GPU acceleration**, and **clean integration of phage-bacteria interaction logic**—we recommend transitioning to a class-based design in Python.

#### ✅ Why Class-Based?

- **Encapsulates simulation state** (`MatP`, `MatB`, `Time`, etc.)
- Enables **parallelism and multiprocessing**
- Organizes logic into clean, testable units
- Supports reusable features like logging, plotting, and data export
- Simplifies future conversions to GPU (Numba/CUDA) or distributed environments

#### 🧱 Recommended Structure

```python
class PhageSimulator:
    def __init__(self, B, P, params):
        self.B = B
        self.P = P
        self.params = params
        self._initialize_system()

    def _initialize_system(self):
        # Setup MatB, MatP
        pass

    def update_positions(self):
        pass

    def generate_noise(self):
        pass

    def resolve_overlap(self):
        pass

    def run(self, nsteps):
        for _ in range(nsteps):
            self.update_positions()
            self.resolve_overlap()
            self.generate_noise()
```

You can define this in a separate file (e.g., `phage_simulator.py`) and import it into your main script for clean separation between **simulation engine** and **execution control**.
