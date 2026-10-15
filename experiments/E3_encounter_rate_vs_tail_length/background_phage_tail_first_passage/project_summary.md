
# Phage-Host Morphology and Adsorption Rates

## Research Objective

To investigate how physical traits of bacteriophages (e.g., tail length, capsid size) and bacterial hosts (e.g., geometry, motility, receptor density) impact the **adsorption rate** of phage to host cells. We aim to quantify these effects using statistical modeling, encounter rate theory, and molecular-scale stochastic analysis.

---

## Data Collection & Processing

### Adsorption Rate Data
- Compiled from peer-reviewed publications across **three major tailed phage families**: Siphoviridae, Myoviridae, and Podoviridae.
- Adsorption constants were **reconstructed** from percent-adsorbed data using the formula:

  \[
  k = -rac{1}{tB} \log(1 - \Delta \hat{P})
  \]

  where:
  - \( t \): incubation time
  - \( B \): bacterial concentration (or derived via MOI)
  - \( \Delta \hat{P} \): percent of phage adsorbed

### Bacterial Receptor Data
- Extracted adsorption vs. receptor density data (e.g., from Schwartz 1976) using MATLAB's `grabit.m`.

---

## Theoretical Framework

### Diffusion of Phage
Phage are approximated as **Brownian spheres**, with diffusion constant:

\[
D_P = rac{k_B T}{6 \pi \eta R_P}
\]

Where:
- \( R_P \): phage radius
- \( \eta \): fluid viscosity
- \( k_B T \): thermal energy

### Bacterial Motility
Bacteria swim at speed \( v_B \), with directional persistence time \( 	au_R \), giving:

\[
D_B = rac{v_B^2 	au_R}{3(1 - lpha)}
\]

Where \( lpha \) is the average cosine of turn angles.

---

## Adsorption Rate Models

### Ballistic Limit:
When bacterial motion dominates:
\[
k_{	ext{ball}} = \pi v_B R_C^2
\]

### Diffusion Limit:
When phage diffusion dominates:
\[
k_{	ext{diff}} = 4 \pi R_C D_P
\]

### With Steric Factor:
Accounts for **phage tail orientation stochasticity**:
\[
k = \psi \cdot k_{	ext{max}} \cdot rac{1 - \gamma}{1 - \gamma \psi}
\]

Where:
- \( \psi \): probability of tail alignment
- \( \gamma \): probability of missing receptors per collision

---

## Advanced Corrections

This framework is conceptually aligned with **first-passage time theory**, where adsorption is modeled as a stochastic boundary-crossing process influenced by both molecular orientation and receptor encounter probability.

- **Sherwood correction** for fluid shear flow around bacteria:

  \[
  	ext{Sh} = \left(rac{R_C v_B}{D_P}ight)^{1/3}, \quad k_{	ext{max}} = k_{	ext{diff}} \cdot 	ext{Sh}
  \]

- **Stochastic escape probability** of phage tail before successful adsorption modeled via:

  \[
  P_{	ext{esc}} = rac{1 - \psi}{1 - \gamma \psi}, \quad k = k_{	ext{max}} (1 - P_{	ext{esc}})
  \]

---

## Assumptions & Uncertainties

### Hard Uncertainties
- Unknown bacterial strain lineage
- Unreported receptor density
- Instrument calibration differences

### Soft Uncertainties
- Missing bacterial motility/geometry
- Incomplete phage morphological data
- Variability in experimental fluid properties

---

## References

- Schwartz, 1976; Moldovan, 2007; Joiner et al., 2019; Berg, 2018
- Data repositories: [Phage-host adsorption data](https://luquelab.github.io/ThunderPhage/approach/phage-host%20adsorption%20data/data.html)

---

## Next Steps

- Validate predicted rates against experimental adsorption datasets.
- Extend model to include environmental variables (viscosity, host density).
- Use Inkscape diagrams and merged project files to generate visualizations.
