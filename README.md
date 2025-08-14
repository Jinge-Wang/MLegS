# MLegS
Updated Mapped Legendre Polynomial spectral CFD code with python interface

## Code Deployment

After downloading the code package, grant execution permission to the deployment script:

```bash
chmod +x deploy.sh
```

Then run the script to build all required libraries and set up the Conda environment for Jupyter Notebook data analysis:

```bash
./deploy.sh
```

The deployment system automatically detects whether you are on macOS or Linux (including supercomputers) and builds the programs accordingly. **Windows is no longer supported.**

## Compiling Programs

To compile the Fortran programs, run `make program_name` from the root directory. For example:

```bash
make bsnsq_test
```

You can speed up compilation on multicore systems by using the `-j` flag:

```bash
make -j bsnsq_test
```

This will build the specified Fortran executable in the `f90` directory.

## Code Architecture and Numerical Methods

### Overview

The code package implements a pseudo-spectral solver for the three-dimensional Boussinesq equations in a cylindrical geometry. It is designed for studying rotating and stratified flows, such as those found in geophysical and astrophysical contexts. The velocity field is decomposed into poloidal and toroidal components, which is a standard technique for incompressible flows. The spatial discretization uses a pseudo-spectral method with Fourier series in the azimuthal and axial directions, and mapped associated Legendre polynomials for the radial direction. The code includes multiple time-stepping schemes, including a second-order Adams-Bashforth/Crank-Nicolson (AB2-CN) scheme and a more advanced second-order Exponential Time Differencing (ETD) scheme.

### Governing Equations

The code solves the incompressible Boussinesq equations for a fluid with constant rotation **Ω** = Ω **ẑ** and stable stratification characterized by a constant Brunt-Väisälä frequency N_BV. The equations for the velocity perturbation **u**, pressure perturbation p, and buoyancy perturbation b are:

```
∂**u**/∂t + (**u** ⋅ ∇)**u** = -∇p + b**ẑ** + 2Ω (**u** × **ẑ**) + ν∇²**u** - ν_p(-∇²)^P **u**
∂b/∂t + (**u** ⋅ ∇)b = -N_BV² (**u** ⋅ **ẑ**) + κ∇²b - κ_p(-∇²)^P b
∇ ⋅ **u** = 0
```

where:
- ν is the kinematic viscosity.
- ν_p is the hyperviscosity coefficient with power P.
- κ is the thermal diffusivity.
- κ_p is the hyperdiffusivity coefficient.

### Numerical Method

#### Spatial Discretization

The code employs a sophisticated pseudo-spectral method.

*   **Poloidal-Toroidal Decomposition**: To satisfy the incompressibility condition (∇ ⋅ **u** = 0) automatically, the velocity field is decomposed into poloidal and toroidal components using scalar potentials χ (poloidal) and ψ (toroidal):
    ```
    **u** = ∇ × ∇ × (χ **ẑ**) + ∇ × (ψ **ẑ**)
    ```
    The evolution equations are solved for χ, ψ, and the buoyancy b.

*   **Spectral Basis**:
    *   **Azimuthal (θ) and Axial (z)**: The fields are represented by Fourier series, which is optimal for periodic boundary conditions in these directions.
    *   **Radial (r)**: A mapped coordinate system is used to handle the semi-infinite radial domain [0, ∞). The mapping is of the form `x = (r - L)/(r + L)` or similar, where `L` is the map parameter `ELL` from `read.input`. This transforms the radial domain to `x ∈ [-1, 1]`. The basis functions in this mapped coordinate are Associated Legendre Polynomials. This combination is powerful for capturing boundary layers and resolving features in a large domain.

*   **Pseudo-Spectral Method and Dealiasing**: The nonlinear terms, such as (**u** ⋅ ∇)**u**, are calculated in physical space to avoid expensive convolution sums in spectral space. The code transforms the fields from spectral space to physical space (`FFF_SPACE` to `PPP_SPACE`), computes the products, and then transforms back. To prevent aliasing errors from the quadratic nonlinearities, the code uses a 3/2 padding rule, indicated by `CHOPSET(3)`. This means the physical space grid is 50% larger than what is required to represent the truncated spectral modes.

#### Time Discretization

The code provides two main time-advancement schemes, controlled by the `BSNSQ%ADAMS` flag.

*   **Scheme 1: Adams-Bashforth Crank-Nicolson (AB2-CN)** (when `BSNSQ%ADAMS = .TRUE.`)
    This is a classic semi-implicit scheme.
    *   **Nonlinear and Linear Wave Terms**: The advection term and the linear Coriolis and buoyancy terms are treated explicitly with the second-order Adams-Bashforth method.
    *   **Diffusion Terms**: The viscosity and diffusivity terms are treated implicitly with the second-order, unconditionally stable Crank-Nicolson method to avoid severe time-step restrictions at high resolution. This is implemented in `VISC2`.

*   **Scheme 2: Exponential Time Differencing (ETD2)** (when `BSNSQ%ADAMS = .FALSE.`)
    This is a more advanced scheme that is particularly effective when stiff linear terms (like fast waves) are present.
    *   **Linear Wave Terms**: The fast linear wave dynamics due to rotation (inertial waves) and stratification (internal gravity waves) are solved exactly using an exponential integrator. The core of this method is in `CALC_BOUSS_DIAG`, which diagonalizes the linear operator `L`, and `ETD_INIT`, which pre-computes the exponential and related matrices.
    *   **Nonlinear Terms**: The nonlinear terms are handled explicitly using a second-order Adams-Bashforth-like formula (ETD2AB).
    *   **Diffusion Terms**: As with the AB2-CN scheme, diffusion is handled separately and implicitly using the Crank-Nicolson method (`VISC2`).

### Code Workflow and Module Analysis

The simulation is driven by a main program (e.g., `bsnsq_test.f90`) that uses the following modules.

*   **`mod_init`**: This module is responsible for the setup phase. It reads the `read.input` file, initializes all parameters (grid sizes, physical constants, time step), allocates the main field arrays (`PSI`, `CHI`, `B`), and creates the initial conditions.

*   **`mod_march`**: This is the heart of the simulation, controlling the time evolution. It contains the `SOLVER_T` derived type which encapsulates the state and methods for time-stepping, including the AB2-CN and ETD2 schemes. It handles the main time loop, calls diagnostic routines, and manages data output.

*   **`mod_boussinesq`**: This module contains the core physics of the Boussinesq equations.
    *   `BOUSSINESQ_FULL`: This routine calculates the full right-hand-side (RHS) for the evolution equations, orchestrating the pseudo-spectral calculation of nonlinear terms and the addition of linear terms.
    *   `CALC_BOUSS_DIAG`: This routine is crucial for the ETD scheme. It analytically computes the eigenvalues and eigenvectors of the linear Boussinesq operator (Coriolis and buoyancy forces), allowing the fast wave dynamics to be integrated exactly.

*   **Diagnostics**: Diagnostic routines are present in `mod_march` and `mod_boussinesq`.
    *   `DIAGNOST`: Called periodically to monitor the simulation. It calculates conserved quantities like total energy, helicity, and momentum.
    *   `CALC_ENERGY`: Calculates the kinetic and potential energy spectra.

### Analysis and Potential Issues

The code is a high-quality, research-grade simulation tool. The numerical methods are sophisticated and well-suited for the target problem.

*   **Potential Energy Calculation**: In `CALC_ENERGY`, the potential energy (PE) is calculated based on the vertical center of mass of the buoyancy perturbation (`INT[Z*B']`), not the standard available potential energy proportional to `∫ B² dV`. This diagnostic is valid, but it's not the energy component that is directly converted from kinetic energy.

*   **Hyperviscosity Implementation**: The code contains two methods for applying hyperviscosity: a standard implicit solver (`HYPERV`) and a more modern spectral filter (`HYPERV3`) that offers better scale selectivity.

*   **Legacy Code**: The `mod_march` module contains some older, likely unused time-stepping subroutines (`ADAMSB`, `RICH`, `EULER`) alongside the newer `SOLVER_T` implementation. Refactoring to remove this legacy code would improve maintainability.
