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

Subroutines for the incompressible Navier-Stokes equations are also included.

### Governing Equations

The code solves the incompressible Boussinesq equations for a fluid with constant rotation $\mathbf{\Omega} = \Omega \mathbf{\hat{z}}$ and stable stratification characterized by a constant Brunt-Väisälä frequency $N_{BV}$. The equations for the total velocity $\mathbf{u}$, pressure $p$, and buoyancy perturbation $b$ are:

```math
\begin{aligned}
\frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u} \cdot \nabla)\mathbf{u} &= -\nabla p - b\mathbf{\hat{z}} - 2\Omega (\mathbf{\hat{z}} \times \mathbf{u}) + \nu \nabla^2 \mathbf{u} - \nu_p(-\nabla^2)^P \mathbf{u} \\
\frac{\partial b}{\partial t} + (\mathbf{u} \cdot \nabla)b &= N_{BV}^2 (\mathbf{u} \cdot \mathbf{\hat{z}}) + \kappa \nabla^2 b - \kappa_p(-\nabla^2)^P b \\
\nabla \cdot \mathbf{u} &= 0
\end{aligned}
```

where:
- $\nu$ is the kinematic viscosity.
- $\nu_p$ is the hyperviscosity coefficient with power $P$.
- $\kappa$ is the thermal diffusivity.
- $\kappa_p$ is the hyperdiffusivity coefficient.

### Numerical Method

#### Spatial Discretization

The code employs a pseudo-spectral method based on [Matsushima & Marcus (1997)](https://doi.org/10.1006/jcph.1997.5804).

*   **Poloidal-Toroidal Decomposition**: To satisfy the incompressibility condition ($\nabla \cdot \mathbf{u} = 0$) automatically, the velocity field is decomposed into poloidal and toroidal components using scalar potentials $\chi$ (poloidal) and $\psi$ (toroidal):
    ```math
    \mathbf{u} = \nabla \times \nabla \times (\chi \mathbf{\hat{z}}) + \nabla \times (\psi \mathbf{\hat{z}})
    ```
    The evolution equations are solved for $\chi$, $\psi$, and the buoyancy $b$. Note that, the Poloidal-Toroidal projection of a vector field essentially removes any contribution from scalar potentials, and the governing equations become simplified. For example, the projection of $(\mathbf{u} \cdot \nabla)\mathbf{u}$ becomes $\boldsymbol{\omega} \times \mathbf{u}$, where $\boldsymbol{\omega}$ is the vorticity vector, and all terms related to pressure as well as centrifugal acceleration vanish.

*   **Spectral Basis**:
    *   **Azimuthal ($\theta$) and Axial ($z$)**: The fields are represented by Fourier series, which is optimal for periodic boundary conditions in these directions.
    *   **Radial ($r$)**: A mapped coordinate system is used to handle the semi-infinite radial domain $[0, \infty)$. The mapping is of the form $x = (r^2 - L^2)/(r^2 + L^2)$, where $L$ is the map parameter `ELL` from `read.input`. This transforms the radial domain to $x \in [-1, 1]$. The basis functions in this mapped coordinate are Associated Legendre Polynomials. They satisfy the pole condition exactly at $r=0$ and their behavior as $r\rightarrow\infty$ is suitable for expanding smooth functions which decay algebraically or exponentially at far field.
    *   **Logarithmic term**: As discussed in Matsushima & Marcus (1997), logarithmic terms must be included for the poloidal and toroidal fields for completeness, which account for the behavior of the azimuthal velocity at far field and the mean axial components of the velocity and vorticity. It is worth pointing out that the linear terms of the Boussinesq approximations actively contribute to the change of the logarithmic terms of both the poloidal and toroidal fields, which makes its implementation more complicated than the incompressible model.

*   **Pseudo-Spectral Method and Dealiasing**: The nonlinear terms, such as $\mathbf{u} \times \boldsymbol{\omega}$, are calculated in physical space to avoid expensive convolution sums in spectral space. The code transforms the fields from spectral space to physical space (`FFF_SPACE` to `PPP_SPACE`), computes the products, and then transforms back. To prevent aliasing errors, users can either set the number of the spectral modes to 2/3 of the physical collocation points or make use of the hyperviscosity.

#### Time Discretization

The code provides two main time-advancement schemes, controlled by the `BSNSQ%ADAMS` flag.

*   **Scheme 1: Adams-Bashforth Crank-Nicolson (AB2-CN)** (when `BSNSQ%ADAMS = 1`)
    This is a classic semi-implicit scheme.
    *   **Nonlinear and Linear Wave Terms**: The advection term and the linear Coriolis and buoyancy terms are treated explicitly with the second-order Adams-Bashforth method.
    *   **Diffusion Terms**: The viscosity and diffusivity terms are treated implicitly with the second-order, unconditionally stable Crank-Nicolson method to avoid severe time-step restrictions at high resolution. This is implemented in `CALC_BOUSSI_VISC_CN`.

*   **Scheme 2: Exponential Time Differencing (ETD2)** (when `BSNSQ%ADAMS = .FALSE.`)
    This is a more advanced scheme that is particularly effective when stiff linear terms (like fast waves) are present.
    *   **Linear Wave Terms**: The fast linear wave dynamics due to rotation (inertial waves) and stratification (internal gravity waves) are solved exactly using an exponential integrator. The core of this method is in `CALC_BOUSS_DIAG`, which diagonalizes the linear operator `L`, and `CALC_BOUSSI_ETD_OP`, which pre-computes the exponential and related matrix operators.
    *   **Nonlinear Terms**: The nonlinear terms are handled explicitly using a second-order Adams-Bashforth-like formula (ETD2AB).
    *   **Diffusion Terms**: As with the AB2-CN scheme, diffusion is handled separately and implicitly using the Crank-Nicolson method (`CALC_BOUSSI_VISC_CN`).

### Code Workflow and Module Analysis

The simulation is driven by a main program (e.g., `bsnsq_test.f90`) that uses the following modules.

*   **`mod_init`**: This module is responsible for the setup phase. It reads the `read.input` file, initializes all parameters (grid sizes, physical constants, time step), allocates the main field arrays (`PSI`, `CHI`, `B`), and creates the initial conditions.

*   **`mod_march`**: This is the heart of the simulation, controlling the time evolution. It contains the `SOLVER_T` derived type which encapsulates the state and methods for time-stepping, including the AB2-CN and ETD2 schemes. It handles the main time loop, calls diagnostic routines, and manages data output and intermediate data storage.

*   **`mod_boussinesq`**: This module contains the core physics of the Boussinesq equations.
    *   `CALC_BOUSSI_FORCE`: This routine calculates the full right-hand-side (RHS) for the evolution equations, orchestrating the pseudo-spectral calculation of nonlinear terms and the addition of linear terms.
    *   `CALC_BOUSS_DIAG`: This routine is crucial for the ETD scheme. It analytically computes the eigenvalues and eigenvectors of the linear Boussinesq operator (Coriolis and buoyancy forces), allowing the fast wave dynamics to be integrated exactly.

*   **Diagnostics**: Diagnostic routines are present in `mod_march` and `mod_boussinesq`.
    *   `DIAGNOST`: Called periodically to monitor the simulation. It calculates conserved quantities like total energy, helicity, and momentum.
    *   `CALC_BOUSSI_ENERGY`: Calculates the kinetic and potential energy spectra.

### Analysis and Potential Issues

*   **Available Potential Energy Calculation**: In `CALC_BOUSSI_ENERGY`, the potential energy calculated is the available potential energy (APE) defined as $\int B^2 dV / (2N_{BV}^2)$.

*   **Hyperviscosity Implementation**: The code contains two methods for applying hyperviscosity: a standard implicit solver (`HYPERV`) and a more modern spectral filter (`HYPERV3`) that offers better scale selectivity.

*   **Legacy Code**: The `mod_march` module contains time-stepping subroutines for the incompressible model (`STEP_INCOMP_AB_CN`, `RICH_INCOMP_FE_BE`), but they have NOT been ported to the new `SOLVER_T` wrapper.

## Author Information

*   **Name**: Jinge Wang
*   **Affiliation**: University of California, Berkeley
*   **Email**: [jinge@berkeley.edu](mailto:jinge@berkeley.edu)
*   **LinkedIn**: [jinge-wang-cfd](https://www.linkedin.com/in/jinge-wang-cfd/)

*   **Name**: Sangjoon Lee
*   **Affiliation**: Stanford University
