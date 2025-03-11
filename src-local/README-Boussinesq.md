# Boussinesq Approximation for Thermal Convection

## Overview

This module implements the Boussinesq approximation for thermal convection in fluid dynamics simulations using Basilisk C. The Boussinesq approximation assumes that density variations are small and only matter in the buoyancy term of the momentum equation, significantly simplifying the governing equations while still capturing the essential physics of thermal convection.

## Physical Model

The Boussinesq approximation modifies the incompressible Navier-Stokes equations by adding a buoyancy term to the momentum equation:

```
∂u/∂t + (u·∇)u = -∇p/ρ₀ + ν∇²u + g·β·(T-T₀)
```

Where:
- `u` is the velocity field
- `p` is the pressure
- `ρ₀` is the reference density
- `ν` is the kinematic viscosity
- `g` is the gravitational acceleration (implemented as `AccG` in the code)
- `β` is the thermal expansion coefficient
- `T` is the temperature
- `T₀` is the reference temperature

The temperature evolves according to an advection-diffusion equation:

```
∂T/∂t + (u·∇)T = κ∇²T
```

Where `κ` is the thermal diffusivity.

## Key Parameters

- `kappa`: Thermal diffusivity
- `AccG`: Gravitational acceleration (acceleration due to gravity)
- `beta`: Thermal expansion coefficient
- `T0`: Reference temperature

## Implementation Details

### Temperature Diffusion

The temperature diffusion is implemented using Basilisk's implicit diffusion solver. The diffusion equation is solved using:

```c
face vector D[];
foreach_face()
  D.x[] = kappa; // Constant diffusion coefficient
mgd = diffusion(T, dt, D);
```

The `mgd` variable stores multigrid statistics to track convergence of the diffusion solver, which can be useful for debugging or performance analysis.

## Non-dimensional Numbers

Two important non-dimensional parameters characterize thermal convection:

1. **Rayleigh Number (Ra)**: Represents the strength of buoyancy-driven flow
   ```
   Ra = g·β·ΔT·L³/(ν·κ)
   ```
   - When Ra exceeds a critical value (~1708 for a fluid layer heated from below), convection begins
   - Higher Ra values lead to more turbulent convection

2. **Prandtl Number (Pr)**: Ratio of momentum diffusivity to thermal diffusivity
   ```
   Pr = ν/κ
   ```
   - Pr influences flow patterns and heat transfer efficiency
   - Common values: ~0.7 for air, ~7 for water, ~10³ for silicone oils

## Usage

To use this module in your Basilisk code:

1. Include the header file:
   ```c
   #include "convection-Boussinesq.h"
   ```

2. Set parameters before initialization:
   ```c
   event init(t = 0) {
     // Setting physical parameters
     kappa = 0.01;    // Thermal diffusivity 
     AccG = 9.81;     // Gravitational acceleration
     beta = 3e-3;     // Thermal expansion coefficient
     T0 = 0.5;        // Reference temperature
     
     // Other initialization code...
   }
   ```

3. Temperature boundary conditions can be modified:
   ```c
   // Example: Different boundary conditions
   T[top] = dirichlet(0.0);      // Cold top
   T[bottom] = dirichlet(1.0);   // Hot bottom
   T[left] = neumann(0.0);       // Insulated left wall
   T[right] = neumann(0.0);      // Insulated right wall
   ```

## Example: Rayleigh-Bénard Convection

A classic test case is provided in `testCases/2-Rayleigh-Benard.c` that demonstrates Rayleigh-Bénard convection - a fluid layer heated from below and cooled from above.

The example uses:
- Constant viscosity defined by `MU` 
- Thermal diffusivity set to match a Prandtl number of 1.0
- Gravitational acceleration set to `AccG = 1000.0`
- Thermal expansion coefficient of 1.0

To run the example:
```
cd testCases
qcc -I$PWD/../src-local -autolink 2-Rayleigh-Benard.c -o rb -lm
./rb
```

## Output Analysis

The module automatically computes the Nusselt number, which measures the enhancement of heat transfer due to convection:

```
Nu = Total heat flux / Conductive heat flux
```

- Nu = 1: Pure conduction
- Nu > 1: Convection enhances heat transfer

The progress of the simulation can be monitored through:
- Multigrid statistics (`mgd`) for diffusion convergence 
- RMS and maximum velocity values (output in the logfile)
- Nusselt number calculation (output periodically)

## References

1. Chandrasekhar, S. (1961). Hydrodynamic and Hydromagnetic Stability. Oxford University Press.
2. Getling, A. V. (1998). Rayleigh-Bénard Convection: Structures and Dynamics. World Scientific.
3. Ahlers, G., Grossmann, S., & Lohse, D. (2009). Heat transfer and large scale dynamics in turbulent Rayleigh-Bénard convection. Reviews of Modern Physics, 81(2), 503. 