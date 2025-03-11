/**
# Two-Phase Flow with Heat Transfer

This file extends Basilisk's two-phase flow capabilities to include heat transfer,
allowing simulation of flows where thermal effects are important, such as boiling,
condensation, Marangoni flows, and thermocapillary migration.

## Physical Model

The system tracks two immiscible fluids using a Volume-Of-Fluid (VOF) approach:
- `f`: Volume fraction field (f=1 in phase 1, f=0 in phase 2)

Each phase has physical properties:
- Densities: *rho1*, *rho2*
- Dynamic viscosities: *mu1*, *mu2*
- Thermal diffusivities: *D1*, *D2*

## Heat Transfer Equation

The temperature evolves according to the advection-diffusion equation:
$$
\partial_t T + \mathbf{u} \cdot \nabla T = \nabla \cdot (D \nabla T)
$$

where:
- $T$ is the temperature field
- $\mathbf{u}$ is the velocity field (from Navier-Stokes)
- $D$ is the thermal diffusivity, which varies between the two phases

The implementation handles:
1. Advection of temperature with the fluid flow
2. Diffusion of heat with phase-dependent thermal diffusivity
3. Sharp temperature gradients across fluid interfaces

## Implementation Details

The framework combines several key components:
- VOF interface tracking (`vof.h`)
- Scalar field advection (`tracer.h`)
- Diffusion solver (`diffusion.h`)

It provides:
1. Proper averaging of thermal properties at phase interfaces
2. Optional filtering for numerical stability
3. Automatic coupling with Basilisk's adaptive mesh refinement
4. Properly constructed face and cell-centered fields for the solvers

## Usage

This header should be included after the Navier-Stokes solver but before
any problem-specific code:

```c
#include "navier-stokes/centered.h"
#include "two-phase-thermal.h"
```

Then set the physical properties for both phases:
```c
// Phase 1 properties
rho1 = 1.0;  mu1 = 0.01;  D1 = 0.005;
// Phase 2 properties
rho2 = 0.001;  mu2 = 0.0001;  D2 = 0.001;
```

The temperature field `T` is automatically added to the list of tracers
that will be advected with the flow and will have diffusion applied.
*/

#include "tracer.h"
#include "diffusion.h"
#include "vof.h"

scalar f[], * interfaces = {f};

/**
## Global Variables and Fields

We define a scalar field for temperature and add it to the list of tracers
that will be advected with the flow.
*/

scalar T[];          // Temperature field
// Add T to list of tracers if tracers are not already defined
#ifndef tracers
scalar * tracers = {T};
#else
event defaults (i = 0) {
  tracers = list_append (tracers, T);
}
#endif

// Physical properties for both phases
double rho1 = 1., mu1 = 0., rho2 = 1., mu2 = 0.;  // Density and viscosity
double D1 = 0., D2 = 0.;                          // Thermal diffusivity

/**
## Auxiliary Fields

These fields are necessary for the solvers:
- `alphav`: Face-centered specific volume (1/ρ) for momentum equation
- `rhov`: Cell-centered density field
- `D`: Face-centered thermal diffusivity for heat equation
*/

face vector alphav[];  // Face-centered specific volume
scalar rhov[];         // Cell-centered density
face vector D[];       // Face-centered thermal diffusivity

/**
## Initialization

This event sets up the necessary fields during simulation initialization.
*/
event defaults (i = 0) {
  alpha = alphav;     // Connect to Navier-Stokes solver
  rho = rhov;         // Connect to Navier-Stokes solver

  /**
  If viscosity is non-zero, allocate the face-centered viscosity field
  needed by the Navier-Stokes solver.
  */
  mu = new face vector;
}

/**
## Property Calculations

These macros define how physical properties are calculated in mixed cells
(cells containing interfaces). The default implementation uses arithmetic
averaging, but users can define their own averaging methods if needed.
*/

#ifndef rho
# define rho(f) (clamp(f,0.,1.)*(rho1 - rho2) + rho2)
#endif
#ifndef mu
// Arithmetic mean for viscosity
# define mu(f)  (clamp(f,0.,1.)*(mu1 - mu2) + mu2)
#endif
#ifndef D
// Arithmetic mean for thermal diffusivity
# define D(f)  (clamp(f,0.,1.)*(D1 - D2) + D2)
#endif

/**
## Optional Filtering

For improved numerical stability, we can apply filtering (smearing) to the
volume fraction field before calculating properties. This reduces issues
from sharp property jumps at interfaces.
*/

#ifdef FILTERED
scalar sf[];  // Smoothed volume fraction for property calculations
#else
# define sf f  // Without filtering, use original volume fraction
#endif

/**
## Tracer Advection

This event handles the filtering of the volume fraction field for property
calculations. It runs once per timestep before the properties are calculated.
*/
event tracer_advection (i++) {

  /**
  When filtering is enabled, calculate a smoothed volume fraction field
  using a weighted average of neighboring cells. This creates a more
  gradual transition of properties across interfaces.
  
  The stencil applies different weights to:
  - Center cell: highest weight
  - Face neighbors: medium weight 
  - Corner neighbors: lowest weight
  */

#ifndef sf
#if dimension <= 2
  foreach()
    sf[] = (4.*f[] +                                          // Center (weight 4)
	    2.*(f[0,1] + f[0,-1] + f[1,0] + f[-1,0]) +        // Faces (weight 2)
	    f[-1,-1] + f[1,-1] + f[1,1] + f[-1,1])/16.;       // Corners (weight 1)
#else // dimension == 3
  foreach()
    sf[] = (8.*f[] +                                          // Center (weight 8)
	    4.*(f[-1] + f[1] + f[0,1] + f[0,-1] + 
             f[0,0,1] + f[0,0,-1]) +                          // Faces (weight 4)
	    2.*(f[-1,1] + f[-1,0,1] + f[-1,0,-1] + f[-1,-1] +
		f[0,1,1] + f[0,1,-1] + f[0,-1,1] + f[0,-1,-1] +
		f[1,1] + f[1,0,1] + f[1,-1] + f[1,0,-1]) +      // Edges (weight 2)
	    f[1,-1,1] + f[-1,1,1] + f[-1,1,-1] + f[1,1,1] +
	    f[1,1,-1] + f[-1,-1,-1] + f[1,-1,-1] + 
            f[-1,-1,1])/64.;                                  // Corners (weight 1)
#endif
#endif

#if TREE
  // For adaptive mesh refinement, set appropriate interpolation method
  sf.prolongation = refine_bilinear;
  sf.dirty = true; // Mark that boundary conditions need updating
#endif
}

/**
## Properties Update

This event calculates all physical properties at each timestep based on
the current interface position. It runs once per timestep and updates
all the fields needed by the solvers.
*/
event properties (i++) {
  
  /**
  Calculate face-centered properties needed by the solvers:
  1. Get average volume fraction at each face
  2. Compute density and convert to specific volume
  3. Calculate face-centered viscosity
  4. Calculate face-centered thermal diffusivity
  */
  foreach_face() {
    double ff = (sf[] + sf[-1])/2.;          // Average f at the face
    alphav.x[] = fm.x[]/rho(ff);             // Specific volume scaled by face fraction
    face vector muv = mu;
    muv.x[] = fm.x[]*mu(ff);                 // Viscosity scaled by face fraction
    D.x[] = fm.x[]*D(ff);                    // Thermal diffusivity scaled by face fraction
  }

  /**
  Calculate cell-centered density for each cell.
  This is scaled by the cell fraction (cm) which is important for
  embedded boundary methods.
  */
  foreach(){
    rhov[] = cm[]*rho(sf[]);                 // Density scaled by cell fraction
  }

#if TREE
  // Restore proper prolongation method for VOF after filtering
  sf.prolongation = fraction_refine;
  sf.dirty = true; // Mark that boundary conditions need updating
#endif
}

/**
## Temperature Diffusion

This event solves the diffusion term of the temperature equation.
It uses the diffusion solver with the face-centered thermal diffusivity
calculated in the properties event.
*/
event tracer_diffusion (i++) {
  diffusion(T, dt, D);  // Solve ∂T/∂t = ∇·(D∇T) for one timestep
}
