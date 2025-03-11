/**
# Boussinesq Approximation for Thermal Convection

This module implements the Boussinesq approximation for thermal convection
in the classic Rayleigh-Benard configuration. Temperature is treated as a 
passive tracer that advects with the flow and diffuses, while also generating
buoyancy forces according to the Boussinesq approximation.

The momentum equation with the Boussinesq approximation becomes:
$$
\partial_t \mathbf{u} + \mathbf{u} \cdot \nabla \mathbf{u} = 
-\nabla p + \nu \nabla^2 \mathbf{u} + \mathbf{g} \beta (T - T_0)
$$

where $\beta$ is the thermal expansion coefficient, $T_0$ is the reference temperature,
and $\mathbf{g}$ is the gravitational acceleration vector (pointing downward).

The temperature evolves according to the advection-diffusion equation:
$$
\partial_t T + \mathbf{u} \cdot \nabla T = \kappa \nabla^2 T
$$

where $\kappa$ is the thermal diffusivity.
*/

#include "tracer.h"
#include "diffusion.h"

/**
## Global variables

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

/**
## Physical parameters

These parameters can be modified by the user before initialization.
*/

// Thermal diffusivity (default value, can be modified by the user)
double kappa = 1.0e-2;

// Gravitational acceleration (default value, can be modified by the user)
double AccG = 1.0;

// Thermal expansion coefficient (default value, can be modified by the user)
double beta = 1.0;

// Reference temperature (default value, can be modified by the user)
double T0 = 0.5;

// Rayleigh number = g·β·ΔT·L³/(ν·κ), computed for reference
double Ra = 0.0;

/**
## Boundary conditions for temperature

By default, we set boundary conditions for the classical Rayleigh-Benard configuration:
- Top boundary: Cold (T = 0)
- Bottom boundary: Hot (T = 1)
- Left and right boundaries: Insulated (zero-gradient)
*/

// Cold at the top
T[top] = dirichlet(0.0);

// Hot at the bottom 
T[bottom] = dirichlet(1.0);

// Insulated sides (zero gradient)
T[left] = neumann(0.0);
T[right] = neumann(0.0);

/**
## Initialization

We initialize the temperature field with a linear profile plus a small
random perturbation to kick-start the instability.
*/

event init (i = 0) {
  // Initialize temperature field with linear profile + small perturbation
  foreach() {
    // Linear temperature profile from bottom (T=1) to top (T=0)
    double y_norm = (y + 0.5); // Normalized y-coordinate [0,1]
    T[] = 1.0 - y_norm + 1e-3*noise();  // Small random perturbation
  }
}

/**
## Acceleration term

We add the buoyancy term to the acceleration field of the Navier-Stokes solver.
For Rayleigh-Benard convection, the buoyancy force acts in the y-direction,
with gravity pointing downward (negative y-direction in our coordinate system).
*/

event acceleration (i++) {
  // Check if acceleration field 'a' is allocated or needs to be created
  if (is_constant(a.x)) {
    a = new face vector;
    foreach_face()
      a.x[] = 0.;
  }
  
  // Add buoyancy force to y-component of acceleration
  // We use face averaging between neighboring cells for better accuracy
  foreach_face(y)
    a.y[] += 0.5 * AccG * beta * ((T[] - T0) + (T[0,-1] - T0));
}

/**
## Temperature diffusion

The temperature field diffuses according to thermal diffusivity kappa.
*/
mgstats mgd;
event tracer_diffusion (i++) {
  face vector D[];
  foreach_face()
    D.x[] = kappa; // Constant diffusion coefficient of 1.0
  mgd = diffusion(T, dt, D);
}

/**
## Output

This event computes and outputs the Nusselt number, which is a dimensionless
measure of heat transfer across the layer. For pure conduction, Nu = 1.
For convection, Nu > 1.
*/

event nusselt (t += 0.1) {
  // Compute average heat flux at the bottom
  double flux = 0.0;
  double area = 0.0;
  
  foreach(reduction(+:flux) reduction(+:area)) {
    // Calculate heat flux at the bottom boundary using central difference
    if (y < -0.45) { // Only compute near the bottom boundary
      flux += -kappa * (T[0,1] - T[0,-1]) / (2.0 * Delta);
      area += Delta;
    }
  }
  
  // Nusselt number is the ratio of total heat transfer to conductive heat transfer
  // For a layer with T=1 at bottom and T=0 at top, conductive flux would be 1/L
  double Nu = (area > 0.0) ? (flux / area) * L0 : 0.0;
  fprintf(stderr, "t = %g, Nusselt = %g\n", t, Nu);
} 