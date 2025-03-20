/**
# Rayleigh-Benard Convection with Boussinesq Approximation

This simulation models Rayleigh-Benard convection in a square domain using
the Boussinesq approximation. A fluid layer is heated from below and cooled 
from above, leading to convection when the Rayleigh number exceeds a critical value.

## Physics
- The Boussinesq approximation treats density as constant except in the buoyancy term
- Temperature boundary conditions: T=1 (hot) at bottom, T=0 (cold) at top, insulated sides
- Velocity boundary conditions: No-slip at all walls

## Parameters
- Rayleigh number: Ra = g·β·ΔT·L³/(ν·κ)
- Prandtl number: Pr = ν/κ
*/

#include "navier-stokes/centered.h"
#include "convection-Boussinesq.h"

// Constants
#define LEVEL   7       // Grid refinement level
#define MAXDT   (1e-3)  // Maximum timestep

// Global variables
int imax = 1e5;                   // Maximum iterations
double tmax = 10.0;               // Maximum simulation time
double tsnap = 0.1;               // Time interval between snapshots
double end = 10.0;                // End time for simulation

// Scalar field for convergence check
scalar un[], vn[];                // Previous velocity components
#define NU (1e-2)
const face vector muv[] = {NU, NU};  // Face-centered viscosity field

/**
## Boundary Conditions
*/
// No-slip boundary conditions at all walls
u.t[top] = dirichlet(0);
u.t[bottom] = dirichlet(0);
u.t[left]   = dirichlet(0);
u.t[right]  = dirichlet(0);

/**
## Initialization
*/
event init (t = 0) {
  // Set constant viscosity
  mu = muv;
  
  // Set Boussinesq parameters
  kappa = 1e-2;  // Thermal diffusivity 
  AccG = 1e3;         // Gravitational acceleration
  beta = 1e0;         // Thermal expansion coefficient
  T0 = 0.5;           // Reference temperature

  fprintf(ferr, "Prandtl number is Pr = %g and Rayleigh number is Ra = %g\n", NU/kappa, AccG*beta/(NU*kappa));
  
  // Initialize velocity field
  foreach() {
    u.x[] = 0;
    u.y[] = 0;
    un[] = 0;
    vn[] = 0;
  }
  
  dump (file = "start");
}

/**
## Main Function
*/
int main() {
  // Initialize grid and parameters
  init_grid(1<<LEVEL);
  L0 = 1.0;
  origin(-0.5, -0.5);
  DT = MAXDT;
  TOLERANCE = 1e-5;
  CFL = 0.8;
  
  // Store current velocity for convergence check
  foreach() {
    un[] = u.x[];
    vn[] = u.y[];
  }

  // Create a folder named intermediate where all the simulation snapshots are stored.
  char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system(comm);
  
  // Run simulation
  run();
}

/**
## Snapshot Generation
Save snapshots at regular intervals for flow visualization
*/
event writingFiles (t=0.; t += tsnap; t < tmax+tsnap) {
  char filename[100];
  sprintf(filename, "intermediate/snapshot-%5.4f", t);  
  dump(file=filename);
}

/**
## Convergence Monitoring
Log information about simulation progress and convergence
*/
event logfile (i++; i <= imax) {
  foreach() {
    un[] = u.x[];
    vn[] = u.y[];
  }
  
  // Calculate RMS velocity and maximum velocity
  double u_rms = 0.0, v_rms = 0.0;
  double u_max = 0.0, v_max = 0.0;
  double total_cells = 0.0;
  
  foreach(reduction(+:u_rms) reduction(+:v_rms) 
          reduction(+:total_cells) reduction(max:u_max) reduction(max:v_max)) {
    u_rms += sq(u.x[]);
    v_rms += sq(u.y[]);
    if (fabs(u.x[]) > u_max) u_max = fabs(u.x[]);
    if (fabs(u.y[]) > v_max) v_max = fabs(u.y[]);
    total_cells += 1.0;
  }
  
  // Calculate root mean square
  u_rms = sqrt(u_rms/total_cells);
  v_rms = sqrt(v_rms/total_cells);
  
  fprintf(ferr, "i = %d: dt = %g, t = %g, |u|_rms = %g, |v|_rms = %g, |u|_max = %g, |v|_max = %g\n", 
          i, dt, t, u_rms, v_rms, u_max, v_max);
}

/**
## Output & Visualization
Generate final output for post-processing and visualization
*/
event end (t = end) {  
  // Output fields in a format suitable for visualization
  dump(file="results");
}