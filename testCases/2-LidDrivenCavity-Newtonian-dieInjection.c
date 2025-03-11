/**
# Lid-Driven Cavity Flow of a Newtonian Fluid with Die Injection

This simulation models a lid-driven cavity flow for a Newtonian fluid with
constant viscosity and includes die injection for flow visualization. 
This extends the classic benchmark case with a passive tracer to visualize flow patterns.

## Parameters
- Reynolds number: Re = ρUL/μ = 1/μ (with ρ=1, U=1, L=1)
- We use μ = 1.0 by default (Re = 1)
- Die injection at t=0.05 in the upper center of the cavity
*/

#include "navier-stokes/centered.h"
#include "die-injection.h"

// Constants
#define LEVEL   7       // Grid refinement level
#define MAXDT   (1e-4)  // Maximum timestep

// Global variables
int imax = 1e5;                   // Maximum iterations
double tmax = 2.0;                // Maximum simulation time
double tsnap = 0.05;              // Time interval between snapshots
double end = 2.0;                 // End time for simulation

// Scalar field for convergence check
scalar un[];  // Previous x-velocity
const face vector muv[] = {1.0, 1.0};      // Face-centered viscosity field

/**
## Boundary Conditions
*/
// Top moving wall (lid)
u.t[top] = dirichlet(1);
// Other no-slip boundaries
u.t[bottom] = dirichlet(0);
u.t[left]   = dirichlet(0);
u.t[right]  = dirichlet(0);
uf.n[left]   = 0;
uf.n[right]  = 0;
uf.n[top]    = 0;
uf.n[bottom] = 0;

/**
## Initialization
*/
event init (t = 0) {
  // Set constant viscosity for Newtonian fluid
  mu = muv;
  
  // Initialize velocity field
  foreach() {
    u.x[] = 0;
    u.y[] = 0;
    un[] = 0;
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
  CFL = 0.25;
  
  // Store current velocity for convergence check
  foreach() {
    un[] = u.x[];
  }

  // Die injection parameters
  tInjection = 0.05;        // Inject the die after flow is established
  xInjection = 0.00;        // X position (center of cavity)
  yInjection = 0.40;        // Y position (center of cavity)

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
  }
  fprintf(ferr, "i = %d: dt = %g, t = %g, err = %g\n", i, dt, t, change(u.x, un));
}

/**
## Output & Visualization
Generate final output for post-processing and visualization
*/
event end (t = end) {  
  // Output fields in a format suitable for visualization
  dump(file="results");
}