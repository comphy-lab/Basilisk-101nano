/**
 * 2D Transient Heat Conduction Solver
 * 
 * This program solves the transient heat conduction equation in two dimensions:
 * 
 * \frac{\partial T}{\partial t} = \nabla^2 T
 * 
 * Subject to boundary conditions:
 * T = 1 on the top boundary
 * T = 0 on the bottom, left, and right boundaries
 * 
 * Initial condition:
 * T = 0 everywhere
*/

/* Include necessary headers in the correct order for Basilisk */
#include "run.h"
#include "diffusion.h"

// Declare scalar field for temperature
scalar T[];

// Simulation parameters
#define EPS 0.1  // Width of initial temperature peak
#define tmax 1.0
#define tsnap 0.1

mgstats mgd;
char nameOut[80], dumpFile[80], logFile[80];
int main() {
  // Domain setup
  L0 = 10.0;     // Domain length
  X0 = -L0/2;    // Left boundary
  init_grid (1 << 7);      // Number of cells
  
  // We can use a larger timestep with the implicit solver
  // compared to the explicit version which requires dt = dx^2/2
  DT = 0.01;
  
  // Create a folder named intermediate where all the simulation snapshots are stored.
  char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system(comm);
  // Name of the restart file. See writingFiles event.
  sprintf (dumpFile, "restart");
  // Name of the log file. See logWriting event.
  sprintf (logFile, "logData.dat");
  
  // Boundary conditions
  T[top] = dirichlet(1.);
  T[bottom] = dirichlet(0.);
  T[left] = dirichlet(0.);
  T[right] = dirichlet(0.);
  
  // Run simulation
  run();
}

/**
 * Initialize temperature field
 * 
 * Sets up a "Dirac delta" approximated by a thin rectangle
 * centered at x=0 with total integral = 1.
 */
event init (t = 0) {
  foreach()
    T[] = 0.0; 
}

/**
 * Time integration using implicit diffusion solver
 */
event integration (i++) {
  // Get timestep for this iteration
  double dt = dtnext(DT);
  
  // Use the diffusion() function from diffusion.h to solve the equation
  // The heat equation is: ∂T/∂t = ∇²T which corresponds to diffusion with D = 1
  face vector D[];
  foreach_face()
    D.x[] = 1.0; // Constant diffusion coefficient of 1.0
  
  mgd = diffusion(T, dt, D);
}

/**
 * Save snapshots at regular intervals
 */
event writingFiles (t = 0.0; t += tsnap; t < tmax+tsnap) {
  dump (file = dumpFile);
  sprintf (nameOut, "intermediate/snapshot-%5.4f", t);
  dump(file=nameOut);
}

/** 
 * Write logs every timestep about the convergence of the diffusion solver
 * 
 */
event logWriting (i++) {

  if (i == 0) {
    fprintf (stderr, "i t dt mgd.i\n");
    FILE *fp = fopen (logFile, "w");
    fprintf (fp, "i t dt mgd.i\n");
    fclose (fp);
  }
  fprintf (stderr, "%d %g %g %d\n", i, t, dt, mgd.i);
  FILE *fp = fopen (logFile, "a");
  fprintf (fp, "%d %g %g %d\n", i, t, dt, mgd.i);
  fclose (fp);
}