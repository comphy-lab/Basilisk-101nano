/**
 # 2D Transient Heat Conduction Solver
 
 This program solves the transient heat conduction equation in two dimensions:
 
 $$
 \frac{\partial T}{\partial t} = \nabla^2 T
 $$
 
 ## Subject to boundary conditions:

  - T = 1 on the top boundary
  - T = 0 on the bottom, left, and right boundaries
 
 ## Initial condition:
  - T = 0 everywhere

 ## To use make file, do

 ```bash 
 CFLAGS=-DDISPLAY=-1 make 1-conduction-2D.tst
 ```

*/

// Include necessary headers in the correct order for Basilisk
#include "run.h"
#include "diffusion.h"

// Declare scalar field for temperature
scalar T[];
// Boundary conditions
T[top] = dirichlet(1.);
T[bottom] = dirichlet(0.);
T[left] = dirichlet(0.);
T[right] = dirichlet(0.);

// Simulation parameters
#define tmax 10.0
#define tsnap 1.0

mgstats mgd;
char nameOut[80], dumpFile[80], logFile[80];
int main() {
  // Domain setup
  L0 = 8.0;     // Domain length
  X0 = -L0/2;    // Left boundary
  Y0 = -L0/2;    // Bottom boundary
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

  
  // Run simulation
  run();
}

/**
 ## Initialize temperature field
 
 Sets up a "Dirac delta" approximated by a thin rectangle
 centered at x=0 with total integral = 1.
 */
event init (t = 0) {
  foreach()
    T[] = 0.0; 
}

/**
 ## Time integration using implicit diffusion solver
 */
event integration (i++) {
  // Get timestep for this iteration
  double dt = dtnext(DT);
  mgd = diffusion(T, dt);
}

event adapt(i++){
  adapt_wavelet ((scalar *){T},
    (double[]){1e-4},
    10);
}

/**
 ## Save snapshots at regular intervals
 */
event writingFiles (t = 0.0; t += tsnap; t < tmax+tsnap) {
  dump (file = dumpFile);
  sprintf (nameOut, "intermediate/snapshot-%5.4f", t);
  dump(file=nameOut);
}

/** 
 ## Log writing
 
 Write logs every timestep about the convergence of the diffusion solver
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