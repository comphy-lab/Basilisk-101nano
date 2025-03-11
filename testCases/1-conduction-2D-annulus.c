/**
 * 2D Transient Heat Conduction Solver in an Annulus
 * 
 * This program solves the transient heat conduction equation in an annular domain:
 * 
 * \frac{\partial T}{\partial t} = \nabla^2 T
 * 
 * Subject to boundary conditions:
 * T = 1 on the inner boundary (r = 1)
 * T = 0 on the outer boundary (r = 4)
 * 
 * We use the embedded boundary method to define the annular domain.
*/

/* Include necessary headers in the correct order for Basilisk */
#include "run.h"
#include "embed.h"
#include "diffusion.h"

// Declare scalar field for temperature
scalar T[];

// Define inner and outer radii
#define INNER_RADIUS 1.0
#define OUTER_RADIUS 4.0

// Simulation parameters
#define tmax 100.0
#define tsnap 1.0

mgstats mgd;
char nameOut[80], dumpFile[80], logFile[80];

int main() {
  // Domain setup
  // Make the domain large enough to contain the outer circle
  L0 = 2.0 * OUTER_RADIUS;  
  // Center the domain at origin
  origin (-L0/2, -L0/2);    
  // Initialize grid (adjust resolution as needed)
  init_grid (1 << 7);      
  
  // We use the implicit solver for stability
  DT = 0.01;
  
  // Create a folder for simulation snapshots
  char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system(comm);
  sprintf (dumpFile, "restart");
  sprintf (logFile, "logData.dat");
  
  // Run simulation
  run();
}

/**
 * Initialize temperature field and boundary conditions
 */
event init (t = 0) {
  if (!restore(file = dumpFile)) {
      /**
   * Define the geometry of the embedded boundary as an annulus:
   * - outer circle: sq(OUTER_RADIUS) - sq(x) - sq(y) (positive inside circle of radius 4)
   * - inner circle: sq(INNER_RADIUS) - sq(x) - sq(y) (positive inside circle of radius 1)
   * - difference gives the annular region (positive inside the annulus)
   */
  solid (cs, fs, difference (sq(OUTER_RADIUS) - sq(x) - sq(y),
                           sq(INNER_RADIUS) - sq(x) - sq(y)));
    
  // Set boundary conditions using dirichlet()
  // T = 1 on inner boundary (r = 1)
  // T = 0 on outer boundary (r = 4)
  T[embed] = dirichlet (sq(x) + sq(y) < sq(INNER_RADIUS + 0.1) ? 1.0 : 0.0);
  }
}

/**
 * Time integration using implicit diffusion solver
 */
event integration (i++) {
  // Get timestep for this iteration
  double dt = dtnext(DT);
  
  // Define diffusion coefficient field (constant = 1.0)
  face vector D[];
  foreach_face()
    D.x[] = fs.x[] > 0. ? 1.0 : 0.0; // Diffusion only in fluid cells
  
  // Solve the diffusion equation while respecting embedded boundaries
  mgd = diffusion (T, dt, D);
}

event adapt(i++){
  adapt_wavelet ((scalar *){T, cs},
    (double[]){1e-4, 1e-3},
    9);
  
  solid (cs, fs, difference (sq(OUTER_RADIUS) - sq(x) - sq(y),
                           sq(INNER_RADIUS) - sq(x) - sq(y)));
  T[embed] = dirichlet (sq(x) + sq(y) < sq(INNER_RADIUS + 0.1) ? 1.0 : 0.0);
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