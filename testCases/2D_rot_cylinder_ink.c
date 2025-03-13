#include "grid/multigrid.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "die-injection.h"

// Define inner and outer radii
#define INNER_RADIUS 1.0
#define OUTER_RADIUS 1.5

// Velocidad angular de rotación (rad/s)
#define OMEGA 1.0
#define tsnap 1.0
#define tmax 10.0
const face vector muv[] = {1, 1}; // Ej. nu = 1e-3

char logFile[80];
char dumpFile[80];
char nameOut[80];

int main() {
  // Ajustamos el tamaño del dominio para que sea [-1,1] x [-1,1] (por ejemplo)
 
    L0 = 2.0 * OUTER_RADIUS;
    // Center the domain at origin
    origin (-L0/2, -L0/2);
    // Initialize grid (adjust resolution as needed)
    init_grid (1 << 7);
    
    TOLERANCE = 1e-5;
    CFL = 0.25;
    
    // We use the implicit solver for stability
    DT = 0.01;
    
    // Die injection parameters
    tInjection = 0.05;        // Inject the die after flow is established
    xInjection = 0.00;        // X position (center of cavity)
    yInjection = 1.25;        // Y position (center of cavity)
    
    // Create a folder for simulation snapshots
    char comm[80];
    sprintf (comm, "mkdir -p intermediate");
    system(comm);
    sprintf (dumpFile, "restart");
    sprintf (logFile, "logData.dat");
    
    //mu = muv;
    

  // Iniciamos la simulación
  run();
}

/**
 * En t=0, definimos la geometría del cilindro y su condición de contorno.
 */
event init (t = 0) {
    mu = fm;
    // Inicializamos el campo de velocidad u = 0
    
    
    solid (cs, fs, difference (sq(OUTER_RADIUS) - sq(x) - sq(y),
                             sq(INNER_RADIUS) - sq(x) - sq(y)));
    
    u.n[embed] = dirichlet (x*x + y*y > sq(OUTER_RADIUS-0.1) ? 0. : - y);
    u.t[embed] = dirichlet (x*x + y*y > sq(OUTER_RADIUS-0.1) ? 0. :   x);
    
}

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
    fprintf (stderr, "i t dt\n");
    FILE *fp = fopen (logFile, "w");
    fprintf (fp, "i t dt\n");
    fclose (fp);
  }
  fprintf (stderr, "%d %g %g\n", i, t, dt);
  FILE *fp = fopen (logFile, "a");
  fprintf (fp, "%d %g %g\n", i, t, dt);
  fclose (fp);
}

event timeReversal(t=1){
    u.n[embed] = dirichlet (x*x + y*y > sq(OUTER_RADIUS-0.1) ? 0. : y);
    u.t[embed] = dirichlet (x*x + y*y > sq(OUTER_RADIUS-0.1) ? 0. :  -x);
    
}
