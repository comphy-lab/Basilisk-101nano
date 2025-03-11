/* Title: Drop Impact
*/

#include "axi.h"
#include "navier-stokes/centered.h"
// #define FILTERED
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "reduced.h"
#include "distance.h"
#include "adapt_wavelet_limited.h"

#define tsnap (0.0025)
// Error tolerancs
#define fErr (1e-3)                                 // error tolerance in f1 VOF
#define KErr (1e-4)                                 // error tolerance in f2 VOF
#define VelErr (1e-4)                            // error tolerances in velocity

#define Hint 1.5
// Ratios:
#define RHO21 (1e-3)
#define MU21 (1e-2)
// Numbers:
#define Ga (1e4) // Gallelli number
#define Bo (1e1) // Bond number
#define Ldomain 8

// boundary conditions: right part of the domain is outflow
u.n[right] = neumann(0.);
p[right] = dirichlet(0.); //p is reduced pressure

int MAXlevel;
double tmax, Fr;
int  main(int argc, char const *argv[]) {
  dtmax = 1e-3;
  L0 = Ldomain;
  origin (-L0/2., 0.);
  init_grid (1 << (6));
  // Values taken from the terminal
  MAXlevel = atoi(argv[1]);
  tmax = atof(argv[2]);
  Fr = atof(argv[3]);
  fprintf(ferr, "Level %d, tmax %f, Fr %f, Ga %g, Bo %g\n", MAXlevel, tmax, Fr, Ga, Bo);
  rho1 = 1., rho2 = RHO21;
  mu1 = 1./sqrt(Ga), mu2 = MU21/sqrt(Ga);
  f.sigma = 1./Bo;
  G.x = -1.;
  run();
}

int refRegion(double x, double y, double z){
  return (y < 2 ? MAXlevel :  y < 4 ? MAXlevel-1 : y < 6 ? MAXlevel-2 : y < 8 ? MAXlevel-3: y < 10 ? MAXlevel-5: MAXlevel-5);
}

event init (t = 0) {
  if (!restore (file = "dump")){
    refine (sq(x-Hint) + sq(y) < 1.1 && level < MAXlevel);
    fraction(f, 1 - sq(x-Hint) - sq (y));
    scalar m[];
    fraction(m, -x);// and the pool
    foreach(){
      u.x[] = -f[]*sqrt(Fr);
      f[] += m[];
    }
    boundary((scalar *){f, u.x, u.y});
  }
}


event adapt(i++){
  scalar KAPPA[], ux[], uy[];
  curvature(f, KAPPA);
  foreach(){
    ux[] = rho(f[])*u.x[];
    uy[] = rho(f[])*u.y[];
  }
  boundary((scalar*){f, u.x, u.y, ux, uy, KAPPA});
  adapt_wavelet_limited ((scalar *){f, ux, uy, KAPPA},
     (double[]){fErr, VelErr, VelErr, KErr},
      refRegion, 0);
}

// Outputs
event writingFiles (t = 0; t += tsnap; t <= tmax) {
  dump (file = "dump");
  char nameOut[60];
  sprintf (nameOut, "intermediate/snapshot-%5.4f", t);
  dump(file=nameOut);
}

event logWriting (i++) {
  double ke = 0.;
  foreach (reduction(+:ke)){
    ke += (2*pi*y)*(0.5*(f[])*rho1*(sq(u.x[]) + sq(u.y[])))*sq(Delta);
  }
  static FILE * fp;
  if (i == 0) {
    fprintf (ferr, "i dt t ke\n");
    fp = fopen ("log", "w");
    fprintf (fp, "i dt t ke\n");
    fprintf (fp, "%d %g %g %g\n", i, dt, t, ke);
    fclose(fp);
  } else {
    fp = fopen ("log", "a");
    fprintf (fp, "%d %g %g %g\n", i, dt, t, ke);
    fclose(fp);
  }
  fprintf (ferr, "%d %g %g %g\n", i, dt, t, ke);
}
