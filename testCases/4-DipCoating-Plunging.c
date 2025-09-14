/**

@file 4-DipCoating-Withdrawal.c
@brief Simulation of dip coating withdrawal - Landau-Levich problem

This code simulates the classical Landau-Levich dip coating process where
a vertical plate is withdrawn from a liquid bath at constant velocity.
The simulation captures the formation of the liquid film on the plate
and the meniscus dynamics.

Problem scales: 
 - Length: L = l_c = \sqrt{\sigma/\rho g}
 - Velocity: \sigma/\mu
 - Time: l_c/(\sigma/\mu)
 - Pressure: \rho g l_c

Physical parameters:
 - Capillary number (Ca): Viscous vs surface tension forces
 - Density ratio: Liquid vs air density
 - Viscosity ratio: Liquid vs air viscosity
 - Outer length scale: L0/l_c
 
Key physics:
 - Contact line dynamics at the three-phase contact line
 - Meniscus formation and shape
 - Steady-state film formation


@author Vatsal Sanjay (vatsal.sanjay@comphy-lab.org)
https://comphy-lab.org
@date 2025-09-14
@version 1.0

*/

// Phase labeling: 1 is liquid, 2 is air

// ======= Include necessary Basilisk modules =======
#include "navier-stokes/centered.h"    // NS solver with centered discretization
#define FILTERED 1                     // Use filtered VOF advection
#include "two-phase.h"                 // Two-phase interface tracking
#include "navier-stokes/conserving.h"  // Conservative momentum advection
#include "tension.h"                   // Surface tension model
#include "reduced.h"                   // Reduced gravity approach

// ======= Numerical parameters for adaptivity =======
#define fErr (1e-3)             // Error tolerance in Volume of Fluid (interface position)
#define KErr (1e-6)             // Error tolerance in curvature (KAPPA)
#define VelErr (1e-2)           // Error tolerances in velocity fields
#define DissErr (1e-3)          // Error tolerances in dissipation rate

// ======= Geometry parameters =======
#define BathLevel (4.0)         // Initial liquid level in the bath

// ======= Global parameters =======
int MAXlevel;                        // Maximum refinement level
double tmax, Ca, Rho21, Mu21, Ldomain;
#define MINlevel 3               // Minimum refinement level
#define tsnap (0.1)             // Time interval for snapshot outputs

// ======= Boundary conditions =======
// Left boundary: moving plate with withdrawal velocity
u.t[left] = dirichlet(-Ca);      // Tangential velocity = withdrawal velocity
u.n[left] = dirichlet(0.0);         // Normal velocity = 0 (no-slip)
f[left] = neumann(0.0);           // No-flux for volume fraction (90 degree contact angle)

// Right boundary: symmetric boundary
// Bottom boundary: symmetric boundary
// Top boundary: outflow condition
u.n[top] = neumann(0.);             // Zero gradient for normal velocity
u.t[top] = neumann(0.);            // Zero gradient for tangential velocity
p[top] = dirichlet(0.0);            // Reference pressure = 0 (outflow)


int main(int argc, char const *argv[]) {

  // Default parameter values
  MAXlevel = 8;                  // Maximum grid refinement level
  tmax = 100.0;                   // Maximum simulation time (dimensionless)
  Ca = 0.1;                      // Capillary number (Ca = \mu U/\sigma)
  Mu21 = 1e-2;
  Rho21 = 1e-2;
  Ldomain = 1e1;                 // Domain size (dimensionless)

  // Log the simulation parameters
  fprintf(ferr, "Level %d tmax %g. Ca %g, Mu21 %g, Rho21 %g, Ldomain %g\n",
          MAXlevel, tmax, Ca, Mu21, Rho21, Ldomain);

  // ======= Set up the computational domain =======
  L0 = Ldomain;                  // Domain size
  X0 = 0.; Y0 = 0.;              // Domain origin
  init_grid(1 << (4));           // Start with a 16�16 base grid (2^4)

  // Create directory for intermediate results
  char comm[80];
  sprintf(comm, "mkdir -p intermediate");
  system(comm);

  // ======= Set physical properties =======
  // Note: All quantities are in dimensionless form
  rho1 = 1.0;                    // Density of liquid (normalized)
  rho2 = Rho21;                  // Density of air
  mu1 = 1e0;                   // Viscosity of liquid (from Ca and Re)
  mu2 = Mu21;                // Viscosity of air
  f.sigma = 1.0;                 // Surface tension coefficient (normalized)
  G.y = -1e0;                     // Gravity in y-direction (from Bond number)

  run();
}

// ======= Initialize the simulation =======
event init(t = 0){
  if(!restore(file = "dump")){
    // Refine mesh near the initial interface (meniscus region)
    refine((y < BathLevel + 0.2) && (x < 1.0) && (level < MAXlevel));

    // Initialize the liquid bath
    // f=1 in liquid, f=0 in air
    fraction(f, BathLevel-y);

    // Set initial velocity field
    foreach () {
      u.x[] = 0.0;     // Liquid moves with plate initially
      u.y[] = 0.0;               // No initial y-velocity
    }

  }
}


// ======= Adaptive mesh refinement =======
scalar KAPPA[], D2c[];    // Curvature field and dissipation field
event adapt(i++){
  // Calculate curvature for interface refinement
  curvature(f, KAPPA);

  // Calculate local dissipation rate as refinement criteria
  foreach(){
    // Calculate velocity gradient components
    double D11 = (u.y[0,1] - u.y[0,-1])/(2*Delta);    // u_y/y
    double D22 = (u.x[1,0] - u.x[-1,0])/(2*Delta);    // u_x/x
    double D12 = 0.5*((u.y[1,0] - u.y[-1,0] + u.x[0,1] - u.x[0,-1])/(2*Delta)); // Shear rate

    // Calculate dissipation rate (sum of squares of strain rates)
    double D2 = 2.0*(sq(D11) + sq(D22) + 2.0*sq(D12));
    D2c[] = f[]*D2;    // Dissipation rate in the liquid phase
  }

  // Adapt mesh based on multiple criteria:
  // - Interface position (f)
  // - Interface curvature (KAPPA)
  // - Velocity field components (u.x, u.y)
  // - Dissipation rate (D2c)
  adapt_wavelet((scalar *){f, KAPPA, u.x, u.y, D2c},
               (double[]){fErr, KErr, VelErr, VelErr, DissErr},
               MAXlevel, MINlevel);

  // Prevent unnecessary refinement far from the interface
  unrefine(x > 0.8*Ldomain && f[] < 1e-3);
}

// ======= Output snapshots at regular intervals =======
event writingFiles(t = 0, t += tsnap; t <= tmax) {
  p.nodump = false;    // Include pressure in output for post-processing
  dump(file = "dump"); // Save state for possible restart

  // Save numbered snapshots for visualization
  char nameOut[80];
  sprintf(nameOut, "intermediate/snapshot-%5.4f", t);
  dump(file = nameOut);
}

// ======= Log simulation data =======
event logWriting(i += 10) {
  // Calculate kinetic energy and other diagnostics
  double ke = 0., vol_liquid = 0., max_height = 0.;

  foreach (reduction(+:ke) reduction(+:vol_liquid) reduction(max:max_height)){
    ke += 0.5*rho(f[])*(sq(u.x[]) + sq(u.y[]))*sq(Delta);
    vol_liquid += f[]*sq(Delta);
    if (f[] > 0.5 && x < 0.2) {  // Track maximum height near the plate
      max_height = max(max_height, y);
    }
  }

  static FILE * fp;

  if (pid() == 0){    // Execute only on main process for parallel runs
    if (i == 0) {
      // Initialize log file with header
      fprintf(ferr, "i dt t ke vol_liquid max_height\n");
      fp = fopen("log", "w");
      
      fprintf(fp, "Level %d tmax %g. Ca %g, Mu21 %g, Rho21 %g, Ldomain %g\n",
              MAXlevel, tmax, Ca, Mu21, Rho21, Ldomain);
      fprintf(fp, "i dt t ke vol_liquid max_height\n");
      fprintf(fp, "%d %g %g %g %g %g\n", i, dt, t, ke, vol_liquid, max_height);
      fclose(fp);
    } else {
      // Append data to log file
      fp = fopen("log", "a");
      fprintf(fp, "%d %g %g %g %g %g\n", i, dt, t, ke, vol_liquid, max_height);
      fclose(fp);
    }
    fprintf(ferr, "%d %g %g %g %g %g\n", i, dt, t, ke, vol_liquid, max_height);
  }
}
