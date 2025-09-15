/**
 * @file 3-DropImpactOnSolids.c
 * @brief Simulation of a liquid drop impacting a solid surface
 * 
 * This code simulates the axisymmetric impact of a liquid drop on a solid surface.
 * The simulation tracks the interface between two phases (liquid drop and surrounding air)
 * and captures the dynamics of spreading and possible rebound.
 * 
 * Physical parameters:
 * - Weber number (We): Inertial vs surface tension forces
 * - Ohnesorge numbers (Ohd, Ohs): Viscous vs inertial and surface tension forces
 * - Bond number (Bo): Gravitational vs surface tension forces
 * 
 * @author Vatsal Sanjay (vatsal.sanjay@comphy-lab.org)
 * https://comphy-lab.org
 * @date 2025-03-11
 * @version 1.0 
*/

// Phase labeling: 1 is drop, 2 is surrounding fluid (air)

// ======= Include necessary Basilisk modules =======
#include "axi.h"                // Axisymmetric coordinates
#include "navier-stokes/centered.h"  // NS solver with centered discretization
#define FILTERED 1              // Use filtered VOF advection
#include "two-phase.h"          // Two-phase interface tracking
#include "navier-stokes/conserving.h"  // Conservative momentum advection
#include "tension.h"            // Surface tension model
#include "reduced.h"            // Reduced gravity approach

// ======= Numerical parameters for adaptivity =======
// Error tolerances for adaptive mesh refinement
#define fErr (1e-3)             // Error tolerance in Volume of Fluid (interface position)
#define KErr (1e-6)             // Error tolerance in curvature (KAPPA)
#define VelErr (1e-2)           // Error tolerances in velocity fields
#define DissErr (1e-5)          // Error tolerances in dissipation rate

// ======= Physical parameters =======
#define Rho21 (1e-3)            // Density ratio (air/water)
// Drop positioning parameters
#define SPdist (0.02)           // Distance parameter for drop placement
#define R2Drop(x,y) (sq(x - 1e0 - SPdist) + sq(y))  // Function to define circular drop shape

// ======= Boundary conditions =======
// Left boundary: solid wall (no-slip and no flux)
u.t[left] = dirichlet(0.0);     // Tangential velocity = 0 (no-slip)
f[left] = dirichlet(0.0);       // Volume fraction = 0 (solid wall)

// Right boundary: outflow condition
u.n[right] = neumann(0.);       // Zero gradient for normal velocity
p[right] = dirichlet(0.0);      // Reference pressure = 0

// Top boundary: outflow condition
u.n[top] = neumann(0.);         // Zero gradient for normal velocity
p[top] = dirichlet(0.0);        // Reference pressure = 0

// ======= Global parameters =======
int MAXlevel;                   // Maximum refinement level
double tmax, We, Ohd, Ohs, Bo, Ldomain;
#define MINlevel 2              // Minimum refinement level
#define tsnap (0.01)            // Time interval for snapshot outputs

int main(int argc, char const *argv[]) {

  // Default parameter values
  MAXlevel = 8;                 // Maximum grid refinement level
  tmax = 2.0;                   // Maximum simulation time
  We = 1.0;                     // Weber number (We = ρU²D/σ)
  Ohd = 1e-2;                   // Ohnesorge number for drop (Oh = μ/√(ρσD))
  Ohs = 1e-4;                   // Ohnesorge number for surrounding fluid
  Bo = 0.5;                     // Bond number (Bo = ρgD²/σ)
  Ldomain = 4.0;                // Domain size (dimensionless)

  /**
   * Parameter input from command line
   * Uncomment the following to pass the arguments from the command line:
   * ./executable MAXlevel tmax We Ohd Ohs Bo Ldomain
   */
  // if (argc < 8){
  //   fprintf(ferr, "Lack of command line arguments. Check! Need %d more arguments\n",8-argc);
  //   return 1;
  // }
  // MAXlevel = atoi(argv[1]);
  // tmax = atof(argv[2]);
  // We = atof(argv[3]); // We is 1 for 0.22 m/s <1250*0.22^2*0.001/0.06>
  // Ohd = atof(argv[4]); // <\mu/sqrt(1250*0.060*0.001)>
  // Ohs = atof(argv[5]); //\mu_r * Ohd
  // Bo = atof(argv[6]);
  // Ldomain = atof(argv[7]); // size of domain. must keep Ldomain \gg 1

  // Log the simulation parameters
  fprintf(ferr, "Level %d tmax %g. We %g, Ohd %3.2e, Ohs %3.2e, Bo %g, Lo %g\n", 
          MAXlevel, tmax, We, Ohd, Ohs, Bo, Ldomain);

  // ======= Set up the computational domain =======
  L0 = Ldomain;                 // Domain size
  X0 = 0.; Y0 = 0.;             // Domain origin
  init_grid(1 << (4));          // Start with a 16×16 base grid (2^4)

  // Create directory for intermediate results
  char comm[80];
  sprintf(comm, "mkdir -p intermediate");
  system(comm);

  // ======= Set physical properties =======
  // Note: All quantities are in dimensionless form
  rho1 = 1.0;                   // Density of drop (normalized)
  mu1 = Ohd/sqrt(We);           // Viscosity of drop derived from Oh and We
  rho2 = Rho21;                 // Density of surrounding fluid
  mu2 = Ohs/sqrt(We);           // Viscosity of surrounding fluid
  f.sigma = 1.0/We;             // Surface tension coefficient derived from We
  G.x = -Bo/We;                 // Gravity in x-direction (from Bond number)

  run();
}

// ======= Initialize the simulation =======
event init(t = 0){
  if(!restore(file = "dump")){
    // Refine mesh near the drop interface
    refine((R2Drop(x,y) < 1.05) && (level < MAXlevel));
    
    // Initialize the drop shape using volume fraction
    // f=1 inside drop, f=0 outside
    fraction(f, 1. - R2Drop(x,y));
    
    // Set initial velocity field (drop approaching the surface)
    foreach () {
      u.x[] = -1.0*f[];  // Initial velocity in x-direction (toward wall)
                         // Only the drop has initial velocity
      u.y[] = 0.0;       // No initial velocity in y-direction
    }
  }

  // Checkpoint for restarting the simulation
  // Uncomment to enable auto-restart functionality
  // dump (file = "dump");
  // return 1;
}

// ======= Adaptive mesh refinement =======
scalar KAPPA[], D2c[];  // Curvature field and dissipation field
event adapt(i++){
  // Calculate curvature for interface refinement
  curvature(f, KAPPA);
  
  // Calculate local dissipation rate as refinement criteria
  foreach(){
    // Calculate velocity gradient components in cylindrical coordinates
    double D11 = (u.y[0,1] - u.y[0,-1])/(2*Delta);  // ∂u_y/∂y
    double D22 = (u.y[]/max(y,1e-20));              // u_y/y (azimuthal strain rate)
    double D33 = (u.x[1,0] - u.x[-1,0])/(2*Delta);  // ∂u_x/∂x
    double D13 = 0.5*((u.y[1,0] - u.y[-1,0] + u.x[0,1] - u.x[0,-1])/(2*Delta)); // Shear rate
    
    // Calculate dissipation rate (sum of squares of strain rates)
    double D2 = (sq(D11)+sq(D22)+sq(D33)+2.0*sq(D13));
    D2c[] = f[]*D2;  // Dissipation rate in the drop phase
  }
  
  // Adapt mesh based on multiple criteria:
  // - Interface position (f)
  // - Interface curvature (KAPPA)
  // - Velocity field components (u.x, u.y)
  // - Dissipation rate (D2c)
  adapt_wavelet((scalar *){f, KAPPA, u.x, u.y, D2c},
               (double[]){fErr, KErr, VelErr, VelErr, DissErr},
               MAXlevel, MINlevel);
  
  // Prevent unnecessary refinement near outflow boundary
  // This helps reduce computational cost
  unrefine(x > 0.95*Ldomain);
}

// ======= Output snapshots at regular intervals =======
event writingFiles(t = 0, t += tsnap; t <= tmax) {
  p.nodump = false;  // Include pressure in output for post-processing
  dump(file = "dump");  // Save state for possible restart
  
  // Save numbered snapshots for visualization
  char nameOut[80];
  sprintf(nameOut, "intermediate/snapshot-%5.4f", t);
  dump(file = nameOut);
  output_ppm (f, linear = true, spread = 2, file = "f.mp4", n = 200);
}

// ======= Log simulation data =======
event logWriting(i += 10) {
  // Calculate kinetic energy of the system
  // For axisymmetric simulations, we integrate 2πy(...) to account for volume element
  double ke = 0.;
  foreach (reduction(+:ke)){
    ke += 2*pi*y*(0.5*rho(f[])*(sq(u.x[]) + sq(u.y[])))*sq(Delta);
  }

  static FILE * fp;

  if (pid() == 0){  // Execute only on main process for parallel runs
    if (i == 0) {
      // Initialize log file with header
      fprintf(ferr, "i dt t ke p\n");
      fp = fopen("log", "w");
      fprintf(fp, "Level %d tmax %g. We %g, Ohd %3.2e, Ohs %3.2e, Bo %g\n", 
              MAXlevel, tmax, We, Ohd, Ohs, Bo);
      fprintf(fp, "i dt t ke\n");
      fprintf(fp, "%d %g %g %g\n", i, dt, t, ke);
      fclose(fp);
    } else {
      // Append data to log file
      fp = fopen("log", "a");
      fprintf(fp, "%d %g %g %g\n", i, dt, t, ke);
      fclose(fp);
    }
    fprintf(ferr, "%d %g %g %g\n", i, dt, t, ke);
  }
}