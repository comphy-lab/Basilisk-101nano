/**
# Cold Drop Impact on a Heated Substrate

This simulation models the thermal-fluid dynamics of a cold liquid drop impacting a heated solid surface.
The simulation captures both the fluid dynamics (spreading, rebound) and the heat transfer processes
that occur when a cold drop comes into contact with a hot surface. This physical phenomenon is
relevant to numerous applications including spray cooling, thermal management systems, 
hot surface deposition, and Leidenfrost effects.

## Physical Processes Modeled

The simulation captures several coupled physical phenomena:
1. **Fluid dynamics**: Drop deformation, spreading, and possible rebound
2. **Interface dynamics**: Evolution of the liquid-air interface with surface tension
3. **Heat transfer**: Conduction from the heated substrate into the cold drop
4. **Thermal boundary layer**: Development of temperature gradients near the solid surface
5. **Potential thermal effects**: Temperature-dependent surface tension effects may emerge

## Dimensionless Parameters

The simulation uses dimensionless numbers to characterize the physical regime:
- **Weber number (We)**: Ratio of inertial to surface tension forces
  We = ρU²D/σ, where U is impact velocity and D is drop diameter
- **Ohnesorge numbers (Oh)**: Ratio of viscous to inertial and surface tension forces
  Oh = μ/√(ρσD), defined separately for the drop (Ohd) and surrounding fluid (Ohs)
- **Bond number (Bo)**: Ratio of gravitational to surface tension forces
  Bo = ρgD²/σ
- **Thermal parameters** (implicit in diffusivity ratios):
  Prandtl number (Pr = ν/α) and thermal diffusivity ratio (D1/D2)

## Governing Equations

The system solves the coupled equations:

1. Navier-Stokes equations (conservation of momentum):
   ρ(∂u/∂t + u·∇u) = -∇p + ∇·(μ∇u) + ρg + σκδₛn

2. Continuity equation (conservation of mass):
   ∇·u = 0

3. Heat transfer equation:
   ∂T/∂t + u·∇T = ∇·(D∇T)

4. Interface advection equation:
   ∂f/∂t + u·∇f = 0

Where ρ is density, μ is viscosity, D is thermal diffusivity, σ is surface tension,
κ is interface curvature, f is the volume fraction field, and T is temperature.

## Numerical Implementation

The simulation employs several advanced numerical techniques:
- **Volume-of-Fluid (VOF)** method for interface tracking
- **Filtered interface** treatment for stability in property transitions
- **Adaptive Mesh Refinement (AMR)** based on multiple criteria (interface position, 
  curvature, velocity gradients, and thermal gradients)
- **Conservative momentum advection** for improved accuracy in high-density-ratio flows
- **Reduced gravity** approach for efficiency

## Boundary Conditions

- Left boundary: Heated solid wall with constant heat flux (Neumann boundary condition)
- Right/top boundaries: Outflow conditions
- No-slip velocity condition at the solid wall
- Temperature gradient at wall: Fixed positive value (heating the initially cold drop)

## Initial Configuration

A cold liquid drop is initialized with:
- Position: Slightly above the heated substrate
- Initial velocity: Directed toward the substrate
- Initial temperature: Uniform cold temperature (T=0), contrasting with the heated substrate

To run this simulation: `CFLAGS=-DDISPLAY=-1 make 3-HeatedDropImpact.tst`
*/

// Phase labeling: 1 is drop, 2 is surrounding fluid (air)

// ======= Include necessary Basilisk modules =======
#include "axi.h"                // Axisymmetric coordinates
#include "navier-stokes/centered.h"  // NS solver with centered discretization
#define FILTERED 1              // Use filtered VOF advection
#include "two-phase-thermal.h"  // Two-phase flow with heat transfer
#include "navier-stokes/conserving.h"  // Conservative momentum advection
#include "tension.h"            // Surface tension model
#include "reduced.h"            // Reduced gravity approach

// ======= Numerical parameters for adaptivity =======
// Error tolerances for adaptive mesh refinement
#define fErr (1e-3)             // Error tolerance in Volume of Fluid (interface position)
#define KErr (1e-6)             // Error tolerance in curvature (KAPPA)
#define VelErr (1e-2)           // Error tolerances in velocity fields
#define DissErr (1e-5)          // Error tolerances in dissipation rate
#define TErr (1e-3)             // Error tolerance in temperature

// ======= Physical parameters =======
#define Rho21 (1e-3)            // Density ratio (air/water)
// Drop positioning parameters
#define SPdist (0.02)           // Distance parameter for drop placement
#define R2Drop(x,y) (sq(x - 1e0 - SPdist) + sq(y))  // Function to define circular drop shape

// ======= Boundary conditions =======
// Left boundary: Heated solid wall (no-slip and fixed temperature gradient)
u.t[left] = dirichlet(0.0);     // Tangential velocity = 0 (no-slip)
f[left] = dirichlet(0.0);       // Volume fraction = 0 (solid wall)

// Right boundary: outflow condition
u.n[right] = neumann(0.);       // Zero gradient for normal velocity
p[right] = dirichlet(0.0);      // Reference pressure = 0

// Top boundary: outflow condition
u.n[top] = neumann(0.);         // Zero gradient for normal velocity
p[top] = dirichlet(0.0);        // Reference pressure = 0

// Thermal boundary condition: constant heat flux from the hot substrate
T[left] = neumann(10.0);        // Fixed positive temperature gradient (heating the cold drop)

// ======= Global parameters =======
int MAXlevel;                   // Maximum refinement level
double tmax, We, Ohd, Ohs, Bo, Ldomain;
#define MINlevel 2              // Minimum refinement level
#define tsnap (0.01)            // Time interval for snapshot outputs

/**
## Main Function

Sets up the simulation parameters, domain, and physical properties.
*/
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
   * This section allows passing arguments from the command line in the format:
   * ./executable MAXlevel tmax We Ohd Ohs Bo Ldomain
   * 
   * It is commented out by default but can be enabled for parameter studies.
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

  // Log the simulation parameters for reference
  fprintf(ferr, "Level %d tmax %g. We %g, Ohd %3.2e, Ohs %3.2e, Bo %g, Lo %g\n", 
          MAXlevel, tmax, We, Ohd, Ohs, Bo, Ldomain);

  // ======= Set up the computational domain =======
  L0 = Ldomain;                 // Domain size
  X0 = 0.; Y0 = 0.;             // Domain origin (axisymmetric, y=0 is symmetry axis)
  init_grid(1 << (4));          // Start with a 16×16 base grid (2^4)

  // Create directory for intermediate results
  char comm[80];
  sprintf(comm, "mkdir -p intermediate");
  system(comm);

  // ======= Set physical properties =======
  // Note: All quantities are in dimensionless form
  
  // Drop phase (1) properties
  rho1 = 1.0;                   // Density of drop (reference value)
  mu1 = Ohd/sqrt(We);           // Viscosity derived from Ohnesorge number
  D1 = 1.0;                     // Thermal diffusivity of drop (reference value)

  // Surrounding phase (2) properties
  rho2 = Rho21;                 // Density ratio (air/water)
  mu2 = Ohs/sqrt(We);           // Viscosity derived from Ohnesorge number
  D2 = 1e-3;                    // Thermal diffusivity ratio (air/water ~ 10^-3)

  // Interface and body force properties
  f.sigma = 1.0/We;             // Surface tension coefficient from Weber number
  G.x = -Bo/We;                 // Gravity from Bond number (vertical direction)

  // Start the simulation
  run();
}

/**
## Initialization
Sets up the initial drop position, shape, velocity, and temperature field.
*/
event init(t = 0){
  if(!restore(file = "dump")){
    // Refine mesh near the drop interface for better resolution
    refine((R2Drop(x,y) < 1.05) && (level < MAXlevel));
    
    // Initialize the drop shape using volume fraction
    // Circle equation: (x-x₀)² + y² = R², where x₀ = 1+SPdist, R = 1
    fraction(f, 1. - R2Drop(x,y));
    
    // Set initial velocity and temperature fields
    foreach () {
      u.x[] = -1.0*f[];  // Initial velocity in x-direction (downward)
                         // Only the drop has initial velocity (-1.0)
      u.y[] = 0.0;       // No initial velocity in y-direction
      T[] = 0.0;         // Initial cold temperature (T=0) for both drop and surrounding air
    }
  }

  // Uncomment to enable auto-restart functionality
  // dump (file = "dump");
  // return 1;
}

/**
## Adaptive Mesh Refinement
Dynamically refines the mesh based on multiple criteria to efficiently 
capture important flow and thermal features.
*/
scalar KAPPA[], D2c[];  // Curvature field and dissipation field
event adapt(i++){
  // Calculate interface curvature for surface tension and refinement
  curvature(f, KAPPA);
  
  // Calculate local viscous dissipation rate as an additional refinement criterion
  foreach(){
    // Calculate velocity gradient components in cylindrical coordinates
    double D11 = (u.y[0,1] - u.y[0,-1])/(2*Delta);  // ∂u_y/∂y
    double D22 = (u.y[]/max(y,1e-20));              // u_y/y (azimuthal strain rate)
    double D33 = (u.x[1,0] - u.x[-1,0])/(2*Delta);  // ∂u_x/∂x
    double D13 = 0.5*((u.y[1,0] - u.y[-1,0] + u.x[0,1] - u.x[0,-1])/(2*Delta)); // Shear rate
    
    // Calculate dissipation rate (sum of squares of strain rates)
    // This identifies regions of high deformation energy
    double D2 = (sq(D11)+sq(D22)+sq(D33)+2.0*sq(D13));
    D2c[] = f[]*D2;  // Dissipation rate in the drop phase
  }
  
  // Refine mesh based on multiple criteria:
  // 1. Interface position (f) - captures the drop boundary
  // 2. Interface curvature (KAPPA) - captures small surface features 
  // 3. Velocity field (u.x, u.y) - captures flow structures
  // 4. Dissipation rate (D2c) - captures regions of high deformation
  // 5. Temperature field (T) - captures thermal gradients
  adapt_wavelet((scalar *){f, KAPPA, u.x, u.y, D2c, T},
               (double[]){fErr, KErr, VelErr, VelErr, DissErr, TErr},
               MAXlevel, MINlevel);
  
  // Reduce refinement near outflow boundary to improve efficiency
  unrefine(x > 0.95*Ldomain);
}

/**
## Output Generation
Saves simulation data at regular intervals for post-processing and visualization.
*/
event writingFiles(t = 0, t += tsnap; t <= tmax) {
  p.nodump = false;  // Include pressure in output for post-processing
  dump(file = "dump");  // Save state for possible restart
  
  // Save numbered snapshots for visualization
  char nameOut[80];
  sprintf(nameOut, "intermediate/snapshot-%5.4f", t);
  dump(file = nameOut);
}

/**
## Data Logging
Records simulation progress and calculates integrated quantities for monitoring.
*/
event logWriting(i += 10) {
  // Calculate kinetic energy of the system
  // For axisymmetric simulations, the volume element is 2πy·dxdy
  double ke = 0.;
  foreach (reduction(+:ke)){
    ke += 2*pi*y*(0.5*rho(f[])*(sq(u.x[]) + sq(u.y[])))*sq(Delta);
  }

  // Write log data to file and terminal
  static FILE * fp;

  if (pid() == 0){  // Execute only on main process for parallel runs
    if (i == 0) {
      // Initialize log file with header and parameters
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