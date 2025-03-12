/**
# Drop Impact on a Liquid Pool

This simulation models the axisymmetric impact of a liquid drop onto a pool of the 
same liquid. The phenomenon is common in nature (e.g., raindrops falling into puddles)
and industrial processes (e.g., spray cooling, inkjet printing, atomization).

The simulation captures the complex dynamics including:
- Initial impact and crater formation
- Capillary waves propagating across the pool surface
- Potential crown formation and splashing
- Possible bubble entrapment beneath the surface
- Eventual coalescence and return to equilibrium

## Dimensionless Parameters

The simulation uses three key dimensionless numbers to characterize the physical regime:

- **Froude number (Fr)**: Ratio of inertial to gravitational forces
  Fr = U²/(gD), where U is impact velocity, g is gravity, and D is drop diameter
  
- **Galilei number (Ga)**: Ratio of gravitational to viscous forces
  Ga = ρgD³/μ², characterizes the role of viscosity
  
- **Bond number (Bo)**: Ratio of gravitational to surface tension forces
  Bo = ρgD²/σ, determines how much the droplet/pool deforms under gravity

Additionally, density and viscosity ratios between the liquid and surrounding air
are set with RHO21 and MU21.

## Governing Equations

The simulation solves the incompressible Navier-Stokes equations with variable density
and surface tension:

1. Momentum equation:
   ρ(∂u/∂t + u·∇u) = -∇p + ∇·(μ∇u) + ρg + σκδₛn

2. Continuity equation:
   ∇·u = 0

3. Interface advection:
   ∂f/∂t + u·∇f = 0

## Numerical Implementation

The simulation employs several advanced numerical techniques:
- **Axisymmetric geometry** for computational efficiency
- **Volume-of-Fluid (VOF)** method for interface tracking
- **Adaptive mesh refinement** to resolve impact dynamics accurately
- **Conservative momentum advection** for better handling of density jumps
- **Tracer field** to distinguish between drop and pool liquid
- **Reduced gravity formulation** for better pressure handling

## Initial Configuration

A spherical drop is initialized above a flat liquid pool, with an initial downward
velocity corresponding to the specified Froude number (Fr). Both fluids are surrounded
by a low-density air phase.

## Usage

To run this simulation: `make 3-DropImpactOnPool.tst`
*/

#include "axi.h"                      // Axisymmetric coordinates
#include "navier-stokes/centered.h"   // NS solver with centered discretization
#define FILTERED 1                    // Enable filtered VOF for stability
#include "two-phase.h"                // Two-phase flow solver
#include "navier-stokes/conserving.h" // Conservative momentum advection
#include "tension.h"                  // Surface tension forces
#include "reduced.h"                  // Reduced gravity formulation
#include "distance.h"                 // Distance function utilities

// ======= Numerical parameters =======
#define tsnap (0.0025)                // Time interval for output snapshots

// Error tolerances for adaptive mesh refinement
#define fErr (1e-3)                   // Error tolerance for volume fraction
#define KErr (1e-4)                   // Error tolerance for interface curvature
#define VelErr (1e-4)                 // Error tolerance for velocity field

// ======= Physical parameters =======
#define Hint 1.5                      // Initial height of droplet center above pool
// Phase property ratios
#define RHO21 (1e-3)                  // Density ratio (air/liquid)
#define MU21 (1e-2)                   // Viscosity ratio (air/liquid)
// Dimensionless numbers
#define Ga (1e4)                      // Galilei number (gravitational/viscous forces)
#define Bo (1e1)                      // Bond number (gravitational/surface tension forces)
#define Ldomain 8                     // Domain size (in drop diameters)

// ======= Boundary conditions =======
// Right boundary is outflow condition
u.n[right] = neumann(0.);             // Zero gradient for normal velocity
p[right] = dirichlet(0.);             // Zero reference pressure (note: p is reduced pressure)

// ======= Global parameters =======
int MAXlevel;                         // Maximum refinement level
double tmax, Fr;                      // Max simulation time and Froude number

/**
## Main Function

Sets up the simulation parameters, domain, and physical properties.
*/
int main(int argc, char const *argv[]) {
  dtmax = 1e-3;                       // Maximum allowable timestep
  L0 = Ldomain;                       // Domain size
  origin(-L0/2., 0.);                 // Set origin (pool surface at y=0, centered at x=0)
  init_grid(1 << (6));                // Start with a 64×64 base grid (2^6)

  // Simulation parameters
  MAXlevel = 8;                       // Maximum grid refinement level
  tmax = 4e0;                         // Maximum simulation time
  Fr = 1e0;                           // Froude number (inertial/gravitational forces)
  
  // Log the simulation parameters
  fprintf(ferr, "Level %d, tmax %f, Fr %f, Ga %g, Bo %g\n", MAXlevel, tmax, Fr, Ga, Bo);

  // Create directory for intermediate results
  char comm[80];
  sprintf(comm, "mkdir -p intermediate");
  system(comm);

  // Set physical properties (dimensionless)
  rho1 = 1.0;                         // Liquid density (reference value)
  rho2 = RHO21;                       // Air density
  mu1 = 1./sqrt(Ga);                  // Liquid viscosity derived from Galilei number
  mu2 = MU21/sqrt(Ga);                // Air viscosity
  f.sigma = 1./Bo;                    // Surface tension coefficient from Bond number
  G.x = -1.;                          // Gravitational acceleration in x-direction

  run();                              // Start the simulation
}

/**
## Initialization

Sets up the initial drop and pool configuration. A tracer field is used to
distinguish between the drop and pool liquid, which allows tracking the drop's
evolution even after it merges with the pool.
*/
scalar tagDrop[];                     // Tracer to identify the drop liquid
event init(t = 0) {
  if (!restore(file = "dump")){
    // Refine mesh around the drop for better resolution
    refine(sq(x-Hint) + sq(y) < 1.1 && level < MAXlevel);
    
    // Initialize the fluid configuration:
    // - Combined shape of drop and pool using union of:
    //   - Drop: circle centered at (Hint,0) with radius 1
    //   - Pool: half-plane where x < 0
    fraction(f, union(1. - sq(x-Hint) - sq(y), -x));
    
    // Initialize the tracer to track only the drop (not the pool)
    fraction(tagDrop, 1 - sq(x-Hint) - sq(y));
    
    // Add the drop tracer to the list of fields advected with f
    // Using the correct method to handle the tracer
    f.tracers = list_append(f.tracers, tagDrop);

    // Set initial downward velocity for the drop only
    foreach(){
      u.x[] = -tagDrop[]*sqrt(Fr);    // Velocity proportional to √Fr
    }
  } else {
    // If restoring from a dump file, make sure to add the tracer
    f.tracers = list_append(f.tracers, tagDrop);
  }
}

/**
## Adaptive Mesh Refinement

Dynamically refines the mesh based on the interface position, curvature, and
momentum field to efficiently capture the important flow features during impact.
*/
event adapt(i++){
  scalar KAPPA[], ux[], uy[];         // Fields for refinement criteria
  
  // Calculate interface curvature for surface tension and refinement
  curvature(f, KAPPA);
  
  // Calculate momentum fields (ρu) for refinement criteria
  foreach(){
    ux[] = rho(f[])*u.x[];            // x-momentum
    uy[] = rho(f[])*u.y[];            // y-momentum
  }
  
  // Refine mesh based on multiple criteria:
  // 1. Interface position (f)
  // 2. Momentum fields (ux, uy)
  // 3. Interface curvature (KAPPA)
  adapt_wavelet((scalar *){f, ux, uy, KAPPA},
                (double[]){fErr, VelErr, VelErr, KErr},
                MAXlevel);
}

/**
## Output Generation

Saves simulation data at regular intervals for post-processing and visualization.
*/
event writingFiles(t = 0; t += tsnap; t <= tmax) {
  dump(file = "dump");                // Save full state for possible restart
  
  // Save numbered snapshots for visualization
  char nameOut[60];
  sprintf(nameOut, "intermediate/snapshot-%5.4f", t);
  dump(file=nameOut);
}

/**
## Data Logging

Tracks the kinetic energy of the liquid phase over time, which helps monitor
the energy transfer during impact and subsequent oscillations.
*/
event logWriting(i++) {
  // Calculate kinetic energy of the liquid phase
  // For axisymmetric simulations, we integrate 2πy(...) to account for volume element
  double ke = 0.;
  foreach(reduction(+:ke)){
    ke += (2*pi*y)*(0.5*(f[])*rho1*(sq(u.x[]) + sq(u.y[])))*sq(Delta);
  }
  
  // Write to log file and terminal
  static FILE * fp;
  if (i == 0) {
    fprintf(ferr, "i dt t ke\n");
    fp = fopen("log", "w");
    fprintf(fp, "i dt t ke\n");
    fprintf(fp, "%d %g %g %g\n", i, dt, t, ke);
    fclose(fp);
  } else {
    fp = fopen("log", "a");
    fprintf(fp, "%d %g %g %g\n", i, dt, t, ke);
    fclose(fp);
  }
  fprintf(ferr, "%d %g %g %g\n", i, dt, t, ke);
}
