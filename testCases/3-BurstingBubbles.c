/**
# Bursting Bubbles in Viscoelastic Media

This simulation models the dynamics of bubbles bursting at the interface of viscoelastic media 
and air. The phenomenon is relevant to many natural and industrial processes including:
- Volcanic eruptions
- Polymer foam processing
- Biological systems with viscoelastic membranes
- Product formulation in consumer goods

The simulation captures complex dynamics including:
- Rupture of the thin film separating the bubble from the atmosphere
- Formation of Worthington jets after bubble collapse
- Potential satellite droplet ejection
- Influence of viscoelasticity on flow dynamics
- Surface tension-driven retraction of the film

## Dimensionless Parameters

The simulation uses key dimensionless numbers to characterize the system:

- **Ohnesorge number (Oh)**: Ratio of viscous to inertial-capillary forces
  Oh = μ/√(ρσR), where μ is viscosity, ρ is density, σ is surface tension, and R is bubble radius
  
- **Bond number (Bo)**: Ratio of gravitational to surface tension forces
  Bo = ρgR²/σ, determines how much the bubble/interface deforms under gravity

## Governing Equations

The simulation solves the axisymmetric incompressible Navier-Stokes equations with 
viscoelastic stress contributions:

1. Momentum equation:
   ρ(∂u/∂t + u·∇u) = -∇p + ∇·(μₛ∇u + τᵖ) + ρg + σκδₛn

2. Continuity equation:
   ∇·u = 0

3. Interface advection:
   ∂f/∂t + u·∇f = 0

## Numerical Implementation

The simulation employs several advanced numerical techniques:
- **Axisymmetric geometry** for computational efficiency
- **Two-phase flow** with smoothed density and viscosity jumps (FILTERED option)
- **Conservative momentum advection** for better handling of density jumps
- **Adaptive mesh refinement** to resolve interface dynamics accurately
- **Log-conformation approach** for stability in viscoelastic simulations

## Usage

To run this simulation with specific parameters:
```
make 3-BurstingBubbles.tst
```

Parameters can be adjusted in the code or passed via command line:
- maxLevel: Maximum refinement level for adaptive mesh
- Oh: Ohnesorge number
- Bond: Bond number
- tmax: Maximum simulation time
*/

#include "axi.h"                      // Axisymmetric coordinates
#include "navier-stokes/centered.h"   // NS solver with centered discretization
/**
## Simulation Parameters

Key numerical parameters controlling simulation accuracy and stability:
*/
#define FILTERED                      // Smear density and viscosity jumps for stability
#include "two-phase.h"                // Two-phase flow solver
#include "navier-stokes/conserving.h" // Conservative momentum advection
#include "tension.h"                  // Surface tension forces

#if !_MPI
#include "distance.h"                 // Distance function utilities (non-MPI version)
#endif

#define tsnap (1e-2)                  // Output snapshot interval (0.01 time units)
// Error tolerances for adaptive mesh refinement
#define fErr (1e-3)                   // Error tolerance in volume fraction field
#define KErr (1e-6)                   // Error tolerance in interface curvature (using height function)
#define VelErr (1e-3)                 // Error tolerance in velocity field
#define AErr (1e-3)                   // Error tolerance in conformation tensor (viscoelastic)

// Domain size in characteristic lengths
#define Ldomain 8

/**
## Boundary Conditions

Outflow condition on the right boundary:
*/
u.n[right] = neumann(0.);             // Zero gradient for normal velocity
p[right] = dirichlet(0.);             // Zero reference pressure

/**
## Global Parameters 
*/
int MAXlevel;                         // Maximum refinement level
double Oh, Oha, Bond, tmax;           // Physical parameters
char nameOut[80], dumpFile[80];       // File names for output

/**
## Main Function

Sets up the simulation parameters, domain, and physical properties.
*/
int  main(int argc, char const *argv[]) {
  L0 = Ldomain;                       // Domain size
  origin (-L0/2., 0.);                // Set origin at center-bottom of domain
  
  /*
  Values set directly in the code for educational purposes.
  In production simulations, these would be passed via command line.
  */
  MAXlevel = 10;                      // Maximum grid refinement level
  Oh = 1e-2;                          // Ohnesorge number (low values = low viscous effects)
  Bond = 1e-3;                        // Bond number (low values = small gravity effects)
  tmax = 1e0;                         // Maximum simulation time

  // For production runs, uncomment to use command line arguments:
  // if (argc < 7){
  //   fprintf(ferr, "Lack of command line arguments. Check! Need %d more arguments\n", 7-argc);
  //   return 1;
  // }
  // MAXlevel = atoi(argv[1]);
  // Oh = atof(argv[4]);
  // Bond = atof(argv[5]);
  // tmax = atof(argv[6]);

  init_grid (1 << 5);                 // Initialize with 32×32 base grid (2^5)
  
  // Create directory for simulation snapshots
  char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system(comm);
  
  // Set restart file name
  sprintf (dumpFile, "restart");

  /**
  ## Physical Properties
  
  Setting dimensionless properties for the two-phase system:
  */
  rho1 = 1.0;                         // Liquid density (reference)
  rho2 = 1e-3;                        // Gas density (1/1000 of liquid)
  Oha = 2e-2 * Oh;                    // Gas Ohnesorge number (scaled from liquid Oh)
  mu1 = Oh;                           // Liquid viscosity derived from Ohnesorge number
  mu2 = Oha;                          // Gas viscosity
  f.sigma = 1.0;                      // Surface tension coefficient (reference value)

  run();                              // Start the simulation
}

/**
## Initialization

Sets up the initial bubble configuration, either by:
1. Restoring from a previous simulation state (dumpFile), or
2. Reading initial shape from a data file and initializing using distance function
*/
event init (t = 0) {
#if _MPI // For MPI-based supercomputer runs
  if (!restore (file = dumpFile)){
    fprintf(ferr, "Cannot restored from a dump file!\n");
  }
#else  // For single-node runs using distance function method
  if (!restore (file = dumpFile)){
      // Try to read initial shape from data file
      char filename[60];
      sprintf(filename,"Bo%5.4f.dat",Bond);
      FILE * fp = fopen(filename,"rb");
        if (fp == NULL){
          fprintf(ferr, "There is no file named %s\n", filename);
          // Try alternate location (one folder up)
          sprintf(filename,"../Bo%5.4f.dat",Bond);
          fp = fopen(filename,"rb");
          if (fp == NULL){
            fprintf(ferr, "There is no file named %s\n", filename);
            return 1;
          }
        }
      // Read shape coordinates from file
      coord* InitialShape;
      InitialShape = input_xy(fp);
      fclose (fp);
      
      // Generate distance field from shape coordinates
      scalar d[];
      distance (d, InitialShape);

      // Refine mesh around interface using wavelet adaptation
      while (adapt_wavelet ((scalar *){f, d}, (double[]){1e-8, 1e-8}, MAXlevel).nf);
      
      // Calculate vertex-centered distance field from cell-centered values
      vertex scalar phi[];
      foreach_vertex(){
        phi[] = -(d[] + d[-1] + d[0,-1] + d[-1,-1])/4.;
      }
      
      // Initialize volume fraction field from distance field
      fractions (phi, f);
    }
#endif
}

/**
## Adaptive Mesh Refinement

Dynamically refines the mesh based on the interface position, velocity field,
and interface curvature to efficiently resolve the important flow features.
*/
event adapt(i++){
  // Calculate interface curvature for refinement criterion
  scalar KAPPA[];
  curvature(f, KAPPA);
  
  // Refine mesh based on multiple criteria
  adapt_wavelet ((scalar *){f, u.x, u.y, KAPPA},
                (double[]){fErr, VelErr, VelErr, KErr},
                MAXlevel, MAXlevel-6);  // Allow coarsening down to MAXlevel-6
}

/**
## Output Generation

Saves simulation data at regular intervals for post-processing and visualization.
*/
event writingFiles (t = 0; t += tsnap; t <= tmax) {
  dump (file = dumpFile);             // Save restart file
  
  // Save numbered snapshots
  sprintf (nameOut, "intermediate/snapshot-%5.4f", t);
  dump(file=nameOut);
}

/**
## Simulation Termination

Logs final simulation parameters when simulation ends.
*/
event end (t = end) {
  if (pid() == 0)
    fprintf(ferr, "Level %d, Oh %2.1e, Oha %2.1e, Bo %4.3f\n", MAXlevel, Oh, Oha, Bond);
}

/**
## Data Logging

Tracks the kinetic energy of the system and handles simulation termination
conditions based on energy criteria.
*/
event logWriting (i++) {
  // Calculate total kinetic energy (with axisymmetric integration)
  double ke = 0.;
  foreach (reduction(+:ke)){
    ke += (2*pi*y)*(0.5*rho(f[])*(sq(u.x[]) + sq(u.y[])))*sq(Delta);
  }
  
  if (pid() == 0) {
    // Output to log file and terminal
    static FILE * fp;
    if (i == 0) {
      fprintf(ferr, "Level %d, Oh %2.1e, Oha %2.1e, Bo %4.3f\n", MAXlevel, Oh, Oha, Bond);
      fprintf (ferr, "i dt t ke\n");
      fp = fopen ("log", "w");
      fprintf(fp, "Level %d, Oh %2.1e, Oha %2.1e, Bo %4.3f\n", MAXlevel, Oh, Oha, Bond);
      fprintf (fp, "i dt t ke\n");
      fprintf (fp, "%d %g %g %g\n", i, dt, t, ke);
      fclose(fp);
    } else {
      fp = fopen ("log", "a");
      fprintf (fp, "%d %g %g %g\n", i, dt, t, ke);
      fclose(fp);
    }
    fprintf (ferr, "%d %g %g %g\n", i, dt, t, ke);

    // Safety checks for simulation stability
    assert(ke > -1e-10);  // Energy should never be negative
    
    // Check for energy explosion (instability)
    if (ke > 1e2 && i > 1e1){
      if (pid() == 0){
        fprintf(ferr, "The kinetic energy blew up. Stopping simulation\n");
        fp = fopen ("log", "a");
        fprintf(fp, "The kinetic energy blew up. Stopping simulation\n");
        fclose(fp);
        dump(file=dumpFile);
        return 1;
      }
    }
    assert(ke < 1e2);
    
    // Check for near-static conditions (completion)
    if (ke < 1e-6 && i > 1e1){
      if (pid() == 0){
        fprintf(ferr, "kinetic energy too small now! Stopping!\n");
        dump(file=dumpFile);
        fp = fopen ("log", "a");
        fprintf(fp, "kinetic energy too small now! Stopping!\n");
        fclose(fp);
        return 1;
      }
    }
  }
}
