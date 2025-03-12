/**
 * # Bursting Bubbles in Viscoelastic Media
 * @file burstingBubbleVE.c
 * @brief Simulation of bursting bubbles in viscoelastic media using the Basilisk framework
 * @author Vatsal Sanjay and Ayush Dixit
 * @version 1.0
 * @date Nov 23, 2024
 * 
 * ## Overview
 * This simulation demonstrates the complex dynamics of bubbles bursting in viscoelastic fluids
 * (like polymer solutions). When bubbles burst at interfaces, they can create spectacular 
 * "Worthington jets" that shoot upward and sometimes break into droplets. In viscoelastic media,
 * these dynamics are dramatically altered by the fluid's memory effects.
 * 
 * This simulation is a perfect example of multiphase flows with complex rheology, showcasing:
 * - Two-phase flow modeling (liquid polymer + air)
 * - Viscoelasticity through log-conformation approach
 * - Adaptive mesh refinement for capturing interfacial dynamics
 * - Axisymmetric formulation for computational efficiency
 * 
 * ## Physics Background
 * Viscoelastic fluids exhibit both viscous (fluid-like) and elastic (solid-like) properties.
 * The dimensionless numbers controlling this simulation include:
 * - **Deborah number (De)**: Ratio of relaxation time to flow time. Higher values mean more elastic behavior.
 * - **Elasto-capillary number (Ec)**: Ratio of elastic forces to surface tension.
 * - **Ohnesorge number (Oh)**: Ratio of viscous forces to inertial and capillary forces.
 * - **Bond number (Bond)**: Ratio of gravitational forces to surface tension.
 * 
 * ## Usage
 * ```
 * ./program maxLevel De Ec Oh Bond tmax
 * ```
 * where:
 *   - **maxLevel**: Maximum refinement level for adaptive mesh (8-12 recommended)
 *   - **De**: Deborah number (try 0.1-10 for interesting behaviors)
 *   - **Ec**: Elasto-capillary number (0.001-0.1)
 *   - **Oh**: Ohnesorge number (1e-3 to 1 covers many polymer solutions)
 *   - **Bond**: Bond number (controls initial bubble shape)
 *   - **tmax**: Maximum simulation time (1-10 is typically sufficient)
 * 
 * ## Learning Objectives
 * - Understand how to set up multiphase simulations with complex rheology
 * - Learn to implement and handle non-Newtonian fluid behavior
 * - Analyze the effects of dimensionless parameters on flow dynamics
 * - Practice with adaptive mesh refinement for interface problems
 * 
 * ## EXERCISES
 * 1. Try varying De while keeping other parameters fixed. What happens to the jet as elasticity increases?
 * 2. Compare Newtonian (De=0, Ec=0) vs viscoelastic behavior. What are the key differences?
 * 3. Investigate mesh dependency by changing maxLevel. When does the solution converge?
 * 4. Modify the code to visualize the stress tensor components during jet formation
 * 
 * ## References
 * - Sanjay, V. (2024). Basilisk Implementation of Viscoelastic Flows. Zenodo, DOI: 10.5281/zenodo.14210635
 * - Worthington, A.M. (1908). A Study of Splashes. Longmans, Green and Company
 * 
*/

#include "axi.h"                  // Axisymmetric formulation (r,z) coordinates
#include "navier-stokes/centered.h"  // Main Navier-Stokes solver

/**
 * ## Viscoelastic Model Implementation
 * We use the log-conformation approach for numerical stability when simulating viscoelastic fluids.
 * This method transforms the conformation tensor logarithmically to preserve positive-definiteness
 * during the numerical solution.
 * 
 * Reference: V. Sanjay, Zenodo, DOI: 10.5281/zenodo.14210635 (2024)
 */
// #define _SCALAR // uncomment to use the scalar version of the viscoelastic code
#if !_SCALAR
#include "log-conform-viscoelastic.h"  // Tensor implementation (full anisotropic stresses)
#else 
#include "log-conform-viscoelastic-scalar-2D.h"  // Scalar implementation (faster but less accurate)
#endif

/**
 * ## Simulation Parameters
 * These parameters control the numerical accuracy and computational demands.
 * Adjust them based on your computational resources and required accuracy.
 * 
 * - **FILTERED**: Enable density and viscosity jump smoothing for numerical stability
 * - **tsnap**: Time interval between snapshots (default: 1e-2)
 * - **fErr**: Error tolerance for volume fraction tracking (1e-3)
 * - **KErr**: Error tolerance for curvature calculation (1e-6)
 * - **VelErr**: Error tolerance for velocity field (1e-3)
 * - **AErr**: Error tolerance for conformation tensor (1e-3)
 * - **Ldomain**: Domain size in characteristic lengths (8)
 */
#define FILTERED                           // Smear density and viscosity jumps for stability
#include "two-phaseVE.h"                   // Two-phase flow with viscoelasticity
#include "navier-stokes/conserving.h"      // Conservative form of N-S for better mass conservation
#include "tension.h"                       // Surface tension implementation

#if !_MPI
#include "distance.h"                      // Distance function (not compatible with MPI)
#endif

#define tsnap (1e-2)                       // Time interval between output snapshots
// Error tolerances for adaptive mesh refinement
#define fErr (1e-3)                        // Error tolerance in volume fraction f
#define KErr (1e-6)                        // Error tolerance in curvature calculation
#define VelErr (1e-3)                      // Error tolerances in velocity field
#define AErr (1e-3)                        // Error tolerances in conformation tensor

// Domain size
#define Ldomain 8

// Boundary conditions - outflow on the right boundary
u.n[right] = neumann(0.);                  // Zero normal gradient for velocity at right boundary
p[right] = dirichlet(0.);                  // Zero pressure at right boundary (reference pressure)

// Global parameters
int MAXlevel;                              // Maximum refinement level
double Oh, Oha, De, Ec, Bond, tmax;        // Dimensionless numbers
char nameOut[80], dumpFile[80];            // Output file names

/**
 * ## Main Function
 * Sets up the simulation parameters, initializes the domain, and runs the simulation.
 */
int main(int argc, char const *argv[]) {
  // Maximum allowed timestep (important for stability)
  dtmax = 1e-5;                            // BEWARE: Critical for numerical stability!

  // Set up domain size and origin
  L0 = Ldomain;                            // Set domain size
  origin (-L0/2., 0.);                     // Center domain at (-L0/2, 0) - needed for axisymmetric flows
  
  /**
   * ### Parameter Setup
   * In a real run, these values would be passed from command line.
   * For this example, we use fixed representative values.
   * 
   * EXPERIMENT: Try changing these parameters and observe the effect on jet dynamics!
   */
  MAXlevel = 10;                           // Maximum refinement level (computational cost ~ 4^MAXlevel)
  De = 0.1;                                // Deborah number (elasticity)
  Ec = 0.01;                               // Elasto-capillary number 
  Oh = 1e-2;                               // Ohnesorge number (viscosity)
  Bond = 1e-3;                             // Bond number (gravity)
  tmax = 1e0;                              // Maximum simulation time

  // To use command line arguments, uncomment these:
  // MAXlevel = atoi(argv[1]);
  // De = atof(argv[2]);  // Use a value of 1e30 to simulate the De \to \infty limit
  // Ec = atof(argv[3]);
  // Oh = atof(argv[4]);
  // Bond = atof(argv[5]);
  // tmax = atof(argv[6]);
  
  // Verify command line arguments (commented out for this example)
  // if (argc < 7){
  //   fprintf(ferr, "Lack of command line arguments. Check! Need %d more arguments\n", 7-argc);
  //   return 1;
  // }
  
  // Initialize grid with resolution 2^5 = 32 cells per dimension
  init_grid (1 << 5);
  
  // Create output directory and set up restart file name
  char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system(comm);
  sprintf (dumpFile, "restart");

/**
 * ### Physical Properties
 * Here we define the physical properties of both fluids in the simulation.
 * All properties are non-dimensionalized.
 * 
 * - **rho1, rho2**: Density ratio between liquid and gas (1:0.001)
 * - **mu1, mu2**: Viscosity ratio controlled by Ohnesorge numbers
 * - **lambda1, lambda2**: Relaxation times (polymer has memory, gas does not)
 * - **G1, G2**: Elastic moduli (polymer is elastic, gas is not)
 * - **f.sigma**: Surface tension coefficient (normalized to 1.0)
 * 
 * NOTE: For viscoelastic fluids, the total stress is:
 * σ = -pI + μs(∇u + ∇u^T) + G(A-I)
 * where A is the conformation tensor, μs is solvent viscosity, and G is elastic modulus.
 */
  // Density ratio (liquid : gas = 1000:1)
  rho1 = 1., rho2 = 1e-3;
  
  // Gas phase Ohnesorge number (typically lower than liquid)
  Oha = 2e-2 * Oh;
  
  // Viscosity from Ohnesorge numbers
  mu1 = Oh, mu2 = Oha;
  
  // Polymer relaxation time (liquid has elasticity, gas does not)
  lambda1 = De; lambda2 = 0.;
  
  // Elastic modulus (liquid is viscoelastic, gas is not)
  G1 = Ec; G2 = 0.;

  // Surface tension coefficient (normalized to 1.0)
  f.sigma = 1.0;

  // Numerical parameters for stability
  TOLERANCE = 1e-4;                        // Tolerance for iterative pressure solver
  CFL = 1e-1;                              // Courant-Friedrichs-Lewy number (stability criterion)

  // Start the simulation
  run();
}

/**
 * ## Initialization Event
 * Sets up the initial condition for the simulation at t=0.
 * 
 * This function either:
 * 1. Restores from a previous dump file, OR
 * 2. Initializes the bubble shape based on a pre-computed shape file
 * 
 * The shape file contains the interface coordinates for a bubble at equilibrium
 * under the given Bond number, calculated from a separate simulation.
 */
event init (t = 0) {
#if _MPI // Special handling for MPI parallel runs
  if (!restore (file = dumpFile)){
    fprintf(ferr, "Cannot restored from a dump file!\n");
  }
#else  // Serial run with distance function (not compatible with MPI)
  if (!restore (file = dumpFile)){
      // Try to load the initial bubble shape based on Bond number
      char filename[60];
      sprintf(filename,"Bo%5.4f.dat",Bond);
      FILE * fp = fopen(filename,"rb");
        if (fp == NULL){
          fprintf(ferr, "There is no file named %s\n", filename);
          // Try in folder one level up
          sprintf(filename,"../Bo%5.4f.dat",Bond);
          fp = fopen(filename,"rb");
          if (fp == NULL){
            fprintf(ferr, "There is no file named %s\n", filename);
            return 1;
          }
        }
      // Load the interface coordinates
      coord* InitialShape;
      InitialShape = input_xy(fp);
      fclose (fp);
      
      // Compute signed distance function from the interface
      scalar d[];
      distance (d, InitialShape);

      // Refine mesh based on distance function
      while (adapt_wavelet ((scalar *){f, d}, (double[]){1e-8, 1e-8}, MAXlevel).nf);
      
      // Convert distance function to vertex values for interface reconstruction
      vertex scalar phi[];
      foreach_vertex(){
        phi[] = -(d[] + d[-1] + d[0,-1] + d[-1,-1])/4.;
      }
      
      // Initialize volume fraction from distance function
      fractions (phi, f);
    }
#endif
}

/**
 * ## Adaptive Mesh Refinement
 * 
 * This event dynamically refines the mesh based on solution features.
 * It runs after each timestep (i++) to adapt the mesh according to error
 * estimators for multiple fields:
 * 
 * - Interface location (f)
 * - Velocity field (u.x, u.y)
 * - Conformation tensor components
 * - Interface curvature (KAPPA)
 * 
 * The mesh can refine up to MAXlevel in regions of interest and 
 * coarsen down to (MAXlevel-6) in smooth regions to save computation.
 * 
 * EXERCISE: Try changing the error tolerances and observe the effect on
 * simulation accuracy vs. computational cost.
 */
event adapt(i++){
  // Calculate interface curvature for refinement criterion
  scalar KAPPA[];
  curvature(f, KAPPA);

  // Adapt mesh based on solution gradients
  #if !_SCALAR
   // For tensor viscoelastic model, refine based on all tensor components
   adapt_wavelet ((scalar *){f, u.x, u.y, conform_p.x.x, conform_p.y.y, conform_p.y.x, conform_qq, KAPPA},
      (double[]){fErr, VelErr, VelErr, AErr, AErr, AErr, AErr, KErr},
      MAXlevel, MAXlevel-6);
  #else
   // For scalar viscoelastic model, refine based on scalar components
   adapt_wavelet ((scalar *){f, u.x, u.y, A11, A22, A12, AThTh, KAPPA},
      (double[]){fErr, VelErr, VelErr, AErr, AErr, AErr, AErr, KErr},
      MAXlevel, MAXlevel-6);
  #endif
}

/**
 * ## Simulation Output
 * 
 * This event saves snapshots of the simulation at regular intervals.
 * Two files are created:
 * 1. A restart file (updated continuously)
 * 2. Timestamped snapshots for visualization and analysis
 * 
 * You can use these dumps with Basilisk's view.c to visualize the results
 * or process them with custom analysis tools.
 * 
 * VISUALIZATION TIP: Try running:
 * qcc -Wall -O2 view.c -o view -lm
 * ./view snapshot-*.dat
 */
event writingFiles (t = 0; t += tsnap; t <= tmax) {
  dump (file = dumpFile);                  // Create/update restart file
  sprintf (nameOut, "intermediate/snapshot-%5.4f", t);  // Create snapshot files
  dump(file=nameOut);
}

/**
 * ## Simulation Completion
 * 
 * This event runs at the end of the simulation (t = end) and prints
 * a summary of the simulation parameters to standard error.
 */
event end (t = end) {
  if (pid() == 0)
    fprintf(ferr, "Level %d, De %2.1e, Ec %2.1e, Oh %2.1e, Oha %2.1e, Bo %4.3f\n", 
            MAXlevel, De, Ec, Oh, Oha, Bond);
}

/**
 * ## Simulation Logging
 * 
 * This event runs after each timestep to:
 * 1. Calculate the total kinetic energy in the system
 * 2. Log simulation progress and diagnostics
 * 3. Check for simulation stability and stop if necessary
 * 
 * The kinetic energy is a good indicator of simulation health:
 * - If it grows too large, the simulation is becoming unstable
 * - If it becomes too small, the flow has essentially stopped
 * 
 * IMPORTANT: For axisymmetric simulations, volume integrals must include
 * the 2πy factor in cylindrical coordinates.
 */
event logWriting (i++) {
  // Calculate total kinetic energy in the domain
  double ke = 0.;
  foreach (reduction(+:ke)){
    // For axisymmetric coordinates, include 2πy factor for volume element
    ke += (2*pi*y)*(0.5*rho(f[])*(sq(u.x[]) + sq(u.y[])))*sq(Delta);
  }
  
  // Write log information (only on main process in parallel runs)
  if (pid() == 0) {
    static FILE * fp;
    if (i == 0) {
      // Write header information at start
      fprintf(ferr, "Level %d, De %2.1e, Ec %2.1e, Oh %2.1e, Oha %2.1e, Bo %4.3f\n", 
              MAXlevel, De, Ec, Oh, Oha, Bond);
      fprintf (ferr, "De Ec Oh i dt t ke\n");
      fp = fopen ("log", "w");
      fprintf(fp, "Level %d, De %2.1e, Ec %2.1e, Oh %2.1e, Oha %2.1e, Bo %4.3f\n", 
              MAXlevel, De, Ec, Oh, Oha, Bond);
      fprintf (fp, "i dt t ke\n");
      fprintf (fp, "%d %g %g %g\n", i, dt, t, ke);
      fclose(fp);
    } else {
      // Append data in subsequent timesteps
      fp = fopen ("log", "a");
      fprintf (fp, "%d %g %g %g\n", i, dt, t, ke);
      fclose(fp);
    }
    fprintf (ferr, "%d %g %g %g\n", i, dt, t, ke);

    // Check for negative kinetic energy (physical impossibility, numerical error)
    assert(ke > -1e-10);

    // Check for kinetic energy blowup (simulation instability)
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
    
    // Check if simulation has essentially stopped
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
