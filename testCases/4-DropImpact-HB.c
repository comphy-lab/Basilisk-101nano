/**
 * # Viscoplastic Drop Impact on a Solid Surface
 * 
 * Version 2.0
 * Author: Vatsal Sanjay
 * Last updated: Oct 11, 2024
 * 
 * ## Overview
 * This simulation demonstrates the classical fluid dynamics problem of a viscoplastic (VP) 
 * drop impacting a solid surface. This phenomenon is relevant to numerous applications,
 * including:
 * 
 * - Inkjet printing with non-Newtonian inks
 * - Spray coating with yield-stress fluids
 * - Food processing (e.g., droplets of ketchup, mayonnaise)
 * - Blood splatter analysis
 * - Additive manufacturing with complex materials
 * 
 * ## Physics Background
 * When a viscoplastic drop impacts a surface, several competing forces determine its behavior:
 * 
 * - **Inertia**: Drives the drop to spread upon impact
 * - **Surface tension**: Works to minimize surface area
 * - **Viscous forces**: Dissipate energy, slowing down the spreading
 * - **Yield stress**: Prevents deformation until a threshold stress is exceeded
 * 
 * The interaction of these forces creates complex dynamics that differ significantly from
 * Newtonian drop impacts. For instance, viscoplastic drops may exhibit partial rebound,
 * asymmetric spreading, or solidification mid-deformation.
 * 
 * ## Numerical Implementation
 * In this simulation, we use a two-phase approach where:
 * - **Phase 1**: Viscoplastic liquid drop (using Herschel-Bulkley model)
 * - **Phase 2**: Newtonian gas (surrounding air)
 */

#include "axi.h"
#include "navier-stokes/centered.h"
#define FILTERED 1 // Smear density and viscosity jumps
/**
To model Viscoplastic liquids, we use a modified version of [two-phase.h](http://basilisk.fr/src/two-phase.h). [two-phaseVP-HB.h](../src-local/two-phaseVP-HB.h) contains these modifications.
*/
#include "two-phaseVP-HB.h"
/**
 You can use: conserving.h as well. Even without it, I was still able to conserve the total energy (also momentum?) of the system if I adapt based on curvature and vorticity/deformation tensor norm (see the adapt even). I had to smear the density and viscosity anyhow because of the sharp ratios in liquid (Bingham) and the gas.
*/
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "reduced.h"

#define tsnap (0.01)

// Error tolerancs
#define fErr (1e-3)                                 // error tolerance in f1 VOF
#define VelErr (1e-2)                               // error tolerances in velocity -- Use 1e-2 for low Oh and 1e-3 to 5e-3 for high Oh/moderate to high J
#define D2Err (1e-2)
#define KAPPAErr (1e-2)

// gas properties!
#define RHO21 (1e-3)
#define MU21 (1e-2)

// domain properties!
#define Ldomain 4

// Distance and radius of drop calculations
#define Xdist (1.02)
#define R2Drop(x,y) (sq(x - Xdist) + sq(y))

// boundary conditions
u.t[left] = dirichlet(0.0);
f[left] = dirichlet(0.0);
// all symmetry planes

int MAXlevel;
double We, Oh, J, tmax;
char nameOut[80], dumpFile[80];

int  main(int argc, char const *argv[]) {

  // Ensure that all the variables were transferred properly from the terminal or job script.
  // if (argc < 6){
  //   fprintf(ferr, "Lack of command line arguments. Check! Need %d more arguments, Level, tauy, We, Oh, tmax\n",6-argc);
  //   return 1;
  // }

  L0 = Ldomain;
  origin (0., 0.);
  init_grid (1 << 8);
  // Values taken from the terminal
  MAXlevel = 9; // atoi(argv[1]);
  J = 1e-1; // atof(argv[2]); // plasto-capillary number
  We = 1e1; // atof(argv[3]); // Weber number
  Oh = 1e-2; // atof(argv[4]); // Ohnesorge number
  tmax = 1e0; // atof(argv[5]);

  fprintf(ferr, "Level %d, We %2.1e, Oh %2.1e, J %4.3f\n", MAXlevel, We, Oh, J);

  // Create a folder named intermediate where all the simulation snapshots are stored.
  char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system(comm);
  // Name of the restart file. See writingFiles event.
  sprintf (dumpFile, "dump");

  /**
   * ## Mathematical Formulation
   * 
   * ### Governing Equations
   * 
   * We consider the drop impacting a solid surface, employing non-dimensionalization based on:
   * - Length scale: Initial drop radius $R_0$
   * - Velocity scale: Inertia-capillary velocity $V_\gamma = \sqrt{\gamma/(\rho_lR_0)}$
   * - Pressure scale: Capillary pressure, $\gamma/R_0$
   * 
   * The dimensionless mass and momentum conservation equations for the liquid phase are:
   * 
   * $$
   * \nabla \cdot \boldsymbol{u} = 0
   * $$
   * 
   * $$
   * \frac{\partial\boldsymbol{u}}{\partial t} + \nabla\boldsymbol{\cdot}\left(\boldsymbol{uu}\right) = -\nabla p + \nabla\boldsymbol{\cdot}\boldsymbol{\tau} - \mathcal{B}o\,\hat{\boldsymbol{e}}_{\boldsymbol{\mathcal{Z}}}
   * $$
   * 
   * where:
   * - $\boldsymbol{u}$ is the velocity vector
   * - $t$ is time
   * - $p$ is the pressure
   * - $\boldsymbol{\tau}$ represents the deviatoric stress tensor
   * 
   * ### Herschel-Bulkley Constitutive Equation
   * 
   * For the viscoplastic behavior, we employ the regularized Herschel-Bulkley model:
   * 
   * $$
   * \boldsymbol{\tau} = 2\biggl[\frac{\mathcal{J}}{2\|\boldsymbol{\mathcal{D}}\|+\epsilon} + \mathcal{O}h\,(2\|\boldsymbol{\mathcal{D}}\|+\epsilon)^{n-1}\biggr]\boldsymbol{\mathcal{D}}
   * $$
   * 
   * where:
   * - $\|\boldsymbol{\mathcal{D}}\|$ is the second invariant of the deformation rate tensor
   * - $\epsilon$ is the regularization parameter
   * - $n$ is the power law index (in this simulation, we use $n=1$ for Bingham model)
   * 
   * ### Dimensionless Parameters
   * 
   * The simulation is governed by several dimensionless groups:
   * 
   * 1. **Capillary-Bingham number** ($\mathcal{J}$): Ratio of yield stress to capillary pressure
   *    $$\mathcal{J} = \frac{\tau_yR_0}{\gamma}$$
   * 
   * 2. **Ohnesorge number** ($\mathcal{O}h$): Ratio of viscous forces to inertial and capillary forces
   *    $$\mathcal{O}h = \frac{\mu_l}{\sqrt{\rho_l\gamma R_0}}$$
   * 
   * 3. **Weber number** ($\mathcal{W}e$): Ratio of inertial to capillary forces
   *    $$\mathcal{W}e = \frac{\rho_l V^2 R_0}{\gamma}$$
   *    
   * 4. **Bond number** ($\mathcal{B}o$): Ratio of gravitational to capillary forces
   *    $$\mathcal{B}o = \frac{\rho_l gR_o^2}{\gamma}$$
   * 
   * In this simulation, we primarily vary $J$, $We$, and $Oh$ while keeping density and viscosity 
   * ratios between phases constant:
   * - Density ratio: $\rho_r = \rho_g/\rho_l = 10^{-3}$
   * - Viscosity ratio: $\mu_r = \mu_g/\mu_l = 2 \times 10^{-2}$
   */

  // epsilon = t < tsnap ? 1e-1 : 1e-3;  // epsilon regularisation value of effective viscosity
  epsilon = 1e-2;  // epsilon regularisation value of effective viscosity
  rho1 = 1., rho2 = RHO21;
  mu1 = Oh/sqrt(We), mu2 = MU21*Oh/sqrt(We);
  f.sigma = 1.0/We;
  tauy = J/We;
  CFL = 1e-1;
  run();
}

event init (t = 0) {
  if (!restore (file = dumpFile)){
    refine((R2Drop(x, y) < 1.05) && (level < MAXlevel));
    fraction(f, 1. - R2Drop(x, y));
    foreach() {
      u.x[] = -1.0 * f[];
      u.y[] = 0.0;
    }
  }
}

/**
 * ## Adaptive Mesh Refinement Strategy
 * 
 * The simulation employs a dynamic adaptive mesh refinement (AMR) approach that
 * intelligently allocates computational resources to critical regions. The strategy
 * changes based on the simulation time:
 * 
 * ### Initial Stage (t < 0.01)
 * During the early impact phase, refinement focuses on tracking:
 * - The interface between fluids (volume fraction field f)
 * - Velocity fields for accurate momentum capture
 * 
 * ### Later Stages (t ≥ 0.01)
 * As the simulation progresses, we introduce additional refinement criteria:
 * 
 * 1. **Interface Curvature (KAPPA)**: Ensures the interface shape is resolved accurately,
 *    particularly important for capturing surface tension effects
 * 
 * 2. **Deformation Tensor Magnitude (D2c)**: Concentrates resolution where the material
 *    is actively deforming, especially near the yield surfaces where viscoplastic
 *    behavior is most pronounced
 * 
 * This adaptive strategy significantly improves computational efficiency while
 * maintaining accuracy in the physically important regions of the flow.  
 */
event adapt(i++){
  if (t < 1e-2){
    adapt_wavelet ((scalar *){f, u.x, u.y},
    (double[]){fErr, VelErr, VelErr},
    MAXlevel);
  } else {
    scalar KAPPA[], D2c[];
    curvature(f, KAPPA);
    foreach() {
      double D11 = (u.y[0,1] - u.y[0,-1])/(2*Delta);
      double D22 = (u.y[]/y);
      double D33 = (u.x[1,0] - u.x[-1,0])/(2*Delta);
      double D13 = 0.5*( (u.y[1,0] - u.y[-1,0] + u.x[0,1] - u.x[0,-1])/(2*Delta) );
      double D2 = (sq(D11)+sq(D22)+sq(D33)+2.0*sq(D13));
      D2c[] = f[]*(D2);
    }
    adapt_wavelet ((scalar *){f, u.x, u.y, KAPPA, D2c},
    (double[]){fErr, VelErr, VelErr, KAPPAErr, D2Err},
    MAXlevel);
  }
}

/**
 * ## Output Data Management
 * 
 * The simulation generates two types of output files to enable comprehensive analysis:  
 * 
 * ### 1. Snapshots
 * Periodic simulation states are saved at regular intervals (controlled by `tsnap`),
 * allowing for detailed visualization and post-processing. These files contain complete
 * information about all fields (velocity, pressure, volume fraction, etc.) at each
 * specified time point.  
 * 
 * ### 2. Restart File
 * A single `dump` file is continuously updated, storing the latest state of the simulation.
 * This serves as both a checkpointing mechanism (allowing restart from interruptions) and
 * a backup of the final state.  
 * 
 * All outputs are organized in the `intermediate` directory for easy management and to
 * prevent clutter in the main simulation directory.  
 */
event writingFiles (t = 0; t += tsnap; t <= tmax+tsnap) {
  dump (file = dumpFile);
  sprintf (nameOut, "intermediate/snapshot-%5.4f", t);
  dump(file=nameOut);
}

/**
 * ## Simulation Termination
 * 
 * The simulation concludes when one of these conditions is met:  
 * 
 * 1. **Time Limit**: When the simulation reaches the specified maximum time (`tmax`).  
 * 2. **Energy Threshold**: When the kinetic energy drops below 10^-6 of its maximum value,
 *    indicating the system has essentially reached equilibrium.  
 * 
 * Upon completion, the code outputs a summary of the simulation parameters to the error
 * stream, providing essential information for post-processing and analysis.  
 */
event end (t = end) {
  fprintf(ferr, "Done: Level %d, We %2.1e, Oh %2.1e, J %4.3f\n", MAXlevel, We, Oh, J);
}

/**
 * ## Diagnostic Data Collection and Monitoring
 * 
 * The simulation continuously tracks and logs key diagnostic information:  
 * 
 * ### 1. Energy Conservation Monitoring  
 * - Calculates total kinetic energy at each timestep
 * - Writes values to a log file for post-processing analysis
 * - Outputs real-time energy values to the error stream for monitoring
 * 
 * ### 2. Simulation Integrity Checks  
 * - Verifies kinetic energy remains physically reasonable (non-negative and bounded)  
 * - Automatically terminates if energy falls below threshold (system at rest)  
 * - Prevents wasting computational resources on essentially completed simulations  
 * 
 * This diagnostic approach helps ensure physical correctness of the simulation and
 * provides valuable data for analyzing energy dissipation rates—a critical characteristic
 * of viscoplastic fluid behavior.  
 */
event logWriting (i++) {
  double ke = 0.;
  foreach (reduction(+:ke)){
    ke += (2*pi*y)*(0.5*rho(f[])*(sq(u.x[]) + sq(u.y[])))*sq(Delta);
  }
  if (pid() == 0){
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
  // Ensure that the cut-off Kinetic energy is smaller than or equal to 1e-6 times the maximum kinetic energy of the system.
  assert(ke > -1e-10);
  assert(ke < 1e3);

  if (ke < 1e-6){
    if (i > 1e2){
      fprintf(ferr, "Kinetic energy is too small. Exiting...\n");
      return 1;
    }
  }

}

/**
## Running the code
~~~bash
#!/bin/bash
qcc -fopenmp -Wall -O2 bounce_VP.c -o bounce_VP -lm -disable-dimensions
export OMP_NUM_THREADS=8
./bounce_VP 10 0.25 1e-2 5.0
~~~
**/