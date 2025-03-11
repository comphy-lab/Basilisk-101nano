/**
# Three-phase interfacial flows

This file provides functionality for simulating flows involving three immiscible fluid phases,
each separated by interfaces. It extends Basilisk's two-phase capabilities to handle
three distinct phases, which is essential for many complex multiphase flow problems
such as compound droplets, liquid lenses, and multiple emulsions.

## Conceptual Overview

The method tracks three phases using two Volume-Of-Fluid (VOF) fields:
- `f1`: First volume fraction field (f1=1 in phase 1, f1=0 elsewhere)
- `f2`: Second volume fraction field (f2=1 in phase 2, but only meaningful where f1=1)

The three phases are defined as:
1. Phase 1: regions where f1=1 and f2=0
2. Phase 2: regions where f1=1 and f2=1 
3. Phase 3: regions where f1=0 (regardless of f2)

## Physical Properties

Each phase has its own physical properties:
- Densities: *rho1*, *rho2*, *rho3*
- Dynamic viscosities: *mu1*, *mu2*, *mu3*

These properties are combined using appropriate averaging at interfaces to calculate
the effective properties needed by the Navier-Stokes solver.

## Surface Tension Relationships

In three-phase systems, the surface tensions between phase pairs must follow physical
constraints. When three interfaces meet at a contact line, they form specific angles
dictated by the balance of surface tensions.

For an ideal system with a precursor film, the surface tensions follow:
σ₂₃ = σ₁₃ + σ₁₂

Where:
- σ₁₂: Surface tension between phases 1 and 2
- σ₁₃: Surface tension between phases 1 and 3
- σ₂₃: Surface tension between phases 2 and 3

This relationship (known as Neumann's triangle) determines the equilibrium configuration
of three-phase systems such as liquid lenses at interfaces.

## Usage

This header is typically used in combination with a Navier-Stokes solver like 
`navier-stokes/centered.h`. To use it:

1. Include this file after the Navier-Stokes solver
2. Set the density and viscosity values for all three phases
3. Initialize the volume fraction fields f1 and f2
4. For surface tension, set separate tension coefficients for each interface

## Implementation Details

The file implements:
1. Interface tracking through VOF fields
2. Property averaging at interfaces (density and viscosity)
3. Optional filtering/smearing of properties to improve numerical stability
4. Proper calculation of face-centered and cell-centered properties

This implementation is fully compatible with Basilisk's adaptive mesh refinement.
*/

#include "vof.h"
scalar f1[], f2[], *interfaces = {f1, f2};      // Define the two VOF fields and group them
double rho1 = 1., mu1 = 0., rho2 = 1., mu2 = 0., rho3 = 1., mu3 = 0.;  // Default properties

/**
## Auxiliary Fields

These fields are necessary for the Navier-Stokes solver:
- `alphav`: Face-centered specific volume (1/ρ) field
- `rhov`: Cell-centered density field

These allow the solver to handle variable density and viscosity across interfaces.
*/

face vector alphav[];    // Face-centered specific volume
scalar rhov[];           // Cell-centered density

/**
## Default Event

This event sets up the necessary fields during initialization.
*/
event defaults (i = 0) {
  alpha = alphav;        // Connect Navier-Stokes alpha to our specific volume field
  rho = rhov;            // Connect Navier-Stokes rho to our density field

  /**
  Create a face-centered viscosity field if needed (when viscosity is non-zero).
  This is used by the Navier-Stokes solver for the viscous term in the momentum equation.
  */
  mu = new face vector;  // Allocate memory for face-centered viscosity
}

/**
## Property Calculations

These macros define how physical properties are calculated in mixed cells 
(cells containing interfaces). The default implementation uses volume-weighted
arithmetic averages, but users can define their own averaging methods.

### How the mixing rules work:
- `clamp(f1*(1-f2), 0., 1.)` gives the volume fraction of phase 1
- `clamp(f1*f2, 0., 1.)` gives the volume fraction of phase 2
- `clamp((1-f1), 0., 1.)` gives the volume fraction of phase 3

The clamp function ensures values stay within [0,1] for numerical stability.
*/

#ifndef rho
#define rho(f1, f2) (clamp(f1*(1-f2), 0., 1.) * rho1 + clamp(f1*f2, 0., 1.) * rho2 + clamp((1-f1), 0., 1.) * rho3)
#endif
#ifndef mu
#define mu(f1, f2) (clamp(f1*(1-f2), 0., 1.) * mu1 + clamp(f1*f2, 0., 1.) * mu2 + clamp((1-f1), 0., 1.) * mu3)
#endif

/**
## Filtering Option

For improved numerical stability, we can apply filtering (smearing) to the 
volume fraction fields before calculating properties. This is activated
by defining the `FILTERED` macro before including this header.

When filtering is enabled, smoothed fields (sf1, sf2) are used for property 
calculations while the original fields (f1, f2) are used for interface tracking.
*/

#ifdef FILTERED
scalar sf1[], sf2[], *smearInterfaces = {sf1, sf2};   // Smoothed fields for property calculations
#else
#define sf1 f1                                        // Without filtering, use original fields
#define sf2 f2
scalar *smearInterfaces = {sf1, sf2};
#endif

/**
## Properties Event

This event recalculates physical properties at each timestep based on the 
current interface positions. It is called automatically by Basilisk's
event system.
*/
event properties (i++) {
  /**
  When filtering is enabled, we calculate smoothed volume fractions by 
  weighted averaging of neighboring cells. This reduces numerical instabilities
  caused by sharp property jumps at interfaces.
  
  The filtering stencil uses different weights for:
  - Center cell: highest weight
  - Face neighbors: medium weight
  - Corner neighbors: lowest weight
  
  This creates a smooth transition of properties across interfaces.
  */
#ifdef FILTERED
  int counter1 = 0;
  for (scalar sf in smearInterfaces){
    counter1++;
    int counter2 = 0;
    for (scalar f in interfaces){
      counter2++;
      if (counter1 == counter2){
        // fprintf(ferr, "%s %s\n", sf.name, f.name);
      #if dimension <= 2
          foreach(){
            sf[] = (4.*f[] +                          // Center cell (weight 4)
        	    2.*(f[0,1] + f[0,-1] + f[1,0] + f[-1,0]) +  // Face neighbors (weight 2)
        	    f[-1,-1] + f[1,-1] + f[1,1] + f[-1,1])/16.;  // Corner neighbors (weight 1)
          }
      #else // dimension == 3
          foreach(){
            sf[] = (8.*f[] +                          // Center cell (weight 8)
        	    4.*(f[-1] + f[1] + f[0,1] + f[0,-1] + f[0,0,1] + f[0,0,-1]) +  // Face neighbors (weight 4)
        	    2.*(f[-1,1] + f[-1,0,1] + f[-1,0,-1] + f[-1,-1] +
        		f[0,1,1] + f[0,1,-1] + f[0,-1,1] + f[0,-1,-1] +
        		f[1,1] + f[1,0,1] + f[1,-1] + f[1,0,-1]) +  // Edge neighbors (weight 2)
        	    f[1,-1,1] + f[-1,1,1] + f[-1,1,-1] + f[1,1,1] +
        	    f[1,1,-1] + f[-1,-1,-1] + f[1,-1,-1] + f[-1,-1,1])/64.;  // Corner neighbors (weight 1)
          }
      #endif
      }
    }
  }
  #endif
  
  /**
  For adaptive mesh refinement, we need to set the prolongation method for
  the smoothed fields. This determines how values are interpolated when 
  the mesh is refined.
  
  Here, we temporarily use bilinear interpolation for the filtering step.
  */
#if TREE
  for (scalar sf in smearInterfaces){
    sf.prolongation = refine_bilinear;  // Use bilinear interpolation for filtering
    boundary ({sf});                     // Update boundary values
  }
#endif

  /**
  Calculate face-centered specific volume and viscosity using the
  (possibly filtered) volume fraction fields. These are used by the
  Navier-Stokes solver.
  
  For each face, we:
  1. Calculate average volume fractions at the face
  2. Compute density using the mixing rule
  3. Convert to specific volume (1/density) and scale by face fraction
  4. Calculate viscosity using the mixing rule and scale by face fraction
  */
  foreach_face() {
    double ff1 = (sf1[] + sf1[-1])/2.;  // Average f1 at the face
    double ff2 = (sf2[] + sf2[-1])/2.;  // Average f2 at the face

    alphav.x[] = fm.x[]/rho(ff1, ff2);  // Specific volume scaled by face fraction
    face vector muv = mu;
    muv.x[] = fm.x[]*mu(ff1, ff2);      // Viscosity scaled by face fraction
  }
  
  /**
  Calculate cell-centered density for each cell using the mixing rule.
  This is scaled by the cell fraction (cm) which is important for 
  embedded boundary methods.
  */
  foreach(){
    rhov[] = cm[]*rho(sf1[], sf2[]);    // Density scaled by cell fraction
  }

  /**
  After property calculations, we restore the proper prolongation method
  for the volume fraction fields. For VOF methods, this should be
  the fraction_refine function which preserves volume fractions during
  mesh refinement.
  */
#if TREE
  for (scalar sf in smearInterfaces){
    sf.prolongation = fraction_refine;  // Restore VOF-appropriate prolongation
    boundary ({sf});                     // Update boundary values
  }
#endif
}
