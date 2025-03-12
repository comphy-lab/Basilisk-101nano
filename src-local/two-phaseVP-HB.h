/**
 * # Two-Phase Flow with Viscoplastic (Herschel-Bulkley) Fluid Model
 * 
 * Author: Vatsal Sanjay 
 * Version 4.0, Jan 1, 2025
 * 
 * ## Purpose
 * This header extends Basilisk's two-phase flow capabilities to include viscoplastic 
 * fluids following the Herschel-Bulkley rheological model. It enables simulations of 
 * complex fluids with yield stress (like mud, toothpaste, some polymers) interacting 
 * with Newtonian fluids (like air or water).
 * 
 * ## Usage
 * Include this header after your Navier-Stokes solver and before other modules:
 * ```c
 * #include "navier-stokes/centered.h"
 * #include "two-phaseVP-HB.h"  // This file
 * ```
 * 
 * ## Key Parameters
 * - `rho1`, `mu1`: Density and viscosity of fluid 1 (can be viscoplastic)
 * - `rho2`, `mu2`: Density and viscosity of fluid 2 (always Newtonian)
 * - `tauy`: Yield stress parameter for fluid 1
 * - `n`: Power-law index for Herschel-Bulkley model
 * - `epsilon`: Regularization parameter for numerical stability
 * 
 * ## Changelog
 * 
 * ### Jan 1, 2025 (v4.0)
 * - Fixed D2 calculation bug that occurred with parallel computing
 * 
 * ### Dec 31, 2024 (v3.5)
 * - Added full Herschel-Bulkley formulation:
 *   - n = 1: Recovers Bingham model
 *   - n = 1, tau_y = 0: Recovers Newtonian fluid
 *   - n < 1: Shear-thinning behavior
 *   - n > 1: Shear-thickening behavior
 * 
 * ### Oct 11, 2024 (v2.0) 
 * - Implemented epsilon regularization for viscoplastic fluids
 * - Unified planar and axisymmetric handling in single codebase
 * - Made only fluid 1 eligible for viscoplastic properties
 * 
 * ### Development History 
 * - v0.0: Original Basilisk two-phase.h for Newtonian fluids
 * - v1.0: Added min(temp, mumax) viscoplastic formulation with separate axisymmetric code
 * - v2.0: Implemented epsilon regularization with unified code for all geometries
 * - v3.5: Extended to Herschel-Bulkley model
 * - v4.0: Bug fixes and optimizations
 *
 * ## Physical Background
 * This implementation models two-phase interfacial flows using the Volume-Of-Fluid method.
 * The volume fraction f=1 in fluid 1 (potentially viscoplastic) and f=0 in fluid 2 (Newtonian).
 * Interface tracking is handled by Basilisk's standard VoF approach.
 */

#include "vof.h"

scalar f[], * interfaces = {f};
scalar D2[];
face vector D2f[];
double rho1 = 1., mu1 = 0., rho2 = 1., mu2 = 0.;
double epsilon = 1e-6, tauy = 0., n = 1.;
/**
Auxilliary fields are necessary to define the (variable) specific
volume $\alpha=1/\rho$ as well as the cell-centered density. */

face vector alphav[];
scalar rhov[];

event defaults (i = 0) {
  alpha = alphav;
  rho = rhov;

  /**
  If the viscosity is non-zero, we need to allocate the face-centered
  viscosity field. */

  mu = new face vector;
}

/**
The density and viscosity are defined using arithmetic averages by
default. The user can overload these definitions to use other types of
averages (i.e. harmonic). */

#ifndef rho
# define rho(f) (clamp(f,0.,1.)*(rho1 - rho2) + rho2)
#endif
#ifndef mu
// for Arithmetic mean, use this
# define mu(muTemp, mu2, f)  (clamp(f,0.,1.)*(muTemp - mu2) + mu2)
// for Harmonic mean, use this
// # define mu(muTemp, mu2, f) (1.0 / ((clamp(f,0.,1.) / muTemp) + ((1.0 - clamp(f,0.,1.)) / mu2)))
#endif

/**
We have the option of using some "smearing" of the density/viscosity
jump. */

#ifdef FILTERED
scalar sf[];
#else
# define sf f
#endif

event tracer_advection (i++) {

  /**
  When using smearing of the density jump, we initialise *sf* with the
  vertex-average of *f*. */

#ifndef sf
#if dimension <= 2
  foreach()
    sf[] = (4.*f[] +
	    2.*(f[0,1] + f[0,-1] + f[1,0] + f[-1,0]) +
	    f[-1,-1] + f[1,-1] + f[1,1] + f[-1,1])/16.;
#else // dimension == 3
  foreach()
    sf[] = (8.*f[] +
	    4.*(f[-1] + f[1] + f[0,1] + f[0,-1] + f[0,0,1] + f[0,0,-1]) +
	    2.*(f[-1,1] + f[-1,0,1] + f[-1,0,-1] + f[-1,-1] +
		f[0,1,1] + f[0,1,-1] + f[0,-1,1] + f[0,-1,-1] +
		f[1,1] + f[1,0,1] + f[1,-1] + f[1,0,-1]) +
	    f[1,-1,1] + f[-1,1,1] + f[-1,1,-1] + f[1,1,1] +
	    f[1,1,-1] + f[-1,-1,-1] + f[1,-1,-1] + f[-1,-1,1])/64.;
#endif
#endif

#if TREE
  sf.prolongation = refine_bilinear;
  sf.dirty = true; // boundary conditions need to be updated
#endif
}

event properties (i++) {
/**
 * ## Viscoplastic Fluid Model Implementation
 * 
 * Here we implement the Herschel-Bulkley rheological model, which represents a significant
 * extension to the standard Newtonian fluid model. This is where the core modifications
 * from the original two-phase.h file are located.
 * 
 * ### Herschel-Bulkley Model: Theory and Components
 * 
 * The Herschel-Bulkley model combines three key rheological behaviors:
 * 
 * 1. **Yield Stress (τ_y)**: The material behaves as a rigid body until stress exceeds threshold
 * 2. **Power-Law Flow (n)**: Controls shear-thinning (n<1) or shear-thickening (n>1) behavior
 * 3. **Consistency (K)**: Relates to the fluid's effective viscosity
 * 
 * The constitutive equation for the stress tensor is:
 * 
 * $$
 * \boldsymbol{\tau} = 
 * \tau_{y}\,\boldsymbol{\mathcal{I}} \;+\; K\left(2\boldsymbol{\mathcal{D}}\right)^{n}
 * $$
 * 
 * ### Numerical Implementation: ε-Regularization
 * 
 * For numerical stability, we use an ε-regularization approach that avoids the mathematical 
 * singularity at zero strain rate. This transforms the equation to:
 * 
 * $$
 * \boldsymbol{\tau} =
 * 2\biggl[\frac{\tau_{y}}{2\|\boldsymbol{\mathcal{D}}\|+\varepsilon}\,\boldsymbol{\mathcal{I}}
 * +
 * K\,\bigl(2\|\boldsymbol{\mathcal{D}}\|+\epsilon\bigr)^{n-1}
 * \biggr]\boldsymbol{\mathcal{D}}
 * $$
 * 
 * After non-dimensionalizing with capillary scaling (stresses with γ/R₀, length with R₀, 
 * and velocity with √(γ/ρₗR₀)), we get:
 * 
 * $$
 * \boldsymbol{\tilde{\tau}} =
 * 2\biggl[\frac{\mathcal{J}}{2\|\boldsymbol{\tilde{\mathcal{D}}}\|+\varepsilon}\,\boldsymbol{\mathcal{I}}
 * +
 * Oh_K\,\bigl(2\|\boldsymbol{\tilde{\mathcal{D}}}\|+\epsilon\bigr)^{n-1}
 * \biggr]\boldsymbol{\tilde{\mathcal{D}}}
 * $$
 * 
 * where J is the dimensionless yield stress (capillary-Bingham number) and Oh_K is an effective Ohnesorge number:
 * 
 * $$
 * Oh_K = \frac{K}{\sqrt{\rho_l^n\gamma^{2-n}R_0^{3n-2}}}
 * $$
 * 
 * ### Special Cases
 * 
 * 1. When n = 1: Reduces to Bingham model with Oh = μₗ/√(ρₗγR₀)
 * 2. When n = 1 and J = 0: Reduces to standard Newtonian fluid
 * 3. When n < 1: Models shear-thinning fluids (viscosity decreases with strain rate)
 * 4. When n > 1: Models shear-thickening fluids (viscosity increases with strain rate)
 * 
 * ### Strain Rate Tensor Calculation
 * 
 * For axisymmetric coordinates (r,z), the deformation tensor components are:
 * 
 * $$\mathcal{D}_{11} = \frac{\partial u_r}{\partial r}$$
 * $$\mathcal{D}_{22} = \frac{u_r}{r}$$
 * $$\mathcal{D}_{13} = \mathcal{D}_{31} = \frac{1}{2}\left(\frac{\partial u_r}{\partial z}+ \frac{\partial u_z}{\partial r}\right)$$
 * $$\mathcal{D}_{33} = \frac{\partial u_z}{\partial z}$$
 * $$\mathcal{D}_{12} = \mathcal{D}_{23} = 0$$
 * 
 * The norm is calculated using the Frobenius norm: $\mathcal{D}_2=\sqrt{\mathcal{D}_{ij}\mathcal{D}_{ij}}$
 * Where $\|\mathcal{D}\| = D_2/\sqrt{2}$ is the magnitude used in the constitutive equation.
 * 
 * ### Important Implementation Notes
 * 
 * 1. This model always produces flow (never true rigid behavior) but approximates yield
 *    with very high viscosity below the yield stress.
 * 
 * 2. The ε parameter controls regularization - smaller values better approximate true yield 
 *    behavior but may create numerical challenges.
 * 
 * 3. The implementation follows Balmforth et al. (2013) approach where 
 *    $\dot{\mathcal{S}}_{ij}=2 D_{ij}$ relates strain rate to deformation tensors.
 * 
 * References:
 * - Balmforth et al. (2013). "Yielding to Stress: Recent Developments in Viscoplastic Fluid Mechanics"
 * - Sanjay et al. (2021). "Bursting Bubble in a Viscoplastic Medium"
 */

  foreach_face(x) {
    double ff = (sf[] + sf[-1])/2.;
    alphav.x[] = fm.x[]/rho(ff);
    double muTemp = mu1;
    
    face vector muv = mu;

    double D2temp = 0.;

    D2temp += sq(0.5*( (u.y[0,1] - u.y[0,-1] + u.y[-1,1] - u.y[-1,-1])/(2.*Delta) )); // D11
#if AXI
    D2temp += sq((u.y[0,0] + u.y[-1, 0])/(2*max(y, 1e-20))); // D22
#endif
    D2temp += sq((u.x[] - u.x[-1,0])/Delta); // D33
    D2temp += 2.0*sq(0.5*( (u.y[] - u.y[-1, 0])/Delta + 0.5*( (u.x[0,1] - u.x[0,-1] + u.x[-1,1] - u.x[-1,-1])/(2.*Delta) ) )); // D13

    D2temp = sqrt(D2temp/2.0);

    if (tauy > 0.){
      muTemp = tauy/(2.0*D2temp + epsilon) + mu1*pow((2.0*D2temp + epsilon), n-1);
    }
    
    muv.x[] = fm.x[]*mu(muTemp, mu2, ff);
    D2f.x[] = D2temp;
  }

  foreach_face(y) {
    double ff = (sf[0,0] + sf[0,-1])/2.;
    alphav.y[] = fm.y[]/rho(ff);
    double muTemp = mu1;
    face vector muv = mu;

    double D2temp = 0.;

    D2temp += sq((u.y[0,0] - u.y[0,-1])/Delta); // D11
#if AXI
    D2temp += sq((u.y[0,0] + u.y[0,-1])/(2*max(y, 1e-20))); // D22
#endif
    D2temp += sq(0.5*( (u.x[1,0] - u.x[-1,0] + u.x[1,-1] - u.x[-1,-1])/(2.*Delta) )); // D33
    D2temp += 2.0*sq(0.5*( (u.x[0,0] - u.x[0,-1])/Delta + 0.5*( (u.y[1,0] - u.y[-1,0] + u.y[1,-1] - u.y[-1,-1])/(2.*Delta) ) )); // D13

    D2temp = sqrt(D2temp/2.0);

    if (tauy > 0.){
      muTemp = tauy/(2.0*D2temp + epsilon) + mu1*pow((2.0*D2temp + epsilon), n-1);
    }

    muv.y[] = fm.y[]*mu(muTemp, mu2, ff);
    D2f.y[] = D2temp;
  }

#if dimension == 3
  error("3D not implemented yet");
#endif

  /**
   * ## Cell-Centered Deformation Rate Storage
   * 
   * We also calculate a cell-centered scalar field D2 that stores the magnitude of the
   * deformation rate tensor $\|\mathbf{\mathcal{D}}\|$ throughout the domain. This serves
   * two important purposes:
   * 
   * 1. **Yield Surface Detection**: Allows identification of "fake-yield surfaces" where 
   *    the strain rate is near the critical threshold (|D| ≈ τ_y/μ)
   * 
   * 2. **Adaptive Mesh Refinement**: Enables more accurate mesh refinement in regions 
   *    where the material transitions between yielded and unyielded states
   *    
   * This field can be used for both visualization and as a criterion for adaptive refinement.
   */
  foreach(){
    rhov[] = cm[]*rho(sf[]);
    D2[] = f[]*(D2f.x[]+D2f.y[]+D2f.x[1,0]+D2f.y[0,1])/4.;
  }
#if TREE
  sf.prolongation = fraction_refine;
  sf.dirty = true; // boundary conditions need to be updated
#endif
}
