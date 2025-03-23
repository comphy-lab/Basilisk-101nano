/* Title: Getting y data from Basilisk file along the x = 0 line
Compile this code using:
qcc -O2 -Wall getDatay.c -o getDatay -lm
# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
*/

#include "utils.h"
#include "output.h"

vector u[];
char filename[80];
scalar shear[], D2p[];


int main(int a, char const *arguments[])
{
  // Boundary condition: periodic right - left
  periodic (right);
  // Slip at the top
  u.t[top] = neumann(0);
  u.n[top] = neumann(0);
  // No slip at the bottom
  u.n[bottom] = dirichlet(0);
  u.t[bottom] = dirichlet(0);
  
  sprintf (filename, "%s", arguments[1]);
  
  restore (file = filename);
  boundary((scalar *){u});

  FILE * fp = ferr;

  foreach(){
    shear[] = (u.x[0,1] - u.x[0,-1])/(2*Delta);
  }

  foreach() {
    // Calculate deformation tensor components at cell centers
    double D11 = (u.x[1,0] - u.x[-1,0])/2.0;
    double D22 = (u.y[0,1] - u.y[0,-1])/2.0;
    double D12 = 0.5*(u.x[0,1] - u.x[0,-1] + u.y[1,0] - u.y[-1,0])/2.0;
    
    // Calculate second invariant
    double D2 = sqrt(sq(D11) + sq(D22) + 2.0 * sq(D12)) / Delta;
    D2p[] = D2;
  }

  for (double y = -0.501; y < 0.499; y += 1e-2){
    fprintf(ferr, "%g %g %g %g\n", y, interpolate(u.x, 0.0, y),
            interpolate(shear, 0.0, y), interpolate(D2p, 0.0, y));
  }
  fflush (fp);
  fclose (fp);
}
