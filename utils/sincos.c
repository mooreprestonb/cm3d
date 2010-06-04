
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define NMAX 1000

int main(){
  int i;
  double angle;
  double dx;

  dx = M_PI/(NMAX-1.0);
  for (i=0;i<NMAX;++i){
    angle = -M_PI/2 + i*dx;
    printf("%g %g\n",angle,1./sin(angle));
  }
  exit(0);
}
