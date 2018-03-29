#include <stdlib.h>

#ifndef FLOAT
#define FLOAT double
#endif

FLOAT ld_linear(int n, FLOAT gam[], FLOAT rho);

FLOAT mag_binext(FLOAT y1, FLOAT y2, FLOAT rho, FLOAT d, FLOAT q,
        FLOAT (*ld_func)(int,FLOAT*,FLOAT), int n, FLOAT gam[],
	FLOAT acc, FLOAT ld_acc);

FLOAT mag_binpt(FLOAT y1, FLOAT y2, FLOAT d, FLOAT q);

double Adaptive_Contouring_Linear(double d, double q, double y1, double y2, 
		double rho, double gamma, double acc, double ld_acc) {
// Wrapper around Martin Dominik's Adaptive Contouring code.
  FLOAT gam[1];
  gam[0] = gamma;

  if (gamma == 0.0) 
    return(mag_binext(y1,y2,rho,d,q,NULL,-1,NULL,acc,ld_acc));
  else 
    return(mag_binext(y1,y2,rho,d,q,ld_linear,1,gam,acc,ld_acc));  
}

double Adaptive_Contouring_Point_Source(double d, double q, double y1, double y2) {
  return mag_binpt(y1, y2, d, q);
}

