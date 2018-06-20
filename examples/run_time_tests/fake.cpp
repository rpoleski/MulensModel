// g++ -lm fake.cpp ../../source/VBBL/VBBinaryLensingLibrary.o -o fake
// ./fake > fake.out
#include <stdio.h>
#include <math.h>
#include "../../source/VBBL/VBBinaryLensingLibrary.h"

int main() {
  VBBinaryLensing VBBL;
  char coordinatefile[256] = "fake.txt";
  char sattabledir[256] = ".";
  const int np = 1001;
  double t_array[np], y1_array[np], y2_array[np];
  double mag_par_array[np];
  double pr[5];

  VBBL.parallaxsystem = 1;
  VBBL.SetObjectCoordinates(coordinatefile, sattabledir);
  pr[0] = 0.1;
  pr[1] = log(150.);
  pr[2] = 6900.; // It seems VBBL uses short JD.
  pr[3] = 1.;
  pr[4] = 2.;
  for (int i = 0; i < np; i++) 
    t_array[i] = 6820. + i * 160. / np;

  VBBL.PSPLLightCurveParallax(pr, t_array, mag_par_array, y1_array, y2_array, np);
// PSPLLightCurveParallax(double *parameters, double *t_array, double *mag_array, double *y1_array, double *y2_array, int np);
// Parameters are {u0, log_tE, t0, pai1, pai2}

  for (int i = 0; i < np; i++) 
    printf("%lf %lf\n", t_array[i], mag_par_array[i]); 
}

