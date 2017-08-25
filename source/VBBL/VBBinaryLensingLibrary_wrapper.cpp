#include <iostream>

#include "VBBinaryLensingLibrary.h"

extern "C" {
  double VBBinaryLensing_BinaryMagDark(double a,double q,double y1,double y2,double RSv,double a1, double Tol) {
    VBBinaryLensing VBBL;
    
    return VBBL.BinaryMagDark(a, q, y1, y2, RSv, a1, Tol);
  }
}

