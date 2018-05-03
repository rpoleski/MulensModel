#include <iostream>

#include "VBBinaryLensingLibrary.h"

extern "C" {
  double VBBinaryLensing_BinaryMagDark(double a,double q,double y1,double y2,double RSv,double a1, double Tol) {
    static VBBinaryLensing VBBL;
    
    return VBBL.BinaryMagDark(a, q, y1, y2, RSv, a1, Tol);
  }
}

extern "C" double* VBBL_SG12_5(double p0, double p1, 
                    double p2, double p3, double p4, double p5, double p6, 
                    double p7, double p8, double p9, double p10, double p11) {
    static VBBinaryLensing VBBL;
    complex complex_poly[6], complex_roots[5];
    static double roots[10];
    int i;

    complex_poly[0] = complex(p0, p6);
    complex_poly[1] = complex(p1, p7);
    complex_poly[2] = complex(p2, p8);
    complex_poly[3] = complex(p3, p9);
    complex_poly[4] = complex(p4, p10);
    complex_poly[5] = complex(p5, p11);

    VBBL.cmplx_roots_gen(complex_roots, complex_poly, 5, true, true);
    
    for (i=0; i<5; i++) {
        roots[i] = complex_roots[i].re;
        roots[i+5] = complex_roots[i].im;
    }

    return roots;
}
