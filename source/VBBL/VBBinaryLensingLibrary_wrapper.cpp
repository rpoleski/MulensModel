#include <iostream>

#include "VBBinaryLensingLibrary.h"

extern "C" {
  double VBBinaryLensing_BinaryMagDark(double a, double q, double y1, double y2, double RSv, double a1, double Tol) {
    static VBBinaryLensing VBBL;
    
    return VBBL.BinaryMagDark(a, q, y1, y2, RSv, a1, Tol);
  }
}

extern "C" {
  double VBBinaryLensing_BinaryMag0(double a,double q,double y1,double y2, double K, double G, double Gi) {
    static VBBinaryLensing VBBL;
    
    return VBBL.BinaryMag0_shear(a, q, y1, y2, K, G, Gi);
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

extern "C" {
  double VBBL_BinaryMag(double a, double q, double y1, double y2) {
    static VBBinaryLensing VBBL;

    return VBBL.BinaryMag0(a, q, y1, y2);
  }
}


extern "C" double* VBBL_SG12_9(double p0, double p1, 
                    double p2, double p3, double p4, double p5, double p6, 
                    double p7, double p8, double p9, double p10, double p11,
                    double p12, double p13, double p14, double p15,
                    double p16, double p17, double p18, double p19) {
    static VBBinaryLensing VBBL;
    complex complex_poly[10], complex_roots[9];
    static double roots[18];
    int i;

    complex_poly[0] = complex(p0, p10);
    complex_poly[1] = complex(p1, p11);
    complex_poly[2] = complex(p2, p12);
    complex_poly[3] = complex(p3, p13);
    complex_poly[4] = complex(p4, p14);
    complex_poly[5] = complex(p5, p15);
    complex_poly[6] = complex(p6, p16);
    complex_poly[7] = complex(p7, p17);
    complex_poly[8] = complex(p8, p18);
    complex_poly[9] = complex(p9, p19);

    VBBL.cmplx_roots_gen(complex_roots, complex_poly, 9, true, true);
    
    for (i=0; i<9; i++) {
        roots[i] = complex_roots[i].re;
        roots[i+9] = complex_roots[i].im;
    }

    return roots;
}
