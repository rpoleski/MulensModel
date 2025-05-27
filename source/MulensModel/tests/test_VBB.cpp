#include <stdio.h>
#include <math.h>
#include "VBBinaryLensingLibrary.h"

int main()
{
	VBBinaryLensing VBBL;

	double pr[15]; // Array of parameters
        double u0, t0, tE, alpha;
        double s, q, Rs;

        s = 1.2; //separation between the two lenses^M
        q = 0.123;
        u0 = 0.555; // Impact parameter^M
        alpha = 0.3054326190991; //17.5; // Angle between a vector pointing to the left and the source velocity^M
        t0 = 5500; // Time of closest approach to the center of mass^M
        tE = 100; // Einstein time^M
        Rs = 0.01; // Source radius in Einstein radii of the total mass.^M

	pr[0] = log(s); // Note that log_s is used as an input parameter to BinaryLightCurve.^M
        pr[1] = log(q);
        pr[2] = u0;
        pr[3] = alpha;
        pr[4] = log(Rs);
        pr[5] = log(tE);
        pr[6] = t0;

        const int np = 101;
        double t_array[np];

        for (int i = 0; i < np; i++) {
                t_array[i] = t0 - 3 * tE + i * (6 * tE / (np - 1));
        }

	//char coordinatefile[256] = "OB151212coords.txt";
        //VBBL.SetObjectCoordinates(coordinatefile);

	double w1, w2, w3;
        w1 = 0.013;
        w2 = -0.2 ;
        w3 = 0.05;

        pr[9] = w1;
        pr[10] = w2;
        pr[11] = w3;

	double mag_orbital_array[np], y1_orbital_array[np], y2_orbital_array[np], sep_array[np];
        VBBL.BinaryLightCurveOrbital(pr, t_array, mag_orbital_array, y1_orbital_array, y2_orbital_array, sep_array, np);

        printf("Trajectory with circular orbital motion\n");
        for (int i = 0; i < np; i++) {
                printf("t: %lf  y1: %lf y2: %lf\n", t_array[i], y1_orbital_array[i], y2_orbital_array[i]);
        }






