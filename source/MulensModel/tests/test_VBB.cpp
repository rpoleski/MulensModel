#include <stdio.h>
#include <math.h>
#include "VBBinaryLensingLibrary.h"
#include <math.h>

int main()
{
	VBBinaryLensing VBBL;

	double pr[15]; // Array of parameters
        double u0, t0, tE, alpha;
        double s, q, Rs;

        s = 1.2; //separation between the two lenses^M
        q = 0.123; //mass ratio of the lens components
        u0 = 0.555; // Impact parameter^M
        alpha = 17.5 * M_PI/180; // Angle between a vector pointing to the left and the source velocity [radians]
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

        const int np = 5;
        double t_array[np];
        // making the array of times
        for (int i = 0; i < np; i++) {
                t_array[i] = t0 - 3 * tE + i * (6 * tE / (np - 1));
        }

	
        VBBL.SetObjectCoordinates("18:00:00 -26:00:00");

	double ds_dt, dalpha_dt, ds_z_dt;
	double s_z, s_prim, conv_factor;
	double ds_dt_VBB, da_dt_VBB, dsz_dt_VBB;

	ds_dt = 20.2; //the rate of change of the separation per year
        dalpha_dt = -30.3; //rate of change of alpha deg per year
        ds_z_dt = 20; // rate of change of the distance along the line of sight per year

	conv_factor = -ds_dt * s;
        s_z = conv_factor/ds_z_dt; //distance along the line of sight
        s_prim = sqrt(pow(s,2) + pow(s_z,2)); 

        ds_dt_VBB = ds_dt * 1/s_prim * 1/365.2422; // the rate of change of the separation per day
        da_dt_VBB = dalpha_dt * M_PI/180 * 1/365.2422; // rate of change of alpha radian per day
        dsz_dt_VBB = ds_z_dt * 1/s_prim * 1/365.2422; // rate of change of the distance along the line of sight per year


        pr[9] = ds_dt_VBB;
        pr[10] = da_dt_VBB;
        pr[11] = dsz_dt_VBB;

	double mag_orbital_array[np], y1_orbital_array[np], y2_orbital_array[np], sep_array[np];
	// calculating the trajectory 
        VBBL.BinaryLightCurveOrbital(pr, t_array, mag_orbital_array, y1_orbital_array, y2_orbital_array, sep_array, np);
        

	printf("[%lf, %lf, %lf, %lf, %lf]\n", y1_orbital_array[0], y1_orbital_array[1], y1_orbital_array[2], y1_orbital_array[3], y1_orbital_array[4]);
	printf("[%lf, %lf, %lf, %lf, %lf]\n", y2_orbital_array[0], y2_orbital_array[1], y2_orbital_array[2], y2_orbital_array[3], y2_orbital_array[4]);
	printf("[%lf, %lf, %lf, %lf, %lf]\n", sep_array[0], sep_array[1], sep_array[2], sep_array[3], sep_array[4]);


}





