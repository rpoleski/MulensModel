// VBBinaryLensing v1.2.1
// This code has been developed by Valerio Bozza, University of Salerno.
// Any use of this code for scientific publications should be acknowledged by a citation to
// V. Bozza, MNRAS 408 (2010) 2188

//#pragma once

#ifndef __binlens
#define __binlens
#define __unmanaged

#define _L1 x1-((x1+a/2.0)/((x1+a/2.0)*(x1+a/2.0)+x2*x2)+q*(x1-a/2.0)/((x1-a/2.0)*(x1-a/2.0)+x2*x2))/(1.0+q) // Used in PlotCrits
#define _L2 x2-(x2/((x1+a/2.0)*(x1+a/2.0)+x2*x2)+q*x2/((x1-a/2.0)*(x1-a/2.0)+x2*x2))/(1.0+q)
#define _LL (y-z)+coefs[21]/(zc-coefs[20])+coefs[22]/zc //Lens equation test
#define _J1 m1/((zc-m2*a)*(zc-m2*a))+m2/((zc+m1*a)*(zc+m1*a)) //#define _J1 m1/((zc-0.5*a)*(zc-0.5*a))+m2/((zc+0.5*a)*(zc+0.5*a))
#define _J1c coefs[21]/((zc-coefs[20])*(zc-coefs[20]))+coefs[22]/(zc*zc) //#define _J1 m1/((zc-0.5*a)*(zc-0.5*a))+m2/((zc+0.5*a)*(zc+0.5*a))
#define _J2 -2.0*(coefs[21]/((z-coefs[20])*(z-coefs[20])*(z-coefs[20]))+coefs[22]/(z*z*z))
#define _skew(p1,p2,q1,q2) p1*q2-p2*q1
#define _NP 200.0

#define _sign(x) ((x>0)? +1 : -1)

//using namespace System;

class _curve;
class _sols;
class _theta;
class complex;

#ifndef __unmanaged
namespace VBBinaryLensingLibrary {

	public ref class VBBinaryLensing
#else
	class VBBinaryLensing
#endif
	{
		int *ndatasat;
		double **tsat,***possat;
		double e,phi,phip,phi0,Om,inc,t0,d3,v3,GM,flagits;
		double Obj[3],rad[3],tang[3];
		double Eq2000[3],Quad2000[3],North2000[3]; 

		_curve *NewImages(complex,complex  *,_theta *);
		void OrderImages(_sols *,_curve *);
		void ComputeParallax(double,double,double *);
		void laguer(complex *, int, complex *, int*,double);
		void zroots(complex *,int, complex *, int,double);

	public: 

		double Tol;
		int satellite,parallaxsystem,nsat;
		int minannuli,nannuli;
		double y_1,y_2,av;

		_sols *PlotCrit(double,double);
		void PrintCau(double,double);

		double BinaryMag0(double,double,double,double, _sols **);
		double BinaryMag0(double,double,double,double);
		double BinaryMag(double,double,double,double,double,double, _sols **);
		double BinaryMag(double,double,double,double,double,double);
		double BinaryMagDark(double,double,double,double,double,double,double);

		void SetObjectCoordinates(char *,char *);

		double PSPLCurve(double *,double);
		double PSPLParallaxCurve(double *,double);
		double ESPLCurve(double *,double);
		double ESPLParallaxCurve(double *,double);

		double BinaryLightCurve(double *,double);
		double BinaryLightCurveW(double *,double);
		double BinaryLightCurveParallax(double *,double);
		double BinaryLightCurveOrbital(double *,double);
		double BinSourceMag(double *,double);
		double BinSourceParallaxMag(double *,double);
		double BinSourceXallarapMag(double *,double);

		VBBinaryLensing();
		~VBBinaryLensing();

	};


	struct annulus{
		double bin;
		double cum;
		double Mag;
		double err;
		double f;
		int nim;
		annulus *prev,*next;
	};


#ifndef __unmanaged
}
#endif


class _theta{
public: 
	double th,maxerr,Mag,errworst;
	_theta *prev,*next;

	_theta(double);

};

class _thetas{
public: 
	_theta *first,*last;
	int length;

	_thetas(void);
	~_thetas(void);
	_theta *insert(double);


};

class complex{
public:
	double re;
	double im;
	complex(double,double);
	complex(double);
	complex(void);
};

class _point{
public:
	double x1;
	double x2;
	double parab,ds,dJ;
	complex d;
	_theta *theta;
	_point(double ,double,_theta *);
	_point *next,*prev;
	double operator-(_point);
};

class _curve{
public:
	int length;
	_point *first,*last;
	_curve *next,*prev;
	_curve *partneratstart,*partneratend;
	double parabstart;

	_curve(_point *);
	_curve(void);
	~_curve(void);

	_curve *divide(_point *);
	void drop(_point *);
	void append(double,double);
	void append(_point *);
	void prepend(double,double);
//	void prepend(_point *);
	_curve *join(_curve *);
	_curve *joinbefore(_curve *);
	_curve *reverse(void);
	double closest(_point *,_point **);
	double closest2(_point *,_point **);
	void complement(_point **,int,_point **,int);
};

class _sols{
public:
	int length;
	_curve *first,*last;

	_sols(void);
	~_sols(void);
	void drop(_curve *);
	void append(_curve *);
	void prepend(_curve *);
	void join(_sols *);	
};


#define MR 8
#define MT 10
#define MAXIT (MT*MR)
#define MAXM 30

//extern double y_1,y_2,e,a,av,phi,phip,phi0,Om,inc,t0,d3,v3,GM,flagits;



double rf(double x, double y, double z);
double rd(double x, double y, double z);
double rj(double x, double y, double z, double p);
double rc(double x, double y);
double ellf(double phi, double ak);
double elle(double phi, double ak);
double ellpi(double phi, double en, double ak);


#endif
