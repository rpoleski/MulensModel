/* binext_adap.c */
/* Calculate extended-source magnification for binary lens
     with adaptive contouring algorithm */
/* MD 29-5-06 */
/* MD 14-7-06 */
/* MD 14-9-06 */
/* RP 19-12-17 */
/* RP 02-02-18 */
/* RP 16-03-18 */
/* RP 15-11-18 */

#define _GNU_SOURCE
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#ifndef FLOAT
#define FLOAT double
#endif

#ifndef boolean
typedef enum {FALSE,TRUE} boolean;
#endif

#define SQ(a) ((a)*(a))

#define BIGG 1e37

typedef struct ptlist {
	FLOAT x1,x2;
	int min_lvl;
	struct ptlist *next;
} PTLIST, *PTLIST_PTR;

	
void sol_binpt(FLOAT y1, FLOAT y2, FLOAT d, FLOAT q, PTLIST_PTR *ptlist);
FLOAT mag_binpt(FLOAT y1, FLOAT y2, FLOAT d, FLOAT q);
boolean inside_circ(FLOAT y1, FLOAT y2, FLOAT *q);
FLOAT rho_circ(FLOAT y1, FLOAT y2, FLOAT *q);
FLOAT ld_linear(int n, FLOAT gam[], FLOAT rho);
void lenseq_bin(FLOAT x1, FLOAT x2, FLOAT *y1, FLOAT *y2, 
	 FLOAT *p);
FLOAT adaptive_contour(FLOAT acc, FLOAT ld_acc, 
        PTLIST_PTR *pointlist, PTLIST_PTR *holelist,
	FLOAT (*ld_func)(int n, FLOAT gam[], FLOAT rho),
	int n, FLOAT gam[], 
	FLOAT (*rho_func)(FLOAT, FLOAT, FLOAT*),
        boolean (*inside_func)(FLOAT, FLOAT, FLOAT*),
        void (*lenseq_func)(FLOAT, FLOAT, FLOAT*, FLOAT*, FLOAT*),
	int ns, int nl, ...);
FLOAT mag_binext(FLOAT y1, FLOAT y2, FLOAT rho, FLOAT d, FLOAT q, 
	FLOAT (*ld_func)(int,FLOAT*,FLOAT), int n, FLOAT gam[],
	FLOAT acc, FLOAT ld_acc);
FLOAT magt_binext(FLOAT t, FLOAT t0, FLOAT tE, FLOAT alpha, FLOAT u0,
	FLOAT rho, FLOAT d, FLOAT q, 
	FLOAT (*ld_func)(int,FLOAT*,FLOAT), int n, FLOAT gam[],
	FLOAT acc, FLOAT ld_acc);

// from erdlcaust.c :
void get_crit_circ(FLOAT d, FLOAT m1, FLOAT y1c, FLOAT y2c, FLOAT rcirc,
	PTLIST_PTR *crit_points);

void lenseq_bin(FLOAT x1, FLOAT x2, FLOAT *y1, FLOAT *y2, 
	 FLOAT *p)
{
	FLOAT d,q;
	FLOAT m1,m2;
	FLOAT x11,x12,d1,d2;
	FLOAT x22, m1d1, m2d2;

	d = p[0];
	q = p[1];

	m1 = 1.0/(1.0+q);
	m2 = q * m1;
	x11 = x1-m2*d;
	x12 = x1+m1*d;
	x22 = x2 * x2;
	d1 = SQ(x11)+x22;
	d2 = SQ(x12)+x22;
	if (d1 == 0.0 || d2 == 0.0) {
          *y1 = -BIGG;
	  *y2 = -BIGG; 
	} else {
	  m1d1 = m1/d1;
	  m2d2 = m2/d2;
	  *y1 = x1-m1d1*x11-m2d2*x12;
	  *y2 = x2-(m1d1+m2d2)*x2;
	}
}


boolean inside_circ(FLOAT y1, FLOAT y2, FLOAT *q)
{
	FLOAT rsq,rhosq,d1,d2;

	d1 = y1-q[0];
	d2 = y2-q[1];
	rsq = d1*d1+d2*d2;
	rhosq = q[2]*q[2];
	if (rsq <= rhosq)
	  return(TRUE);
	else
	  return(FALSE);
}


FLOAT rho_circ(FLOAT y1, FLOAT y2, FLOAT *q)
{
	FLOAT rsq,d1,d2;

	d1 = y1-q[0];
	d2 = y2-q[1];
	rsq = d1*d1+d2*d2;
	if (q[2] == 0.0) return(1e38);
	return(sqrt(rsq)/q[2]);
}



FLOAT ld_linear(int n, FLOAT gam[], FLOAT rho)
{

	if (rho == 0.) 
	  return 1+.5*gam[0];
	else if (rho < 1.)
	  return(1.0+gam[0]*(1.5*sqrt(1.0-SQ(rho))-1.0));
	else 
	  return 1.-gam[0];
}

FLOAT mag_binext(FLOAT y1, FLOAT y2, FLOAT rho, FLOAT d, FLOAT q, 
	FLOAT (*ld_func)(int,FLOAT*,FLOAT), int n, FLOAT gam[],
	FLOAT acc, FLOAT ld_acc)
{
	PTLIST_PTR ptlist,holelist,newhole;
	FLOAT mag;

	if (rho <= 0.0) { /* point source */
	  mag = mag_binpt(y1,y2,d,q);
	} else {
	  if (q > 1.0) { /* require q > 1.0 for erdlcaust */
	    q = 1.0/q;
	    y1 = -y1;
	  }
	  ptlist = NULL;
	  holelist = NULL;
          sol_binpt(y1,y2,d,q,&ptlist);
	  newhole = (PTLIST_PTR) malloc(sizeof(PTLIST));
	  newhole->x1 = q/(1.0+q)*d;
	  newhole->x2 = 0.0;
	  newhole->min_lvl = 0;
	  newhole->next = holelist;
	  holelist = newhole;
	  newhole = (PTLIST_PTR) malloc(sizeof(PTLIST));
	  newhole->x1 = -1.0/(1.0+q)*d;
	  newhole->x2 = 0.0;
	  newhole->min_lvl = 0;
	  newhole->next = holelist;
	  holelist = newhole;
	  /* attention: produce correct m1, what is left and right? */
	  /* need: m1 < 1, if necessary swith x2 and y2 */
          get_crit_circ(d,q/(1+q),y1,y2,rho,&ptlist);
	  mag = adaptive_contour(acc,ld_acc,
	    &ptlist,&holelist,ld_func,n,gam,
	    rho_circ,inside_circ,lenseq_bin,
	    3,2,y1,y2,rho,d,q);
	  mag /= M_PI*SQ(rho);
	}
	return(mag);
}


FLOAT magt_binext(FLOAT t, FLOAT t0, FLOAT tE, FLOAT alpha, FLOAT u0,
	FLOAT rho, FLOAT d, FLOAT q, 
	FLOAT (*ld_func)(int,FLOAT*,FLOAT), int n, FLOAT gam[],
	FLOAT acc, FLOAT ld_acc)
{
	FLOAT p;
	FLOAT y1,y2;
	FLOAT sa,ca;

	p = (t-t0)/tE;

#ifdef __APPLE__
        __sincos(alpha,&sa,&ca);
#else
	sincos(alpha,&sa,&ca);
#endif
	y1 = p*ca-u0*sa;
	y2 = p*sa+u0*ca;

	return(mag_binext(y1,y2,rho,d,q,ld_func,n,gam,acc,ld_acc)); 
}


