/* BINARY_PT.C */
/* Routines for calculating binary lens magnification by solving
	5th order complex polynomial */
/* C-adaption of routines by B.S. Gaudi */
/* uses "complex.c" of Numerical Recipes */
/* Martin Dominik 19.8.99 */
/* adapted from "gl_fifth.c" 6-3-06 */
/* modified 14-09-06 */
/* modified RP 02-02-18 */



#ifndef FLOAT
#define FLOAT double
#endif

#ifndef boolean
typedef int boolean;
#define FALSE 0
#define TRUE 1
#endif

#ifndef MAX
#define MAX(a,b) ((a) > (b)) ? (a) : (b)
#endif

#define BIGG 1.0e38

typedef struct ptlist {
        FLOAT x1,x2;
	int min_lvl;
        struct ptlist *next;
} PTLIST, *PTLIST_PTR;

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "complex.h"

void sol_binpt(FLOAT y1, FLOAT y2, FLOAT d, FLOAT q, PTLIST_PTR *ptlist);
FLOAT mag_binpt(FLOAT y1, FLOAT y2, FLOAT d, FLOAT q);
static fcomplex lenseq(fcomplex z, fcomplex z1, fcomplex z2, FLOAT m1, FLOAT m2);
static void zroots(fcomplex *a, int m, fcomplex *roots, boolean polish);
static void laguer(fcomplex *a, int m, fcomplex *x, int *its);
static FLOAT dmag(fcomplex z);


void sol_binpt(FLOAT y1, FLOAT y2, FLOAT d, FLOAT q, PTLIST_PTR *ptlist)
{
      fcomplex c[7],z[6],z1sq,z1cu,z1fo;
      fcomplex zeta,zetab, z1, z2,zetabsq;
      FLOAT xsi,eta;
      FLOAT m1,m2,dm,mtot,mtotsq,dmsq;
      boolean polish;
      int m,l;
//       fcomplex err,dzeta,zetap,z1b,z2b,dumz,dumzb;
      fcomplex err, zetap, z1b, z2b;
      FLOAT cm;
      FLOAT pi;
//       FLOAT dzetaabs; /* new */
      PTLIST_PTR succ;

      *ptlist = NULL;

      pi=4.0e0*atan(1.0e0);

      m2=q/(1.0+q);
      m1=1.0/(1.0+q);
      dm=(m2-m1)/2.0;
      mtot=(m1+m2)/2.0;
      cm=-d/2.0*(dm/mtot);
            
      xsi=y1+cm;
      eta=y2; 

     
      /* convert to quantities in Witt & Mao paper */

      z1=Complex(d/2.0,0.0);
      z2=Complex(-d/2.0,0.0);
      zeta=Complex(xsi,eta);
      zetab=Complex(xsi,-eta);
      z1sq=Cmul(z1,z1);
      z1cu=Cmul(z1sq,z1);
      z1fo=Cmul(z1cu,z1);
      z1b=z1;
      mtotsq=mtot*mtot;
      dmsq=dm*dm;
      z2b=z2;
      zetabsq=Cmul(zetab,zetab);

      /* Calculate coefficients */

      c[6] = Csub(z1sq,zetabsq);
      c[5] = Cadd(RCmul(-2.0*mtot,zetab),Cmul(zeta,zetabsq));
      c[5] = Csub(c[5],RCmul(2.0*dm,z1));
      c[5] = Csub(c[5],Cmul(zeta,z1sq));
      c[4] = Cmul(RCmul(4.0*mtot,zeta),zetab);
      c[4] = Cadd(c[4],Cmul(RCmul(4.0*dm,zetab),z1));
      c[4] = Cadd(c[4],Cmul(RCmul(2.0,zetabsq),z1sq));
      c[4] = Csub(c[4],RCmul(2.0,z1fo));
      c[3] = Cadd(RCmul(4.0*mtotsq,zeta),RCmul(4.0*mtot*dm,z1));
      c[3] = Csub(c[3],Cmul(Cmul(RCmul(4.0*dm,zeta),zetab),z1));
      c[3] = Csub(c[3],Cmul(Cmul(RCmul(2.0,zeta),zetabsq),z1sq));
      c[3] = Cadd(c[3],RCmul(4.0*dm,z1cu));
      c[3] = Cadd(c[3],Cmul(RCmul(2.0,zeta),z1fo));
      c[2] = Cmul(RCmul(-8.0*mtot*dm,zeta),z1);
      c[2] = Csub(c[2],RCmul(4.0*dmsq,z1sq));
      c[2] = Csub(c[2],RCmul(4.0*mtotsq,z1sq)); 
      c[2] = Csub(c[2],Cmul(Cmul(RCmul(4.0*mtot,zeta),zetab),z1sq));
      c[2] = Csub(c[2],Cmul(RCmul(4.0*dm,zetab),z1cu));
      c[2] = Csub(c[2],Cmul(zetabsq,z1fo));
      c[2] = Cadd(c[2],Cmul(z1fo,z1sq));
      c[1] = Cadd(RCmul(4.0*dmsq,zeta),RCmul(4.0*mtot*dm,z1));
      c[1] = Cadd(c[1],Cmul(Cmul(RCmul(4.0*dm,zeta),zetab),z1));
      c[1] = Cadd(c[1],Cmul(RCmul(2.0*mtot,zetab),z1sq));
      c[1] = Cadd(c[1],Cmul(Cmul(zeta,zetabsq),z1sq));
      c[1] = Csub(c[1],RCmul(2.0*dm,z1cu));
      c[1] = Csub(c[1],Cmul(zeta,z1fo));
      c[1] = Cmul(c[1],z1sq); 

      m=5;


      polish=TRUE;

      zroots(c,m,z,polish);

         

      /* test to see if roots are real */
      for (l=1; l<=5; l++) { 
      	   zetap = lenseq(z[l],z1,z2,m1,m2); 
           err = Csub(zeta,zetap);
      	   if(dmag(err) <= 1.0e-4) { 
	     succ = *ptlist;
	     *ptlist = (PTLIST_PTR) malloc(sizeof(PTLIST));
	     (*ptlist)->x1 = z[l].r-cm;
	     (*ptlist)->x2 = z[l].i;
	     (*ptlist)->min_lvl = 0; /* 14-9-06 */
	     (*ptlist)->next = succ;
	   }
      } 

}


FLOAT mag_binpt(FLOAT y1, FLOAT y2, FLOAT d, FLOAT q)
{
	PTLIST_PTR ptlist,succ;
	FLOAT m1,m2,dm,mtot,cm;
	FLOAT detJ,mag;
	fcomplex z,zb,z1,z1b,dzeta;
	FLOAT dzetaabs;

        m2=q/(1.0+q);
        m1=1.0/(1.0+q);
        dm=(m2-m1)/2.0;
        mtot=(m1+m2)/2.0;
        cm=-d/2.0*(dm/mtot);

	sol_binpt(y1,y2,d,q,&ptlist);

	mag = 0.0;
	while (ptlist != NULL) {
	  z = Complex(ptlist->x1+cm,ptlist->x2);
	  zb = Conjg(z);
	  z1 = Complex(d/2.0,0.0);
	  z1b = Conjg(z1);
	  dzeta = Cdiv(Cdiv(Complex(mtot-dm,0),Csub(z1b,zb)),
			Csub(z1b,zb));
	  dzeta = Cadd(dzeta,Cdiv(Cdiv(Complex(mtot+dm,0),
			Csub(RCmul(-1.0,z1b),zb)),
	 	        Csub(RCmul(-1.0,z1b),zb)));
	  dzetaabs = Cabs(dzeta);
	  detJ = 1.0-dzetaabs*dzetaabs; /* z*z_conj = |z|^2 */
	  if (detJ == 0.0) mag += BIGG;
	   else mag += fabs(1.0/detJ);
	  succ = ptlist->next;
	  free(ptlist);
	  ptlist = succ;
	}
	return(mag);
}


      

static fcomplex lenseq(fcomplex z, fcomplex z1, fcomplex z2, FLOAT m1, FLOAT m2)
{      
      fcomplex z1b,z2b,zb;
      fcomplex zetap;

      z1b = Conjg(z1);
      z2b = Conjg(z2);
      zb = Conjg(z);
      
      zetap = Cadd(z,Cdiv(Complex(m1,0),Csub(z1b,zb)));
      zetap = Cadd(zetap,Cdiv(Complex(m2,0),Csub(z2b,zb)));
      return(zetap);
}

      
     
#define MAXM 501
#define EPS 1.0e-12
     
static void zroots(fcomplex *a, int m, fcomplex *roots, boolean polish)
{
      int i,j,jj,its;
      fcomplex ad[MAXM+1],x,b,c;

      for (j=1; j<=m+1; j++)
           ad[j]=a[j];
      for (j=m; j>=1; j--) {
           x=Complex(0.0,0.0);
           laguer(ad,j,&x,&its);
      	   if(fabs(x.i) <= 2.0*EPS*EPS*fabs(x.r)) x = Complex(x.r,0.0);
           roots[j] = x;
           b = ad[j+1];
           for (jj=j; jj>=1; jj--) {
                c=ad[jj];
                ad[jj]=b;
                b=Cadd(Cmul(x,b),c);
	   }
      }

      if (polish) { 
           for (j=1; j<=m; j++) 
                laguer(a,m,&roots[j],&its);
      }
      for (j=2; j<=m; j++) {
           x = roots[j];
       	   for (i=j-1; i>=1; i--) {
                if(roots[i].r <= x.r) break;
                roots[i+1]=roots[i];
           }
           roots[i+1]=x;
      } 
}

#undef EPS
#undef MAXM

      
#define EPSS 1.0e-12
#define MR 8
#define MT 19000L
#define MAXIT (MT*MR)

static void laguer(fcomplex *a, int m, fcomplex *x, int *its)
{
      int j;
      unsigned long iter;
      FLOAT abx,abp,abm,err;
      fcomplex dx,x1,b,d,f,g,h,sq,gp,gm,g2;
      static FLOAT frac[MR+1] = {0,0.5,0.25,0.75,0.13,0.38,0.62,0.88,1.0};

      for (iter=1; iter<=MAXIT; iter++) { 
           *its=iter;
           b=a[m+1];
           err=Cabs(b);
           d=Complex(0.0,0.0);
           f=Complex(0.0,0.0);
           abx=Cabs(*x);
           for (j=m; j>=1; j--) {
               f=Cadd(Cmul(*x,f),d);
               d=Cadd(Cmul(*x,d),b);
               b=Cadd(Cmul(*x,b),a[j]);
               err=Cabs(b)+abx*err;
           } 
           err=EPSS*err;
           if(Cabs(b) <= err) return;
           g=Cdiv(d,b);
           g2=Cmul(g,g);
           h=Csub(g2,RCmul(2.0,Cdiv(f,b)));
           sq=Csqrt(RCmul(m-1,Csub(RCmul(m,h),g2)));
           gp=Cadd(g,sq);
           gm=Csub(g,sq);
           abp=Cabs(gp);
           abm=Cabs(gm);
           if (abp <= abm) gp=gm;
           if (MAX(abp,abm)>=0) 
                dx=Cdiv(Complex(m,0),gp);
           else
                dx=RCmul(exp(log(1.+abx)),Complex(cos(iter),sin(iter)));
           x1=Csub(*x,dx);
           if(x->r == x1.r && x->i == x1.i) return;
           if ((iter % MT) !=0) 
               *x=x1;
           else
               *x=Csub(*x,RCmul(frac[iter/MT],dx));
      } 
}

#undef MAXIT
#undef MR
#undef MT
#undef EPSS

      
static FLOAT dmag(fcomplex z)
{
      return(fabs(z.r)+fabs(z.i));
}


