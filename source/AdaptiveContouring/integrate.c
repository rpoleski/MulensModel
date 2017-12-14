/* integrate.c */

#include <math.h>

#ifndef FLOAT
#define FLOAT double
#endif


FLOAT *vector(int nl,int nh);
int *ivector(int nl,int nh);
FLOAT **matrix(int nrl,int nrh,int ncl,int nch);
void free_vector(FLOAT *v,int nl,int nh);
void free_ivector(int *v,int nl,int nh);
void free_matrix(FLOAT **m,int nrl,int nrh,int ncl,int nch);

/* Adaptive Simpson integration, based on Algorithm 182 in CACM
	by W.H. McKeeman */
	
/* WARNING: error estimate only valid if integrand is either
 	> 0 in [a,b] or <0 in [a,b] */

/* adapsimpson1f: integrate over first of two variables  */
/* use free limb-darkening function as parameter */ 

#define MAXLEVEL 30
	
FLOAT adapsimpsonf1(FLOAT (*f)(FLOAT,FLOAT,
          FLOAT (*)(), FLOAT*, void (*)(), FLOAT*,
          FLOAT (*)(int,FLOAT*,FLOAT),int,FLOAT*),
        FLOAT a, FLOAT b, FLOAT p2,
        FLOAT (*rho_func)(), FLOAT *ps,
        void (*lenseq_func)(), FLOAT *pl,
        FLOAT (*ld_func)(int,FLOAT*,FLOAT),
        int n, FLOAT *gam, FLOAT tol)
{
	int l;
	FLOAT *dx,*x2,*x3,*f2,*f3,*f4,*fmp,*fbp,*est2,*est3,**pval;
	int *rtrn;
	FLOAT est,fa,fm,fb,da,sx,est1,sum,f1,absarea;
		
	if (b-a == 0.0) return(0);  /* 18.8.95 */	
	
	dx = vector(1,MAXLEVEL);
	x2 = vector(1,MAXLEVEL);
	x3 = vector(1,MAXLEVEL);
	f2 = vector(1,MAXLEVEL);
	f3 = vector(1,MAXLEVEL);
	f4 = vector(1,MAXLEVEL);
	fmp = vector(1,MAXLEVEL);
	fbp = vector(1,MAXLEVEL);
	est2 = vector(1,MAXLEVEL);
	est3 = vector(1,MAXLEVEL);
	rtrn = ivector(1,MAXLEVEL);
	pval = matrix(1,MAXLEVEL,1,3);
	
	
	l = 0;
	est = 0;
	absarea = 0;
	da = b-a;
	fa = (*f)(a,p2,rho_func,ps,lenseq_func,pl,ld_func,n,gam);
	fm = 4.0*(*f)((a+b)/2.0,p2,rho_func,ps,lenseq_func,pl,ld_func,n,gam);
	fb = (*f)(b,p2,rho_func,ps,lenseq_func,pl,ld_func,n,gam);
	
	recur:
		l++;
		dx[l] = da/3.0;
		sx = dx[l]/6.0;
		f1 = 4.0*(*f)(a+dx[l]/2.0,p2,rho_func,ps,lenseq_func,pl,ld_func,n,gam);
		x2[l] = a+dx[l];
		f2[l] = (*f)(x2[l],p2,rho_func,ps,lenseq_func,pl,ld_func,n,gam);
		x3[l] = x2[l]+dx[l];
		f3[l] = (*f)(x3[l],p2,rho_func,ps,lenseq_func,pl,ld_func,n,gam);
		f4[l] = 4.0*(*f)(x3[l]+dx[l]/2.0,p2,rho_func,ps,lenseq_func,pl,ld_func,n,gam);
		fmp[l] = fm;
		fbp[l] = fb;
		est1 = (fa+f1+f2[l])*sx;
		est2[l] = (f2[l]+f3[l]+fm)*sx;
		est3[l] = (f3[l]+f4[l]+fb)*sx;
		sum = est1+est2[l]+est3[l];
		absarea = absarea-fabs(est)+fabs(est1)+fabs(est2[l])
			+fabs(est3[l]);
		if (fabs(est-sum) <= tol*fabs(absarea) || l>=MAXLEVEL) {
			up:
				if (l==1) goto r4; /* 23.9.97 */
				l--;
				pval[l][rtrn[l]] = sum;
				switch(rtrn[l]) {
					case 1: goto r1;
					case 2: goto r2;
					case 3: goto r3;
				}
		}
		rtrn[l] = 1;
		da = dx[l];
		fm = f1;
		fb = f2[l];
		est = est1;
		goto recur;
		
		r1:
		rtrn[l] = 2;
		da = dx[l];
		fa = f2[l];
		fm = fmp[l];
		fb = f3[l];
		est = est2[l];
		a = x2[l];
		goto recur;
		
		r2:
		rtrn[l] = 3;
		da = dx[l];
		fa = f3[l];
		fm = f4[l];
		fb = fbp[l];
		est = est3[l];
		a = x3[l];
		goto recur;
		
		r3:
		sum = pval[l][1]+pval[l][2]+pval[l][3];
		if (l>1) goto up;
		
	r4:  /* 23.9.97 */
	free_vector(dx,1,MAXLEVEL);
	free_vector(x2,1,MAXLEVEL);
	free_vector(x3,1,MAXLEVEL);
	free_vector(f2,1,MAXLEVEL);
	free_vector(f3,1,MAXLEVEL);
	free_vector(f4,1,MAXLEVEL);
	free_vector(fmp,1,MAXLEVEL);
	free_vector(fbp,1,MAXLEVEL);
	free_vector(est2,1,MAXLEVEL);
	free_vector(est3,1,MAXLEVEL);
	free_ivector(rtrn,1,MAXLEVEL);
	free_matrix(pval,1,MAXLEVEL,1,3);
	return(sum);
	
}


/* Adaptive Simpson integration, based on Algorithm 182 in CACM
	by W.H. McKeeman */
	
/* WARNING: error estimate only valid if integrand is either
 	> 0 in [a,b] or <0 in [a,b] */

/* adapsimpson2f: integrate over second of two variables  */
/* use free limb-darkening function as parameter */ 

#define MAXLEVEL 30
	
FLOAT adapsimpsonf2(FLOAT (*f)(FLOAT,FLOAT,
          FLOAT (*)(), FLOAT*, void (*)(), FLOAT*,
          FLOAT (*)(int,FLOAT*,FLOAT),int,FLOAT*),
	FLOAT p1, FLOAT a, FLOAT b,
        FLOAT (*rho_func)(), FLOAT *ps,
        void (*lenseq_func)(), FLOAT *pl,
	FLOAT (*ld_func)(int,FLOAT*,FLOAT),
	int n, FLOAT *gam, FLOAT tol)
{
	int l;
	FLOAT *dx,*x2,*x3,*f2,*f3,*f4,*fmp,*fbp,*est2,*est3,**pval;
	int *rtrn;
	FLOAT est,fa,fm,fb,da,sx,est1,sum,f1,absarea;
		
	if (b-a == 0.0) return(0);  /* 18.8.95 */	
	
	dx = vector(1,MAXLEVEL);
	x2 = vector(1,MAXLEVEL);
	x3 = vector(1,MAXLEVEL);
	f2 = vector(1,MAXLEVEL);
	f3 = vector(1,MAXLEVEL);
	f4 = vector(1,MAXLEVEL);
	fmp = vector(1,MAXLEVEL);
	fbp = vector(1,MAXLEVEL);
	est2 = vector(1,MAXLEVEL);
	est3 = vector(1,MAXLEVEL);
	rtrn = ivector(1,MAXLEVEL);
	pval = matrix(1,MAXLEVEL,1,3);
	
	
	l = 0;
	est = 0;
	absarea = 0;
	da = b-a;
	fa = (*f)(p1,a,rho_func,ps,lenseq_func,pl,ld_func,n,gam);
	fm = 4.0*(*f)(p1,(a+b)/2.0,rho_func,ps,lenseq_func,pl,ld_func,n,gam);
	fb = (*f)(p1,b,rho_func,ps,lenseq_func,pl,ld_func,n,gam);
	
	recur:
		l++;
		dx[l] = da/3.0;
		sx = dx[l]/6.0;
		f1 = 4.0*(*f)(p1,a+dx[l]/2.0,rho_func,ps,lenseq_func,pl,ld_func,n,gam);
		x2[l] = a+dx[l];
		f2[l] = (*f)(p1,x2[l],rho_func,ps,lenseq_func,pl,ld_func,n,gam);
		x3[l] = x2[l]+dx[l];
		f3[l] = (*f)(p1,x3[l],rho_func,ps,lenseq_func,pl,ld_func,n,gam);
		f4[l] = 4.0*(*f)(p1,x3[l]+dx[l]/2.0,rho_func,ps,lenseq_func,pl,ld_func,n,gam);
		fmp[l] = fm;
		fbp[l] = fb;
		est1 = (fa+f1+f2[l])*sx;
		est2[l] = (f2[l]+f3[l]+fm)*sx;
		est3[l] = (f3[l]+f4[l]+fb)*sx;
		sum = est1+est2[l]+est3[l];
		absarea = absarea-fabs(est)+fabs(est1)+fabs(est2[l])
			+fabs(est3[l]);
		if (fabs(est-sum) <= tol*fabs(absarea) || l>=MAXLEVEL) {
			up:
				if (l==1) goto r4; /* 23.9.97 */
				l--;
				pval[l][rtrn[l]] = sum;
				switch(rtrn[l]) {
					case 1: goto r1;
					case 2: goto r2;
					case 3: goto r3;
				}
		}
		rtrn[l] = 1;
		da = dx[l];
		fm = f1;
		fb = f2[l];
		est = est1;
		goto recur;
		
		r1:
		rtrn[l] = 2;
		da = dx[l];
		fa = f2[l];
		fm = fmp[l];
		fb = f3[l];
		est = est2[l];
		a = x2[l];
		goto recur;
		
		r2:
		rtrn[l] = 3;
		da = dx[l];
		fa = f3[l];
		fm = f4[l];
		fb = fbp[l];
		est = est3[l];
		a = x3[l];
		goto recur;
		
		r3:
		sum = pval[l][1]+pval[l][2]+pval[l][3];
		if (l>1) goto up;
		
	r4:  /* 23.9.97 */
	free_vector(dx,1,MAXLEVEL);
	free_vector(x2,1,MAXLEVEL);
	free_vector(x3,1,MAXLEVEL);
	free_vector(f2,1,MAXLEVEL);
	free_vector(f3,1,MAXLEVEL);
	free_vector(f4,1,MAXLEVEL);
	free_vector(fmp,1,MAXLEVEL);
	free_vector(fbp,1,MAXLEVEL);
	free_vector(est2,1,MAXLEVEL);
	free_vector(est3,1,MAXLEVEL);
	free_ivector(rtrn,1,MAXLEVEL);
	free_matrix(pval,1,MAXLEVEL,1,3);
	return(sum);
	
}


