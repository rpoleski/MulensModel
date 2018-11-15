/* erdlcaust.c */
/* V1.00 31-05-06 MD */
/* debugged 30/8/06 MD */
/* modified 14-9-06 MD */
/* modified 1-2-07 MD */
/* modified 10-11-07 MD */
/* modified 02-02-18 RP */
/* modified 24-02-18 RP */
/* modified 15-11-18 RP */


#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define CUSP_PREC 1.0e-10
#define RTPREC 1.0e-14
#define LIMPREC 1.0e-15

#define FLOAT double

#define SQ(a) ((a)*(a))

typedef enum{FALSE,TRUE} boolean;

typedef struct ptlist {
        FLOAT x1,x2;
	int min_lvl;
        struct ptlist *next;
} PTLIST, *PTLIST_PTR;


void get_crit_circ(FLOAT d, FLOAT m1, FLOAT y1c, FLOAT y2c, FLOAT rcirc,
	PTLIST_PTR *crit_points);

static void add_ptlist(PTLIST_PTR *list, FLOAT x1, FLOAT x2);
static FLOAT bisection3(FLOAT (*func)(FLOAT,FLOAT,FLOAT),
  FLOAT x1, FLOAT x2, FLOAT p1, FLOAT p2, FLOAT acc);
static FLOAT bisection2(FLOAT (*func)(FLOAT,FLOAT),
  FLOAT x1, FLOAT x2, FLOAT p1, FLOAT acc);
static FLOAT get_rootdef(FLOAT (*func)(FLOAT,boolean,FLOAT,FLOAT), 
	FLOAT rmin, FLOAT rmax,
	boolean pmflag, FLOAT d, FLOAT m1, int signupper, FLOAT acc);
static FLOAT get_rootdefm(FLOAT (*func)(FLOAT,boolean,FLOAT,FLOAT,FLOAT,FLOAT), 
	FLOAT rmin, FLOAT rmax,
	boolean pmflag, FLOAT d, FLOAT m1, 
	FLOAT y1c, FLOAT y2c,
	int signupper, FLOAT acc);
static void get_rootdef_lu(FLOAT (*func)(FLOAT,boolean,FLOAT,FLOAT), 
	FLOAT rmin, FLOAT rmax,
	boolean pmflag, FLOAT d, FLOAT m1, int signupper, 
	FLOAT *lower, FLOAT *upper, FLOAT acc);
static FLOAT grid_brack(FLOAT (*func)(FLOAT,boolean,FLOAT,FLOAT), 
	FLOAT rmin, FLOAT rmax,
	boolean pmflag, FLOAT d, FLOAT m1, int sign, FLOAT acc);
static FLOAT erdl_dp(FLOAT r, FLOAT d, FLOAT m1);
static FLOAT erdl_pm(FLOAT r, FLOAT d, FLOAT m1);
static FLOAT erdl_pp(FLOAT r, FLOAT d, FLOAT m1);
static void erdl_coeff(FLOAT r, FLOAT d, FLOAT m1, 
	FLOAT *a0, FLOAT *a1, FLOAT *a2);
static FLOAT cosplus(FLOAT r, FLOAT d, FLOAT m1);
static FLOAT cosminus(FLOAT r, FLOAT d, FLOAT m1);
static FLOAT get_cos(FLOAT r, boolean pmflag, FLOAT d, FLOAT m1);
static FLOAT calc_d1(FLOAT m1,FLOAT acc);
static FLOAT d1func(FLOAT d, FLOAT m1);
static boolean close_binary(FLOAT d, FLOAT m1, FLOAT r1);
static boolean wide_binary(FLOAT d, FLOAT m1, FLOAT r2);
static void erdl_limits(FLOAT d, FLOAT m1, FLOAT *r_d1, FLOAT *r_d2,
	FLOAT *r_d3, FLOAT *r_pm, FLOAT *r_pp1, FLOAT *r_pp2, FLOAT *r_pp3);
static FLOAT perp_caust(FLOAT r, boolean pmflag, FLOAT d, FLOAT m1,
        FLOAT y1c, FLOAT y2c);
static FLOAT perp_caustd1(FLOAT r, boolean pmflag, FLOAT d, FLOAT m1,
        FLOAT y1c, FLOAT y2c);
static FLOAT perp_caustd2(FLOAT r, boolean pmflag, FLOAT d, FLOAT m1,
        FLOAT y1c, FLOAT y2c);
static FLOAT curvature_der(FLOAT r, boolean pmflag, FLOAT d, FLOAT m1);
static void add_to_list(PTLIST_PTR *list, FLOAT r, boolean pmflag,
	FLOAT d, FLOAT m1, boolean invert);
static void add_crit(FLOAT rmin, FLOAT rmax, FLOAT ty1min, FLOAT ty2min,
	FLOAT ty1max, FLOAT ty2max, boolean pmflag,
	FLOAT d, FLOAT m1, FLOAT y1c, FLOAT y2c, FLOAT rcirc,
	PTLIST_PTR *crit_points);
static void x1x2 (FLOAT r, FLOAT cval, FLOAT d, FLOAT m1,
   FLOAT *x1, FLOAT *x2);
static FLOAT offdiag(FLOAT r, boolean pmflag, FLOAT d, FLOAT m1);
static FLOAT cspf(FLOAT r, boolean pmflag, FLOAT d, FLOAT m1);
static void normal(FLOAT x1, FLOAT x2, FLOAT d, FLOAT m1, FLOAT *n1, FLOAT *n2);
static void tangent_y(FLOAT x1, FLOAT x2, FLOAT d, FLOAT m1, 
	FLOAT *ty1, FLOAT *ty2);
static void binlens(FLOAT x1, FLOAT x2, FLOAT d, FLOAT m1, 
	FLOAT *y1, FLOAT *y2);
static void binlensd1(FLOAT x1, FLOAT x2, FLOAT d, FLOAT m1, 
	FLOAT *y11, FLOAT *y12, FLOAT *y22);
static void binlensd2(FLOAT x1, FLOAT x2, FLOAT d, FLOAT m1, 
	FLOAT *y111, FLOAT *y112, FLOAT *y122, FLOAT *y222);
static void binlensd3(FLOAT x1, FLOAT x2, FLOAT d, FLOAT m1, 
	FLOAT *y1111, FLOAT *y1112, FLOAT *y1122, FLOAT *y1222, FLOAT *y2222);
static void binlensd4(FLOAT x1, FLOAT x2, FLOAT d, FLOAT m1, 
	FLOAT *y11111, FLOAT *y11112, FLOAT *y11122, 
	FLOAT *y11222, FLOAT *y12222, FLOAT *y22222);
static void binlensd5(FLOAT x1, FLOAT x2, FLOAT d, FLOAT m1, 
	FLOAT *y111111, FLOAT *y111112, FLOAT *y111122, 
	FLOAT *y111222, FLOAT *y112222, FLOAT *y122222, FLOAT *y222222);


#define JMAX 5000

static void add_ptlist(PTLIST_PTR *list, FLOAT x1, FLOAT x2)
{
	PTLIST_PTR new_;

	new_ = (PTLIST_PTR) malloc(sizeof(PTLIST));
	if (new_ == NULL) {
	  fprintf(stderr,"Sorry, not enough memory available\n");
	  exit(1);
	}
	new_->x1 = x1;
	new_->x2 = x2;
	new_->min_lvl = 0;
	new_->next = *list;
	*list = new_;
}


#define MAXCNT 10000

static FLOAT bisection3(FLOAT (*func)(FLOAT,FLOAT,FLOAT),
  FLOAT x1, FLOAT x2, FLOAT p1, FLOAT p2, FLOAT acc)
{
	int j;
	FLOAT val, val2;
	FLOAT diff, mid; 
	FLOAT root;

	val = (*func)(x1,p1,p2);
	val2 = (*func)(x2,p1,p2);
	if (val == 0.0) return(x1);
	if (val2 == 0.0) return(x2);
	if (val*val2 > 0.0) {
	   fprintf(stderr,"Function not bracketed in bisection2");
	   exit(1);
	}
        if (val < 0.0) {
	  root = x1;
	  diff = x2-x1;
        } else {
	  root = x2;
	  diff = x1-x2;
        } 
	for (j=1; j<=MAXCNT; j++) {
	  diff /= 2.0;
	  mid = root+diff; 
	  val2 = (*func)(mid,p1,p2);
	  if (val2 <= 0.0) 
	    root = mid;
	  if (fabs(diff) < acc || val2 == 0.0) 
	    return(root);
	}
	fprintf(stderr,"Too many steps in bisection2");
	exit(1);
}

#undef MAXCNT 

#define MAXCNT 10000

static FLOAT bisection2(FLOAT (*func)(FLOAT,FLOAT),
  FLOAT x1, FLOAT x2, FLOAT p1, FLOAT acc)
{
	int j;
	FLOAT val, val2;
	FLOAT diff, mid; 
	FLOAT root;

	val = (*func)(x1,p1);
	val2 = (*func)(x2,p1);
	if (val == 0.0) return(x1);
	if (val2 == 0.0) return(x2);
	if (val*val2 > 0.0) {
	   fprintf(stderr,"Function not bracketed in bisection2");
	   exit(1);
	}
        if (val < 0.0) {
	  root = x1;
	  diff = x2-x1;
        } else {
	  root = x2;
	  diff = x1-x2;
        } 
	for (j=1; j<=MAXCNT; j++) {
	  diff /= 2.0;
	  mid = root+diff; 
	  val2 = (*func)(mid,p1);
	  if (val2 <= 0.0) 
	    root = mid;
	  if (fabs(diff) < acc || val2 == 0.0) 
	    return(root);
	}
	fprintf(stderr,"Too many steps in bisection2");
	exit(1);
}

#undef MAXCNT 

#define MAXCNT 10000

static FLOAT get_rootdef(FLOAT (*func)(FLOAT,boolean,FLOAT,FLOAT), FLOAT rmin, FLOAT rmax,
	boolean pmflag, FLOAT d, FLOAT m1, int signupper, FLOAT acc)
{
        int j;
	FLOAT val;
	FLOAT mid; 
   
	/* acc *= fabs(rmax-rmin); */ /* fractional accuracy in x */
        for (j=1; j<=MAXCNT; j++) {
	  mid = (rmax+rmin)/2.0;
	  val = (*func)(mid,pmflag,d,m1);
          if (val == 0.0) return(mid);
	  if (val*signupper > 0) 
	    rmax = mid;
	  else
            rmin = mid;
	  if (fabs(rmax-rmin) < acc) 
	    return(mid);
        }
        fprintf(stderr,"Too many steps in get_rootdef");
        exit(1);
}

#undef MAXCNT

#define MAXCNT 10000

static FLOAT get_rootdefm(FLOAT (*func)(FLOAT,boolean,FLOAT,FLOAT,FLOAT,FLOAT), 
	FLOAT rmin, FLOAT rmax,
	boolean pmflag, FLOAT d, FLOAT m1, 
	FLOAT y1c, FLOAT y2c, int signupper, FLOAT acc)
{
        int j;
	FLOAT val;
	FLOAT mid; 

   
	/* acc *= fabs(rmax-rmin); */ /* fractional accuracy in x */
        for (j=1; j<=MAXCNT; j++) {
	  mid = (rmax+rmin)/2.0;
	  val = (*func)(mid,pmflag,d,m1,y1c,y2c);
          if (val == 0.0) return(mid);
	  if (val*signupper > 0) 
	    rmax = mid;
	  else
            rmin = mid;
	  if (fabs(rmax-rmin) < acc) 
	    return(mid);
        }
        fprintf(stderr,"Too many steps in get_rootdefm");
        exit(1);
}

#undef MAXCNT

#define MAXCNT 10000

static void get_rootdef_lu(FLOAT (*func)(FLOAT,boolean,FLOAT,FLOAT), 
	FLOAT rmin, FLOAT rmax,
	boolean pmflag, FLOAT d, FLOAT m1, int signupper, 
	FLOAT *lower, FLOAT *upper, FLOAT acc)
{
        int j;
	FLOAT mid,val;

   
	/* acc *= fabs(rmax-rmin); */ /* fractional accuracy in x */
        for (j=1; j<=MAXCNT; j++) {
	  mid = (rmax+rmin)/2.0;
	  val = (*func)(mid,pmflag,d,m1);
          if (val == 0.0) {
	    *lower = mid-acc*fabs(rmax-rmin);
	    *upper = mid+acc*fabs(rmax-rmin);
	    return;
	  }
	  if (val*signupper > 0) 
	    rmax = mid;
	  else
            rmin = mid;
	  if (fabs(rmax-rmin) < acc) {
	    *lower = rmin;
	    *upper = rmax;
	    return;
	  }
        }
        fprintf(stderr,"Too many steps in get_rootdef_lu");
        exit(1);
}

#undef MAXCNT

#define DEL 1.0e-8 

static FLOAT grid_brack(FLOAT (*func)(FLOAT,boolean,FLOAT,FLOAT), FLOAT rmin, FLOAT rmax,
	boolean pmflag, FLOAT d, FLOAT m1, int sign, FLOAT acc)
{
	FLOAT diff,delta;
	FLOAT rquery;
	FLOAT f;

	/* search for point where function has required sign */
  
	diff = rmax-rmin;
	delta = diff;
	do {
	  rquery = rmin+delta/2.0;    
	  do {
            f=(*func)(rquery,pmflag,d,m1);
	    if (f*sign > 0) return(rquery);
	    rquery += delta;
	  } while (rquery < rmax);
	  delta /= 2.0;
	} while (fabs(delta/diff) > DEL); 
        fprintf(stderr,"Too many trials in grid_brack");
        exit(1);
}

#undef DEL 

static void erdl_coeff(FLOAT r, FLOAT d, FLOAT m1, 
	FLOAT *a0, FLOAT *a1, FLOAT *a2)
{
	FLOAT G,H;
	FLOAT g,h;
	FLOAT rsq,dsq;
	FLOAT m1sq;
	FLOAT m2;

	rsq = SQ(r);
	dsq = SQ(d);
	m1sq = SQ(m1);
	m2 = 1.0-m1;

	if (m2 < 0.0) {
	  m1 = 1.0;
	  m2 = 0.0;
	}
	else if (m2 > 1.0) {
	  m1 = 0.0;
	  m2 = 1.0;
	}
	G = rsq + dsq;	
	H = 2.0*d*r;
	g = 1.0-m1sq/SQ(rsq);
	h = 2.0*m1*m2/rsq;

	*a2 = SQ(H)*g-2.0*dsq*h;
	*a1 = H*(h-2*G*g);
	*a0 = G*(G*g-h)+2.0*dsq*h-SQ(m2);	
}

static FLOAT erdl_dp(FLOAT r, FLOAT d, FLOAT m1)
{
	FLOAT a0,a1,a2;
	
	erdl_coeff(r,d,m1,&a0,&a1,&a2);
	return(SQ(a1)-4.0*a0*a2);
}

static FLOAT erdl_pm(FLOAT r, FLOAT d, FLOAT m1)
{
	FLOAT a0,a1,a2;
	
	erdl_coeff(r,d,m1,&a0,&a1,&a2);
	return(a0-a1+a2);
}

static FLOAT erdl_pp(FLOAT r, FLOAT d, FLOAT m1)
{
	FLOAT a0,a1,a2;
	
	erdl_coeff(r,d,m1,&a0,&a1,&a2);
	return(a0+a1+a2);
}

static FLOAT cosplus(FLOAT r, FLOAT d, FLOAT m1)
{
	FLOAT a0,a1,a2;
	FLOAT det,rtdet;
	FLOAT val;

	erdl_coeff(r,d,m1,&a0,&a1,&a2);
	/* 30/8/06 */
	if (a2 == 0.0) {
	  if (a1 == 0.0)
	    val = 1.0;
	  else 
	    val = -a0/a1;
	}
	/* end 30/8/06 */
	else {
	  det = a1*a1-4.0*a0*a2;
	  if (det < 0.0) det = 0.0;
	  rtdet = sqrt(det);
	  val = -0.5*(a1+rtdet)/a2;
	}
	return(val); 
}

static FLOAT cosminus(FLOAT r, FLOAT d, FLOAT m1)
{
	FLOAT a0,a1,a2;
	FLOAT det,rtdet;
	FLOAT val;

	erdl_coeff(r,d,m1,&a0,&a1,&a2);
	/* 30/8/06 */
	if (a2 == 0.0) {
	  if (a1 == 0.0)
	    val = 1.0;
	  else 
	    val = -a0/a1;
	}
	/* end 30/8/06 */
	else {
	  det = a1*a1-4.0*a0*a2;
	  if (det < 0.0) det = 0.0;
	  rtdet = sqrt(det);
	  val = -0.5*(a1-rtdet)/a2;
	}
	return(val); 
}

static FLOAT get_cos(FLOAT r, boolean pmflag, FLOAT d, FLOAT m1)
{
	FLOAT cval;

        if (pmflag)	
          cval = cosplus(r,d,m1);
        else
          cval = cosminus(r,d,m1);
        if (cval < -1.0) cval = -1.0;
        if (cval > 1.0) cval = 1.0;
	return(cval);
}


static FLOAT d1func(FLOAT d, FLOAT m1)
{
	FLOAT t;

	t = d*d;
	t = t*t;
	t = (1.0-t)/3.0;
	t = t*t*t/(m1*(1.0-m1));
	t = pow(t,0.125);
	return(t-d);
}

static FLOAT calc_d1(FLOAT m1,FLOAT acc)
{
	return(bisection2(d1func,0.0,1.0,m1,acc));
}

static boolean close_binary(FLOAT d, FLOAT m1, FLOAT r1)
{
	FLOAT rtm1,rt4m1;

	rtm1 = sqrt(m1);
	rt4m1 = sqrt(rtm1);

	if (erdl_dp(r1,d,m1) * erdl_dp(rtm1*d,d,m1) >= 0.0
           || erdl_dp(r1,d,m1) * erdl_dp(rt4m1,d,m1) >= 0.0) 
	  return(FALSE);
	else
	  return(TRUE); 
}

static boolean wide_binary(FLOAT d, FLOAT m1, FLOAT r2)
{
	FLOAT rtm1,rt4m1;

	rtm1 = sqrt(m1);
	rt4m1 = sqrt(rtm1);

	if (erdl_pp(r2,d,m1) * erdl_pp(d,d,m1) >= 0.0
           || erdl_pp(r2,d,m1) * erdl_pp(rtm1,d,m1) >= 0.0) 
	  return(FALSE);
	else
	  return(TRUE); 
}


static void erdl_limits(FLOAT d, FLOAT m1, FLOAT *r_d1, FLOAT *r_d2,
	FLOAT *r_d3, FLOAT *r_pm, FLOAT *r_pp1, FLOAT *r_pp2, FLOAT *r_pp3)
{
	/* find roots of erdl_dp, erdl_pm, and erdl_pp */

	FLOAT rtm1,rt4m1;
	FLOAT d1,d2;
	FLOAT r1,r2;

	rtm1 = sqrt(m1);
	rt4m1 = sqrt(rtm1);
	d2 = pow(m1,1/3.0)+pow(1-m1,1/3.0);
	d2 = d2*d2*d2;
	d2 = sqrt(d2);
	r2 = pow(m1*d2,1/3.0);
	d1 = calc_d1(m1,LIMPREC);
	r1 = pow(m1*d1,1/3.0);
	if (d <= 0.0) {
	  *r_d1 = 0.0;
	  *r_d2 = 0.0;
	  *r_d3 = rt4m1;
	  *r_pm = 1.0;
	  *r_pp1 = 1.0;
	  *r_pp2 = 0.0;
	  *r_pp2 = 0.0;
	  return;
	}
	*r_pm = bisection3(erdl_pm,rtm1,1.0,d,m1,1e-16);
	*r_d1 = bisection3(erdl_dp,1e-30,rtm1*d,d,m1,1e-16);
	*r_pp1 = bisection3(erdl_pp,d,d+1,d,m1,1e-16);
	if (close_binary(d,m1,r1)) {
	  *r_d2 = bisection3(erdl_dp,rtm1*d,r1,d,m1,1e-16);
	  *r_d3 = bisection3(erdl_dp,r1,rt4m1,d,m1,1e-16);
	} else {
	  *r_d2 = 0.0;
	  *r_d3 = 0.0;
	}
	if (wide_binary(d,m1,r2)) {
	  *r_pp2 = bisection3(erdl_pp,(d-d2)+r2,d,d,m1,1e-16);
	  *r_pp3 = bisection3(erdl_pp,rtm1,r2,d,m1,1e-16);
	} else {
	  *r_pp2 = 0.0;
	  *r_pp3 = 0.0;
	} 	
}

static void x1x2 (FLOAT r, FLOAT cval, FLOAT d, FLOAT m1,
   FLOAT *x1, FLOAT *x2)
{
	FLOAT m2;
	FLOAT sval;

	m2 = 1.0-m1;
	if (m2 < 0.0) m2 = 0.0;
        if (cval < -1.0) {
          cval = -1.0;
          sval = 0.0;
        }
        else if (cval > 1.0) {
          cval = 1.0;
          sval = 0.0;
        }
	else 
          sval = sqrt(1.0-cval*cval);  
	*x1 = r*cval-m2*d;
	*x2 = r*sval;
}

static FLOAT perp_caust(FLOAT r, boolean pmflag, FLOAT d, FLOAT m1,
        FLOAT y1c, FLOAT y2c)
{         
        /* calculate (y(r)-y_c) \dot t_y */
          
        FLOAT cval,x1,x2,ty1,ty2;
        FLOAT y1,y2;
        FLOAT val;
          
        if (pmflag)
            cval = cosplus(r,d,m1);  
        else
            cval = cosminus(r,d,m1);      
        x1x2(r,cval,d,m1,&x1,&x2);
        binlens(x1,x2,d,m1,&y1,&y2);
        tangent_y(x1,x2,d,m1,&ty1,&ty2);
        val = (y1-y1c)*ty1+(y2-y2c)*ty2;
        return(val);
}


static FLOAT perp_caustd1(FLOAT r, boolean pmflag, FLOAT d, FLOAT m1,
        FLOAT y1c, FLOAT y2c)
{
        /* calculate t_y + (y(r)-y_c) \dot t'_y */

	FLOAT y1,y2;
	FLOAT y11,y12,y22,y111,y112,y122,y222,y1111,y1112,y1122,y1222,y2222;
        FLOAT cval,x1,x2,ty1,ty2;
	FLOAT tty1,tty2;
        FLOAT n1,n2;
	FLOAT n11,n12,n22;
        FLOAT val;
	FLOAT normsq;

        if (pmflag)
            cval = cosplus(r,d,m1);
        else
            cval = cosminus(r,d,m1);
        x1x2(r,cval,d,m1,&x1,&x2);
        normal(x1,x2,d,m1,&n1,&n2);
        binlens(x1,x2,d,m1,&y1,&y2);
        binlensd1(x1,x2,d,m1,&y11,&y12,&y22);
	binlensd2(x1,x2,d,m1,&y111,&y112,&y122,&y222);
	binlensd3(x1,x2,d,m1,&y1111,&y1112,&y1122,&y1222,&y2222);
	n1 = y111*y22+y122*y11-2.0*y112*y12;
	n2 = y112*y22+y222*y11-2.0*y122*y12;
	n11 = y1111*y22-2.0*y1112*y12+y1122*y11+2.0*y111*y122-2.0*y112*y112;
	n12 = y1112*y22-2.0*y1122*y12+y1222*y11+y111*y222-y112*y122;
	n22 = y1122*y22-2.0*y1222*y12+y2222*y11+2.0*y112*y222-2.0*y122*y122;
        ty1 = n1*y12-n2*y11;
        ty2 = n1*y22-n2*y12;
	tty1 = n1*(n12*y12-n22*y11+n1*y122-n2*y112)
		  -n2*(n11*y12-n12*y11+n1*y112-n2*y111);
	tty2 = n1*(n12*y22-n22*y12+n1*y222-n2*y122)
		  -n2*(n11*y22-n12*y12+n1*y122-n2*y112);
	normsq = SQ(ty1)+SQ(ty2);
        val = normsq+(y1-y1c)*tty1+(y2-y2c)*tty2;
        return(val);
}

static FLOAT perp_caustd2(FLOAT r, boolean pmflag, FLOAT d, FLOAT m1,
        FLOAT y1c, FLOAT y2c)
{
        /* calculate 2nd derivative of perp_caust */

	FLOAT y1,y2;
	FLOAT y11,y12,y22,y111,y112,y122,y222,y1111,y1112,y1122,y1222,y2222;
	FLOAT y11111,y11112,y11122,y11222,y12222,y22222;
        FLOAT cval,x1,x2,ty1,ty2;
	FLOAT tty1,tty2;
	FLOAT ttty1, ttty2;
        FLOAT n1,n2;
	FLOAT n11,n12,n22;
	FLOAT n111,n112,n122,n222;
        FLOAT val;
	FLOAT C11,C12,C21,C22,C111,C112,C121,C122,C211,C212,C221,C222;

        if (pmflag)
            cval = cosplus(r,d,m1);
        else
            cval = cosminus(r,d,m1);
        x1x2(r,cval,d,m1,&x1,&x2);
        normal(x1,x2,d,m1,&n1,&n2);
        binlens(x1,x2,d,m1,&y1,&y2);
        binlensd1(x1,x2,d,m1,&y11,&y12,&y22);
	binlensd2(x1,x2,d,m1,&y111,&y112,&y122,&y222);
	binlensd3(x1,x2,d,m1,&y1111,&y1112,&y1122,&y1222,&y2222);
	binlensd4(x1,x2,d,m1,&y11111,&y11112,&y11122,&y11222,&y12222,&y22222);
	n1 = y111*y22+y122*y11-2.0*y112*y12;
	n2 = y112*y22+y222*y11-2.0*y122*y12;
	n11 = y1111*y22-2.0*y1112*y12+y1122*y11+2.0*y111*y122-2.0*y112*y112;
	n12 = y1112*y22-2.0*y1122*y12+y1222*y11+y111*y222-y112*y122;
	n22 = y1122*y22-2.0*y1222*y12+y2222*y11+2.0*y112*y222-2.0*y122*y122;
        n111 = y11111*y22-2.0*y11112*y12+y11122*y11
                 +3.0*y1111*y122-6.0*y1112*y112+3.0*y1122*y111;
        n112 = y11112*y22-2.0*y11122*y12+y11222*y11
                 +y1111*y222-3.0*y1122*y112+2.0*y1222*y111;
        n122 = y11122*y22-2.0*y11222*y12+y12222*y11
                 +2.0*y1112*y222-3.0*y1122*y122+y2222*y111;
        n222 = y22222*y11-2.0*y12222*y12+y11222*y22
                 +3.0*y1122*y222-6.0*y1222*y122+3.0*y2222*y112;
        ty1 = n1*y12-n2*y11;
        ty2 = n1*y22-n2*y12;
	C11 = n12*y12-n22*y11+n1*y122-n2*y112;
	C12 = n11*y12-n12*y11+n1*y112-n2*y111;
	C21 = n12*y22-n22*y12+n1*y222-n2*y122;
	C22 = n11*y22-n12*y12+n1*y122-n2*y112;
	C111 = n112*y12-n122*y11+n11*y122-n12*y112
	         +n12*y112-n22*y111+n1*y1122-n2*y1112;
	C112 = n122*y12-n222*y11+n12*y122-n22*y112
	         +n12*y122-n22*y112+n1*y1222-n2*y1122;
	C121 = n111*y12-n112*y11+n11*y112-n12*y111
		 +n11*y112-n12*y111+n1*y1112-n2*y1111;
	C122 = n112*y12-n122*y11+n12*y112-n22*y111
		 +n11*y122-n12*y112+n1*y1122-n2*y1112;
	C211 = n112*y22-n122*y12+n11*y222-n12*y122
		 +n12*y122-n22*y112+n1*y1222-n2*y1122;
	C212 = n122*y22-n222*y12+n12*y222-n22*y122
		 +n12*y222-n22*y122+n1*y2222-n2*y1222;
	C221 = n111*y22-n112*y12+n11*y122-n12*y112
		 +n11*y122-n12*y112+n1*y1122-n2*y1112;
	C222 = n112*y22-n122*y12+n12*y122-n22*y112
	         +n11*y222-n12*y122+n1*y1222-n2*y1122;
	tty1 = n1*C11-n2*C12;
	tty2 = n1*C21-n2*C22;
	ttty1 = n1*(n12*C11-n22*C12+n1*C112-n2*C122)
	      -n2*(n11*C11-n12*C12+n1*C111-n2*C121);
	ttty2 = n1*(n12*C21-n22*C22+n1*C212-n2*C222)
	      -n2*(n11*C21-n12*C22+n1*C211-n2*C221);
        val = 3.0*ty1*tty1+3.0*ty2*tty2+(y1-y1c)*ttty1+(y2-y2c)*ttty2;
        return(val);
}

static FLOAT curvature_der(FLOAT r, boolean pmflag, FLOAT d, FLOAT m1)
{
	FLOAT y1,y2;
	FLOAT y11,y12,y22,y111,y112,y122,y222,y1111,y1112,y1122,y1222,y2222;
	FLOAT y11111,y11112,y11122,y11222,y12222,y22222;
        FLOAT cval,x1,x2,ty1,ty2;
	FLOAT tty1,tty2;
	FLOAT ttty1, ttty2;
        FLOAT n1,n2;
	FLOAT n11,n12,n22;
	FLOAT n111,n112,n122,n222;
        FLOAT val;
	FLOAT C11,C12,C21,C22,C111,C112,C121,C122,C211,C212,C221,C222;
	FLOAT normsq;
	FLOAT norm5;

        if (pmflag)
            cval = cosplus(r,d,m1);
        else
            cval = cosminus(r,d,m1);
        x1x2(r,cval,d,m1,&x1,&x2);
        normal(x1,x2,d,m1,&n1,&n2);
        binlens(x1,x2,d,m1,&y1,&y2);
        binlensd1(x1,x2,d,m1,&y11,&y12,&y22);
	binlensd2(x1,x2,d,m1,&y111,&y112,&y122,&y222);
	binlensd3(x1,x2,d,m1,&y1111,&y1112,&y1122,&y1222,&y2222);
	binlensd4(x1,x2,d,m1,&y11111,&y11112,&y11122,&y11222,&y12222,&y22222);
	n1 = y111*y22+y122*y11-2.0*y112*y12;
	n2 = y112*y22+y222*y11-2.0*y122*y12;
	n11 = y1111*y22-2.0*y1112*y12+y1122*y11+2.0*y111*y122-2.0*y112*y112;
	n12 = y1112*y22-2.0*y1122*y12+y1222*y11+y111*y222-y112*y122;
	n22 = y1122*y22-2.0*y1222*y12+y2222*y11+2.0*y112*y222-2.0*y122*y122;
        n111 = y11111*y22-2.0*y11112*y12+y11122*y11
                 +3.0*y1111*y122-6.0*y1112*y112+3.0*y1122*y111;
        n112 = y11112*y22-2.0*y11122*y12+y11222*y11
                 +y1111*y222-3.0*y1122*y112+2.0*y1222*y111;
        n122 = y11122*y22-2.0*y11222*y12+y12222*y11
                 +2.0*y1112*y222-3.0*y1122*y122+y2222*y111;
        n222 = y22222*y11-2.0*y12222*y12+y11222*y22
                 +3.0*y1122*y222-6.0*y1222*y122+3.0*y2222*y112;
        ty1 = n1*y12-n2*y11;
        ty2 = n1*y22-n2*y12;
	C11 = n12*y12-n22*y11+n1*y122-n2*y112;
	C12 = n11*y12-n12*y11+n1*y112-n2*y111;
	C21 = n12*y22-n22*y12+n1*y222-n2*y122;
	C22 = n11*y22-n12*y12+n1*y122-n2*y112;
	C111 = n112*y12-n122*y11+n11*y122-n12*y112
	         +n12*y112-n22*y111+n1*y1122-n2*y1112;
	C112 = n122*y12-n222*y11+n12*y122-n22*y112
	         +n12*y122-n22*y112+n1*y1222-n2*y1122;
	C121 = n111*y12-n112*y11+n11*y112-n12*y111
		 +n11*y112-n12*y111+n1*y1112-n2*y1111;
	C122 = n112*y12-n122*y11+n12*y112-n22*y111
		 +n11*y122-n12*y112+n1*y1122-n2*y1112;
	C211 = n112*y22-n122*y12+n11*y222-n12*y122
		 +n12*y122-n22*y112+n1*y1222-n2*y1122;
	C212 = n122*y22-n222*y12+n12*y222-n22*y122
		 +n12*y222-n22*y122+n1*y2222-n2*y1222;
	C221 = n111*y22-n112*y12+n11*y122-n12*y112
		 +n11*y122-n12*y112+n1*y1122-n2*y1112;
	C222 = n112*y22-n122*y12+n12*y122-n22*y112
	         +n11*y222-n12*y122+n1*y1222-n2*y1122;
	tty1 = n1*C11-n2*C12;
	tty2 = n1*C21-n2*C22;
	ttty1 = n1*(n12*C11-n22*C12+n1*C112-n2*C122)
	      -n2*(n11*C11-n12*C12+n1*C111-n2*C121);
	ttty2 = n1*(n12*C21-n22*C22+n1*C212-n2*C222)
	      -n2*(n11*C21-n12*C22+n1*C211-n2*C221);
	normsq = SQ(ty1)+SQ(ty2);
	norm5 = normsq*normsq*sqrt(normsq);
	val = ((ty1*ttty2-ty2*ttty1)*normsq
	  -3.0*(ty1*tty2-ty2*tty1)*(ty1*tty1+ty2*tty2))/norm5;
        return(val);
}



static void add_to_list(PTLIST_PTR *list, FLOAT r, boolean pmflag,
	FLOAT d, FLOAT m1, boolean invert)
{
	FLOAT cval;
	FLOAT x1,x2;

        if (pmflag)
          cval = cosplus(r,d,m1);
        else
          cval = cosminus(r,d,m1);
        x1x2(r,cval,d,m1,&x1,&x2);
	if (invert)
	  add_ptlist(list,x1,-x2);
	else
	  add_ptlist(list,x1,x2);
}

static void add_crit(FLOAT rmin, FLOAT rmax, FLOAT ty1min, FLOAT ty2min,
	FLOAT ty1max, FLOAT ty2max, boolean pmflag,
	FLOAT d, FLOAT m1, FLOAT y1c, FLOAT y2c, FLOAT rcirc,
	PTLIST_PTR *crit_points)
{
	FLOAT cval;
	FLOAT x1,x2,y1,y2;
	FLOAT minperp_p, maxperp_p;
	FLOAT minperp_m, maxperp_m;
	FLOAT mindistsq_p,mindistsq_m;
	FLOAT maxdistsq_p,maxdistsq_m;
	FLOAT rcircsq;
	FLOAT rperp;
	FLOAT minperpd_p,maxperpd_p;
	FLOAT minperpdd_p,maxperpdd_p;
	FLOAT minperpd_m,maxperpd_m;
	FLOAT minperpdd_m,maxperpdd_m;
	FLOAT perpd;
	FLOAT rddperp;
	FLOAT rdperp1,rdperp2;
	FLOAT perp1=0.,perp2=0.;
	

	/* symmetry with respect to x1 (y1)-axis */
	/* check both semi-planes */

	rcircsq = SQ(rcirc);
	if (pmflag)
	    cval = cosplus(rmin,d,m1);	  
	else
	    cval = cosminus(rmin,d,m1);	  
	x1x2(rmin,cval,d,m1,&x1,&x2);
	binlens(x1,x2,d,m1,&y1,&y2);
	mindistsq_p = SQ(y1-y1c)+SQ(y2-y2c);
	mindistsq_m = SQ(y1-y1c)+SQ(y2+y2c);
	minperp_p = (y1-y1c)*ty1min+(y2-y2c)*ty2min;
	minperp_m = (y1-y1c)*ty1min+(y2+y2c)*ty2min;
	if (pmflag)
	    cval = cosplus(rmax,d,m1);	  
	else
	    cval = cosminus(rmax,d,m1);	  
	x1x2(rmax,cval,d,m1,&x1,&x2);
	binlens(x1,x2,d,m1,&y1,&y2);
	maxdistsq_p = SQ(y1-y1c)+SQ(y2-y2c);
	maxdistsq_m = SQ(y1-y1c)+SQ(y2+y2c);
	maxperp_p = (y1-y1c)*ty1max+(y2-y2c)*ty2max;
	maxperp_m = (y1-y1c)*ty1max+(y2+y2c)*ty2max;
	minperpd_p = perp_caustd1(rmin,pmflag,d,m1,y1c,y2c);
	maxperpd_p = perp_caustd1(rmax,pmflag,d,m1,y1c,y2c);
	minperpdd_p = perp_caustd2(rmin,pmflag,d,m1,y1c,y2c);
	maxperpdd_p = perp_caustd2(rmax,pmflag,d,m1,y1c,y2c);
	minperpd_m = perp_caustd1(rmin,pmflag,d,m1,y1c,-y2c);
	maxperpd_m = perp_caustd1(rmax,pmflag,d,m1,y1c,-y2c);
	minperpdd_m = perp_caustd2(rmin,pmflag,d,m1,y1c,-y2c);
	maxperpdd_m = perp_caustd2(rmax,pmflag,d,m1,y1c,-y2c);
	rdperp1 = 0.0;
	rdperp2 = 0.0;
	if (minperpdd_p*maxperpdd_p < 0) {
	  rddperp = get_rootdefm(perp_caustd2,rmin,rmax,pmflag,
			d,m1,y1c,y2c,round(maxperpdd_p/fabs(maxperpdd_p)),RTPREC);
	  perpd = perp_caustd1(rddperp,pmflag,d,m1,y1c,y2c);
	  if (perpd == 0.0) {
	    rdperp1 = rddperp;
	    perp1 = perp_caust(rdperp1,pmflag,d,m1,y1c,y2c);
	  } else {
	    if (minperpd_p*perpd < 0) {
	      rdperp1 = get_rootdefm(perp_caustd1,rmin,rddperp,pmflag,
                        d,m1,y1c,y2c,round(perpd/fabs(perpd)),RTPREC);
	      perp1 = perp_caust(rdperp1,pmflag,d,m1,y1c,y2c);
	    }
	    if (perpd*maxperpd_p < 0) {
	      rdperp2 = get_rootdefm(perp_caustd1,rddperp,rmax,pmflag,
                        d,m1,y1c,y2c,-round(perpd/fabs(perpd)),RTPREC);
	      perp2 = perp_caust(rdperp2,pmflag,d,m1,y1c,y2c);
	    }
	  }
	} else if (minperpd_p*maxperpd_p < 0) {
	  rdperp1 = get_rootdefm(perp_caustd1,rmin,rmax,pmflag,
			d,m1,y1c,y2c,round(maxperpd_p/fabs(maxperpd_p)),RTPREC);
	  perp1 = perp_caust(rdperp1,pmflag,d,m1,y1c,y2c);
	}
	if (rdperp1 > 0.0) {
	  if (perp1 == 0.0) {
	    rperp = rdperp1;
	    add_to_list(crit_points,rperp,pmflag,d,m1,FALSE);
	  }
	  else {
	    if (minperp_p*perp1 < 0) {
	      rperp = get_rootdefm(perp_caust,rmin,rdperp1,pmflag,
		d,m1,y1c,y2c,round(perp1/fabs(perp1)),RTPREC); 
	      add_to_list(crit_points,rperp,pmflag,d,m1,FALSE);
	    }
	    if (rdperp2 == 0.0) {
	      if (maxperp_p*perp1 < 0) {
	        rperp = get_rootdefm(perp_caust,rdperp1,rmax,pmflag,
	 	  d,m1,y1c,y2c,-round(perp1/fabs(perp1)),RTPREC);
	        add_to_list(crit_points,rperp,pmflag,d,m1,FALSE);
	      }
	    }
	  }
	}
	if (rdperp2 > 0.0) {
	  if (perp2 == 0.0) {
	    rperp = rdperp2;
	    add_to_list(crit_points,rperp,pmflag,d,m1,FALSE);
	  } 
	  else {
	    if (maxperp_p*perp2 < 0) {
	      rperp = get_rootdefm(perp_caust,rdperp2,rmax,pmflag,
	 	d,m1,y1c,y2c,-round(perp2/fabs(perp2)),RTPREC); 
	      add_to_list(crit_points,rperp,pmflag,d,m1,FALSE);
	    }
	    if (rdperp1 == 0.0) {
	      if (minperp_p*perp2 < 0) {
	        rperp = get_rootdefm(perp_caust,rmin,rdperp2,pmflag,
	 	  d,m1,y1c,y2c,round(perp2/fabs(perp2)),RTPREC);
	        add_to_list(crit_points,rperp,pmflag,d,m1,FALSE);
	      }
	    }
	  }
	}
	if (rdperp1 > 0.0 && rdperp2 > 0.0) {
	  if (perp1*perp2 < 0) {
	    rperp = get_rootdefm(perp_caust,rdperp1,rdperp2,pmflag,
	 	d,m1,y1c,y2c,round(perp2/fabs(perp2)),RTPREC); 
	    add_to_list(crit_points,rperp,pmflag,d,m1,FALSE);
	  }
	}
	if (rdperp1 == 0.0 && rdperp2 == 0.0) {
	  if (minperp_p*maxperp_p < 0) {
	    rperp = get_rootdefm(perp_caust,rmin,rmax,pmflag,
		 	d,m1,y1c,y2c,round(maxperp_p/fabs(maxperp_p)),RTPREC); 
	    add_to_list(crit_points,rperp,pmflag,d,m1,FALSE);
	  }
	}
	rdperp1 = 0.0;
	rdperp2 = 0.0;
	if (minperpdd_m*maxperpdd_m < 0) {
	  rddperp = get_rootdefm(perp_caustd2,rmin,rmax,pmflag,
			d,m1,y1c,-y2c,round(maxperpdd_m/fabs(maxperpdd_m)),RTPREC);
	  perpd = perp_caustd1(rddperp,pmflag,d,m1,y1c,-y2c);
	  if (perpd == 0.0) {
	    rdperp1 = rddperp;
	    perp1 = perp_caust(rdperp1,pmflag,d,m1,y1c,-y2c);
	  } else {
	    if (minperpd_m*perpd < 0) {
	      rdperp1 = get_rootdefm(perp_caustd1,rmin,rddperp,pmflag,
                        d,m1,y1c,-y2c,round(perpd/fabs(perpd)),RTPREC);
	      perp1 = perp_caust(rdperp1,pmflag,d,m1,y1c,-y2c);
	    }
	    if (perpd*maxperpd_m < 0) {
	      rdperp2 = get_rootdefm(perp_caustd1,rddperp,rmax,pmflag,
                        d,m1,y1c,-y2c,-round(perpd/fabs(perpd)),RTPREC);
	      perp2 = perp_caust(rdperp2,pmflag,d,m1,y1c,-y2c);
	    }
	  }
	} else if (minperpd_m*maxperpd_m < 0) {
	  rdperp1 = get_rootdefm(perp_caustd1,rmin,rmax,pmflag,
			d,m1,y1c,-y2c,round(maxperpd_m/fabs(maxperpd_m)),RTPREC);
	  perp1 = perp_caust(rdperp1,pmflag,d,m1,y1c,-y2c);
	}
	if (rdperp1 > 0.0) {
	  if (perp1 == 0.0) {
	    rperp = rdperp1;
	    add_to_list(crit_points,rperp,pmflag,d,m1,TRUE);
	  }
	  else {
	    if (minperp_m*perp1 < 0) {
	      rperp = get_rootdefm(perp_caust,rmin,rdperp1,pmflag,
		d,m1,y1c,-y2c,round(perp1/fabs(perp1)),RTPREC); 
	      add_to_list(crit_points,rperp,pmflag,d,m1,TRUE);
	    }
	    if (rdperp2 == 0.0) {
	      if (maxperp_m*perp1 < 0) {
	        rperp = get_rootdefm(perp_caust,rdperp1,rmax,pmflag,
	 	  d,m1,y1c,-y2c,-round(perp1/fabs(perp1)),RTPREC);
	        add_to_list(crit_points,rperp,pmflag,d,m1,TRUE);
	      }
	    }
	  }
	}
	if (rdperp2 > 0.0) {
	  if (perp2 == 0.0) {
	    rperp = rdperp2;
	    add_to_list(crit_points,rperp,pmflag,d,m1,TRUE);
	  } 
	  else {
	    if (maxperp_m*perp2 < 0) {
	      rperp = get_rootdefm(perp_caust,rdperp2,rmax,pmflag,
	 	d,m1,y1c,-y2c,-round(perp2/fabs(perp2)),RTPREC); 
	      add_to_list(crit_points,rperp,pmflag,d,m1,TRUE);
	    }
	    if (rdperp1 == 0.0) {
	      if (minperp_m*perp2 < 0) {
	        rperp = get_rootdefm(perp_caust,rmin,rdperp2,pmflag,
	 	  d,m1,y1c,-y2c,round(perp2/fabs(perp2)),RTPREC);
	        add_to_list(crit_points,rperp,pmflag,d,m1,TRUE);
	      }
	    }
	  }
	}
	if (rdperp1 > 0.0 && rdperp2 > 0.0) {
	  if (perp1*perp2 < 0) {
	    rperp = get_rootdefm(perp_caust,rdperp1,rdperp2,pmflag,
	 	d,m1,y1c,-y2c,round(perp2/fabs(perp2)),RTPREC); 
	    add_to_list(crit_points,rperp,pmflag,d,m1,TRUE);
	  }
	}
	if (rdperp1 == 0.0 && rdperp2 == 0.0) {
	  if (minperp_m*maxperp_m < 0) {
	    rperp = get_rootdefm(perp_caust,rmin,rmax,pmflag,
		 	d,m1,y1c,-y2c,round(maxperp_m/fabs(maxperp_m)),RTPREC); 
	    add_to_list(crit_points,rperp,pmflag,d,m1,TRUE);
	  }
	}
}

void get_crit_circ(FLOAT d, FLOAT m1, FLOAT y1c, FLOAT y2c, FLOAT rcirc,
	PTLIST_PTR *crit_points)
{
	FLOAT r_d1,r_d2,r_d3,r_pm,r_pp1,r_pp2,r_pp3;
	FLOAT cval;
	FLOAT rmin,rmax;
	FLOAT d1,r1;
	FLOAT d2,r2;
	FLOAT x1,x2;
	FLOAT rs,rb,ra,rcusp_l,rcusp_u;
	FLOAT rcusp2_l,rcusp2_u;
	FLOAT rcmin;
	FLOAT ty1,ty2,ty1min,ty2min;
	FLOAT ty1c,ty2c;
	FLOAT ty21,ty22;
	FLOAT ty1max,ty2max;
	FLOAT curv_min;

        r_pp3 = 0.0;
        erdl_limits(d,m1,&r_d1,&r_d2,&r_d3,&r_pm,&r_pp1,&r_pp2,&r_pp3);
	d1 = calc_d1(m1,LIMPREC);
	r1 = pow(m1*d1,1/3.0);
        d2 = pow(m1,1/3.0)+pow(1-m1,1/3.0);
        d2 = d2*d2*d2;
        d2 = sqrt(d2);
        r2 = pow(m1*d2,1/3.0);
	if (close_binary(d,m1,r1)) {
	  rmin = r_d1;
	  rmax = r_d2;
          rb = get_rootdef(offdiag,rmin,rmax,TRUE,d,m1,-1,RTPREC);
	  get_rootdef_lu(cspf,rb,rmax,TRUE,d,m1,-1,&rcusp_l,&rcusp_u,CUSP_PREC); 
	  /* intervals: ]rmin,rcusp_l[ and ]rcusp_u,rmax[ */
	  cval = get_cos(rmin,TRUE,d,m1);
	  x1x2(rmin,cval,d,m1,&x1,&x2);
	  tangent_y(x1,x2,d,m1,&ty1min,&ty2min);
	  cval = get_cos(rcusp_l,TRUE,d,m1);
	  x1x2(rcusp_l,cval,d,m1,&x1,&x2);
	  tangent_y(x1,x2,d,m1,&ty1,&ty2);
	  curv_min = curvature_der(rmin,TRUE,d,m1);
	  if (curv_min > 0.0) {
	    rcmin = get_rootdef(curvature_der,rmin,rcusp_l,TRUE,d,m1,-1,RTPREC);
	    if (perp_caust(rcmin,TRUE,d,m1,y1c,y2c) == 0.0) 
	      add_to_list(crit_points,rcmin,TRUE,d,m1,FALSE);
	    if (perp_caust(rcmin,TRUE,d,m1,y1c,-y2c) == 0.0) 
	      add_to_list(crit_points,rcmin,TRUE,d,m1,TRUE);
	    cval = get_cos(rcmin,TRUE,d,m1);
	    x1x2(rcmin,cval,d,m1,&x1,&x2);
	    tangent_y(x1,x2,d,m1,&ty1c,&ty2c);
	    add_crit(rmin,rcmin,ty1min,ty2min,ty1c,ty2c,TRUE,d,m1,
		y1c,y2c,rcirc,crit_points);
	    add_crit(rcmin,rcusp_l,ty1c,ty2c,ty1,ty2,TRUE,d,m1,
		y1c,y2c,rcirc,crit_points);
	  } else { 
	    add_crit(rmin,rcusp_l,ty1min,ty2min,ty1,ty2,TRUE,d,m1,
		y1c,y2c,rcirc,crit_points);
	  }
	  cval = get_cos(rcusp_u,TRUE,d,m1);
	  x1x2(rcusp_u,cval,d,m1,&x1,&x2);
	  tangent_y(x1,x2,d,m1,&ty1,&ty2);
	  cval = get_cos(rmax,TRUE,d,m1);
	  x1x2(rmax,cval,d,m1,&x1,&x2);
	  tangent_y(x1,x2,d,m1,&ty1max,&ty2max);
	  curv_min = curvature_der(rmax,TRUE,d,m1);
	  if (curv_min > 0.0) {
	    rcmin = get_rootdef(curvature_der,rcusp_u,rmax,TRUE,d,m1,-1,RTPREC);
	    if (perp_caust(rcmin,TRUE,d,m1,y1c,y2c) == 0.0) 
	      add_to_list(crit_points,rcmin,TRUE,d,m1,FALSE);
	    if (perp_caust(rcmin,TRUE,d,m1,y1c,-y2c) == 0.0) 
	      add_to_list(crit_points,rcmin,TRUE,d,m1,TRUE);
	    cval = get_cos(rcmin,TRUE,d,m1);
	    x1x2(rcmin,cval,d,m1,&x1,&x2);
	    tangent_y(x1,x2,d,m1,&ty1c,&ty2c);
	    add_crit(rcusp_u,rcmin,ty1,ty2,ty1c,ty2c,TRUE,d,m1,
		y1c,y2c,rcirc,crit_points);
	    add_crit(rcmin,rmax,ty1c,ty2c,ty1max,ty2max,TRUE,d,m1,
		y1c,y2c,rcirc,crit_points);
	  } else {
	    add_crit(rcusp_u,rmax,ty1,ty2,ty1max,ty2max,TRUE,d,m1,
		y1c,y2c,rcirc,crit_points);
	  }
	  rmin = r_d1;
	  rmax = r_d2;
          ra = get_rootdef(offdiag,rmin,rmax,FALSE,d,m1,-1,RTPREC);
          get_rootdef_lu(cspf,rmin,ra,FALSE,d,m1,1,&rcusp_l,&rcusp_u,CUSP_PREC);
          get_rootdef_lu(cspf,ra,rmax,FALSE,d,m1,-1,&rcusp2_l,&rcusp2_u,CUSP_PREC);
	  /* intervals:  ]rmin,rcusp_l[, ]rcusp_u,rcusp2_l[, ]rcusp2_u,rmax[ */
	  cval = get_cos(rmin,FALSE,d,m1);
	  x1x2(rmin,cval,d,m1,&x1,&x2);
	  tangent_y(x1,x2,d,m1,&ty1min,&ty2min);
	  cval = get_cos(rcusp_l,FALSE,d,m1);
	  x1x2(rcusp_l,cval,d,m1,&x1,&x2);
	  tangent_y(x1,x2,d,m1,&ty1,&ty2);
	  curv_min = curvature_der(rmin,FALSE,d,m1);
	  if (curv_min > 0.0) {
	    rcmin = get_rootdef(curvature_der,rmin,rcusp_l,FALSE,d,m1,-1,RTPREC);
	    if (perp_caust(rcmin,FALSE,d,m1,y1c,y2c) == 0.0) 
	      add_to_list(crit_points,rcmin,FALSE,d,m1,FALSE);
	    if (perp_caust(rcmin,FALSE,d,m1,y1c,-y2c) == 0.0) 
	      add_to_list(crit_points,rcmin,FALSE,d,m1,TRUE);
	    cval = get_cos(rcmin,FALSE,d,m1);
	    x1x2(rcmin,cval,d,m1,&x1,&x2);
	    tangent_y(x1,x2,d,m1,&ty1c,&ty2c);
	    add_crit(rmin,rcmin,ty1min,ty2min,ty1c,ty2c,FALSE,d,m1,
		y1c,y2c,rcirc,crit_points);
	    add_crit(rcmin,rcusp_l,ty1c,ty2c,ty1,ty2,FALSE,d,m1,
		y1c,y2c,rcirc,crit_points);
          } else {	
	    add_crit(rmin,rcusp_l,ty1min,ty2min,ty1,ty2,FALSE,d,m1,
		y1c,y2c,rcirc,crit_points);
	  }
	  cval = get_cos(rcusp2_u,FALSE,d,m1);
	  x1x2(rcusp2_u,cval,d,m1,&x1,&x2);
	  tangent_y(x1,x2,d,m1,&ty1,&ty2);
	  cval = get_cos(rmax,FALSE,d,m1);
	  x1x2(rmin,cval,d,m1,&x1,&x2);
	  tangent_y(x1,x2,d,m1,&ty1max,&ty2max);
	  curv_min = curvature_der(rmax,FALSE,d,m1);
	  if (curv_min > 0.0) {
	    rcmin = get_rootdef(curvature_der,rcusp2_u,rmax,FALSE,d,m1,-1,RTPREC);
	    if (perp_caust(rcmin,FALSE,d,m1,y1c,y2c) == 0.0) 
	      add_to_list(crit_points,rcmin,FALSE,d,m1,FALSE);
	    if (perp_caust(rcmin,FALSE,d,m1,y1c,-y2c) == 0.0) 
	      add_to_list(crit_points,rcmin,FALSE,d,m1,TRUE);
	    cval = get_cos(rcmin,FALSE,d,m1);
	    x1x2(rcmin,cval,d,m1,&x1,&x2);
	    tangent_y(x1,x2,d,m1,&ty1c,&ty2c);
	    add_crit(rcusp2_u,rcmin,ty1,ty2,ty1c,ty2c,FALSE,d,m1,
		y1c,y2c,rcirc,crit_points);
	    add_crit(rcmin,rmax,ty1c,ty2c,ty1max,ty2max,FALSE,d,m1,
		y1c,y2c,rcirc,crit_points);
	  } else {
	    add_crit(rcusp2_u,rmax,ty1,ty2,ty1max,ty2max,FALSE,d,m1,
		y1c,y2c,rcirc,crit_points);
	  }
	  /* there must be a absolute curvature minimum */
	  rcmin = get_rootdef(curvature_der,rcusp_u,rcusp2_l,FALSE,d,m1,-1,RTPREC);
	  if (perp_caust(rcmin,FALSE,d,m1,y1c,y2c) == 0.0) 
	    add_to_list(crit_points,rcmin,FALSE,d,m1,FALSE);
	  if (perp_caust(rcmin,FALSE,d,m1,y1c,-y2c) == 0.0) 
	    add_to_list(crit_points,rcmin,FALSE,d,m1,TRUE);
	  cval = get_cos(rcmin,FALSE,d,m1);
	  x1x2(rcmin,cval,d,m1,&x1,&x2);
	  tangent_y(x1,x2,d,m1,&ty1c,&ty2c);
	  cval = get_cos(rcusp_u,FALSE,d,m1);
	  x1x2(rcusp_u,cval,d,m1,&x1,&x2);
	  tangent_y(x1,x2,d,m1,&ty1,&ty2);
	  cval = get_cos(rcusp2_l,FALSE,d,m1);
	  x1x2(rcusp2_l,cval,d,m1,&x1,&x2);
	  tangent_y(x1,x2,d,m1,&ty21,&ty22);
	  add_crit(rcusp_u,rcmin,ty1,ty2,ty1c,ty2c,FALSE,d,m1,
		y1c,y2c,rcirc,crit_points);
	  add_crit(rcmin,rcusp2_l,ty1c,ty2c,ty21,ty22,FALSE,d,m1,
		y1c,y2c,rcirc,crit_points);
	  rmin = r_pm;
	  if (fabs(cosplus((r_pm+r_d3)/2.0,d,m1)) < 1.0) 
            rmin = r_d3;
	  rmax = r_pp1;
          ra = get_rootdef(offdiag,rmin,rmax,TRUE,d,m1,1,RTPREC);
	  get_rootdef_lu(cspf,rmin,ra,TRUE,d,m1,-1,&rcusp_l,&rcusp_u,CUSP_PREC);
	  /* intervals: ]rmin,rcusp_l[ and ]rcusp_u,rmax[ */
	  cval = get_cos(rmin,TRUE,d,m1);
	  x1x2(rmin,cval,d,m1,&x1,&x2);
	  tangent_y(x1,x2,d,m1,&ty1min,&ty2min);
	  cval = get_cos(rcusp_l,TRUE,d,m1);
	  x1x2(rcusp_l,cval,d,m1,&x1,&x2);
	  tangent_y(x1,x2,d,m1,&ty1,&ty2);
	  curv_min = curvature_der(rmin,TRUE,d,m1);
          /* 10-11-07: remove condition */
	  rcmin = get_rootdef(curvature_der,rmin,rcusp_l,TRUE,d,m1,-1,RTPREC);
	  if (perp_caust(rcmin,TRUE,d,m1,y1c,y2c) == 0.0) 
	    add_to_list(crit_points,rcmin,TRUE,d,m1,FALSE);
	  if (perp_caust(rcmin,TRUE,d,m1,y1c,-y2c) == 0.0) 
	    add_to_list(crit_points,rcmin,TRUE,d,m1,TRUE);
	  cval = get_cos(rcmin,TRUE,d,m1);
	  x1x2(rcmin,cval,d,m1,&x1,&x2);
	  tangent_y(x1,x2,d,m1,&ty1c,&ty2c);
	  add_crit(rmin,rcmin,ty1min,ty2min,ty1c,ty2c,TRUE,d,m1,
	    y1c,y2c,rcirc,crit_points);
	  add_crit(rcmin,rcusp_l,ty1c,ty2c,ty1,ty2,TRUE,d,m1,
	    y1c,y2c,rcirc,crit_points);
	  rcmin = get_rootdef(curvature_der,rcusp_u,rmax,TRUE,d,m1,-1,RTPREC);
	  if (perp_caust(rcmin,TRUE,d,m1,y1c,y2c) == 0.0) 
	    add_to_list(crit_points,rcmin,TRUE,d,m1,FALSE);
	  if (perp_caust(rcmin,TRUE,d,m1,y1c,-y2c) == 0.0) 
	    add_to_list(crit_points,rcmin,TRUE,d,m1,TRUE);
	  cval = get_cos(rcmin,TRUE,d,m1);
	  x1x2(rcmin,cval,d,m1,&x1,&x2);
	  tangent_y(x1,x2,d,m1,&ty1c,&ty2c);
	  cval = get_cos(rcusp_u,TRUE,d,m1);
	  x1x2(rcusp_u,cval,d,m1,&x1,&x2);
	  tangent_y(x1,x2,d,m1,&ty1,&ty2);
	  /* at rmax, tangent vector is (-1.0,0.0) */
	  add_crit(rcusp_u,rcmin,ty1,ty2,ty1c,ty2c,TRUE,d,m1,
		y1c,y2c,rcirc,crit_points);
	  add_crit(rcmin,rmax,ty1c,ty2c,-1.0,0.0,TRUE,d,m1,
		y1c,y2c,rcirc,crit_points);
	  if (fabs(cosplus((r_pm+r_d3)/2.0,d,m1)) < 1.0) {
	    rmin = r_d3;
	    rmax = r_pm;
	    cval = get_cos(rmin,FALSE,d,m1);
	    x1x2(rmin,cval,d,m1,&x1,&x2);
	    tangent_y(x1,x2,d,m1,&ty1min,&ty2min);
	    curv_min = curvature_der(rmin,FALSE,d,m1);
	    if (curv_min > 0.0) {
	      rcmin = get_rootdef(curvature_der,rmin,rmax,FALSE,d,m1,-1,RTPREC);
	      if (perp_caust(rcmin,FALSE,d,m1,y1c,y2c) == 0.0) 
	        add_to_list(crit_points,rcmin,FALSE,d,m1,FALSE);
	      if (perp_caust(rcmin,FALSE,d,m1,y1c,-y2c) == 0.0) 
	        add_to_list(crit_points,rcmin,FALSE,d,m1,TRUE);
	      cval = get_cos(rcmin,FALSE,d,m1);
	      x1x2(rcmin,cval,d,m1,&x1,&x2);
	      tangent_y(x1,x2,d,m1,&ty1c,&ty2c);
	      add_crit(rmin,rcmin,ty1min,ty2min,ty1c,ty2c,FALSE,d,m1,
		y1c,y2c,rcirc,crit_points);
	      add_crit(rcmin,rmax,ty1c,ty2c,-1.0,0.0,FALSE,d,m1,
		y1c,y2c,rcirc,crit_points);
	    } else {
	      /* at rmax, tangent vector is (-1.0,0.0) */
	      add_crit(rmin,rmax,ty1min,ty2max,-1.0,0.0,FALSE,d,m1,
		y1c,y2c,rcirc,crit_points);
	    }
	  }
	} 
	else if (wide_binary(d,m1,r2)) {
	  rmin = r_pp2;
	  rmax = r_pp1;
          ra = get_rootdef(offdiag,rmin,rmax,TRUE,d,m1,1,RTPREC);
	  get_rootdef_lu(cspf,rmin,ra,TRUE,d,m1,-1,&rcusp_l,&rcusp_u,CUSP_PREC);
	  /* intervals: ]rmin,rcusp_l[ and ]rcusp_u,rmax[ */
	  /* there must be a absolute curvature minimum */
	  rcmin = get_rootdef(curvature_der,rmin,rcusp_l,TRUE,d,m1,-1,RTPREC);
	  if (perp_caust(rcmin,TRUE,d,m1,y1c,y2c) == 0.0) 
	    add_to_list(crit_points,rcmin,TRUE,d,m1,FALSE);
	  if (perp_caust(rcmin,TRUE,d,m1,y1c,-y2c) == 0.0) 
	    add_to_list(crit_points,rcmin,TRUE,d,m1,TRUE);
	  cval = get_cos(rcmin,TRUE,d,m1);
	  x1x2(rcmin,cval,d,m1,&x1,&x2);
	  tangent_y(x1,x2,d,m1,&ty1c,&ty2c);
	  cval = get_cos(rcusp_l,TRUE,d,m1);
	  x1x2(rcusp_l,cval,d,m1,&x1,&x2);
	  tangent_y(x1,x2,d,m1,&ty1,&ty2);
	  /* at rmin, tangent vector is (-1.0,0.0) */
	  add_crit(rmin,rcmin,-1.0,0.0,ty1c,ty2c,TRUE,d,m1,
		y1c,y2c,rcirc,crit_points);
	  add_crit(rcmin,rcusp_l,ty1c,ty2c,ty1,ty2,TRUE,d,m1,
		y1c,y2c,rcirc,crit_points);
	  /* there must be a absolute curvature minimum */
	  rcmin = get_rootdef(curvature_der,rcusp_u,rmax,TRUE,d,m1,-1,RTPREC);
	  if (perp_caust(rcmin,TRUE,d,m1,y1c,y2c) == 0.0) 
	    add_to_list(crit_points,rcmin,TRUE,d,m1,FALSE);
	  if (perp_caust(rcmin,TRUE,d,m1,y1c,-y2c) == 0.0) 
	    add_to_list(crit_points,rcmin,TRUE,d,m1,TRUE);
	  cval = get_cos(rcmin,TRUE,d,m1);
	  x1x2(rcmin,cval,d,m1,&x1,&x2);
	  tangent_y(x1,x2,d,m1,&ty1c,&ty2c);
	  cval = get_cos(rcusp_u,TRUE,d,m1);
	  x1x2(rcusp_u,cval,d,m1,&x1,&x2);
	  tangent_y(x1,x2,d,m1,&ty1,&ty2);
	  add_crit(rcusp_u,rcmin,ty1,ty2,ty1c,ty2c,TRUE,d,m1,
		y1c,y2c,rcirc,crit_points);
	  /* at rmax, tangent vector is (-1.0,0.0) */
	  add_crit(rcmin,rmax,ty1c,ty2c,-1.0,0.0,TRUE,d,m1,
		y1c,y2c,rcirc,crit_points);
	  rmin = r_d1;
	  rmax = r_pp3;
	  /* interval ]rmin,rmax[ */
	  cval = get_cos(rmin,TRUE,d,m1);
	  x1x2(rmin,cval,d,m1,&x1,&x2);
	  tangent_y(x1,x2,d,m1,&ty1min,&ty2min);
	  curv_min = curvature_der(rmin,TRUE,d,m1);
	  if (curv_min > 0.0) {
	    rcmin = get_rootdef(curvature_der,rmin,rmax,TRUE,d,m1,-1,RTPREC);
	    if (perp_caust(rcmin,TRUE,d,m1,y1c,y2c) == 0.0) 
	      add_to_list(crit_points,rcmin,TRUE,d,m1,FALSE);
	    if (perp_caust(rcmin,TRUE,d,m1,y1c,-y2c) == 0.0) 
	      add_to_list(crit_points,rcmin,TRUE,d,m1,TRUE);
	    cval = get_cos(rcmin,TRUE,d,m1);
	    x1x2(rcmin,cval,d,m1,&x1,&x2);
	    tangent_y(x1,x2,d,m1,&ty1c,&ty2c);
	    add_crit(rmin,rcmin,ty1min,ty2min,ty1c,ty2c,TRUE,d,m1,
		y1c,y2c,rcirc,crit_points);
	    add_crit(rcmin,rmax,ty1c,ty2c,-1.0,0.0,TRUE,d,m1,
		y1c,y2c,rcirc,crit_points);
	  } else {
	    add_crit(rmin,rmax,ty1min,ty2min,-1.0,0.0,TRUE,d,m1,
		y1c,y2c,rcirc,crit_points);
	  }
	  rmin = r_d1;
	  rmax = r_pm;
          ra = get_rootdef(offdiag,rmin,rmax,FALSE,d,m1,-1,RTPREC);
	  get_rootdef_lu(cspf,rmin,ra,FALSE,d,m1,1,&rcusp_l,&rcusp_u,CUSP_PREC);
	  /* intervals: ]rmin,rcusp_l[ and ]rcusp_u,rmax[ */
	  cval = get_cos(rmin,FALSE,d,m1);
	  x1x2(rmin,cval,d,m1,&x1,&x2);
	  tangent_y(x1,x2,d,m1,&ty1min,&ty2min);
	  cval = get_cos(rcusp_l,FALSE,d,m1);
	  x1x2(rcusp_l,cval,d,m1,&x1,&x2);
	  tangent_y(x1,x2,d,m1,&ty1,&ty2);
	  curv_min = curvature_der(rmin,FALSE,d,m1);
	  if (curv_min > 0.0) {
	    rcmin = get_rootdef(curvature_der,rmin,rcusp_l,FALSE,d,m1,-1,RTPREC);
	    if (perp_caust(rcmin,FALSE,d,m1,y1c,y2c) == 0.0) 
	      add_to_list(crit_points,rcmin,FALSE,d,m1,FALSE);
	    if (perp_caust(rcmin,FALSE,d,m1,y1c,-y2c) == 0.0) 
	      add_to_list(crit_points,rcmin,FALSE,d,m1,TRUE);
	    cval = get_cos(rcmin,FALSE,d,m1);
	    x1x2(rcmin,cval,d,m1,&x1,&x2);
	    tangent_y(x1,x2,d,m1,&ty1c,&ty2c);
	    add_crit(rmin,rcmin,ty1min,ty2min,ty1c,ty2c,FALSE,d,m1,
		y1c,y2c,rcirc,crit_points);
	    add_crit(rcmin,rcusp_l,ty1c,ty2c,ty1,ty2,FALSE,d,m1,
		y1c,y2c,rcirc,crit_points);
	  } else {
	    add_crit(rmin,rcusp_l,ty1min,ty2min,ty1,ty2,FALSE,d,m1,
		y1c,y2c,rcirc,crit_points);
	  }
	  /* there must be a absolute curvature minimum */
	  rcmin = get_rootdef(curvature_der,rcusp_u,rmax,FALSE,d,m1,-1,RTPREC);
	  if (perp_caust(rcmin,FALSE,d,m1,y1c,y2c) == 0.0) 
	    add_to_list(crit_points,rcmin,FALSE,d,m1,FALSE);
	  if (perp_caust(rcmin,FALSE,d,m1,y1c,-y2c) == 0.0) 
	    add_to_list(crit_points,rcmin,FALSE,d,m1,TRUE);
	  cval = get_cos(rcmin,FALSE,d,m1);
	  x1x2(rcmin,cval,d,m1,&x1,&x2);
	  tangent_y(x1,x2,d,m1,&ty1c,&ty2c);
	  cval = get_cos(rcusp_u,FALSE,d,m1);
	  x1x2(rcusp_u,cval,d,m1,&x1,&x2);
	  tangent_y(x1,x2,d,m1,&ty1,&ty2);
	  add_crit(rcusp_u,rcmin,ty1,ty2,ty1c,ty2c,FALSE,d,m1,
		y1c,y2c,rcirc,crit_points);
	  /* at rmax, tangent vector is (-1.0,0.0) */
	  add_crit(rcmin,rmax,ty1c,ty2c,-1.0,0.0,FALSE,d,m1,
		y1c,y2c,rcirc,crit_points);
	}
	else {   
	  rmin = r_d1;
	  rmax = r_pp1;
	  rs = grid_brack(offdiag,rmin,rmax,TRUE,d,m1,-1,RTPREC);
	  ra = get_rootdef(offdiag,rs,rmax,TRUE,d,m1,1,RTPREC);
	  rb = get_rootdef(offdiag,rmin,rs,TRUE,d,m1,-1,RTPREC);
	  get_rootdef_lu(cspf,rb,ra,TRUE,d,m1,-1,&rcusp_l,&rcusp_u,CUSP_PREC);
	  /* intervals: ]rmin,rcusp_l[ and ]rcusp_u,rmax[ */
	  cval = get_cos(rmin,TRUE,d,m1);
	  x1x2(rmin,cval,d,m1,&x1,&x2);
	  tangent_y(x1,x2,d,m1,&ty1min,&ty2min);
	  cval = get_cos(rcusp_l,TRUE,d,m1);
	  x1x2(rcusp_l,cval,d,m1,&x1,&x2);
	  tangent_y(x1,x2,d,m1,&ty1,&ty2);
	  curv_min = curvature_der(rmin,TRUE,d,m1);
	  if (curv_min > 0.0) {
	    rcmin = get_rootdef(curvature_der,rmin,rcusp_l,TRUE,d,m1,-1,RTPREC);
	    if (perp_caust(rcmin,TRUE,d,m1,y1c,y2c) == 0.0) 
	      add_to_list(crit_points,rcmin,TRUE,d,m1,FALSE);
	    if (perp_caust(rcmin,TRUE,d,m1,y1c,-y2c) == 0.0) 
	      add_to_list(crit_points,rcmin,TRUE,d,m1,TRUE);
	    cval = get_cos(rcmin,TRUE,d,m1);
	    x1x2(rcmin,cval,d,m1,&x1,&x2);
	    tangent_y(x1,x2,d,m1,&ty1c,&ty2c);
	    add_crit(rmin,rcmin,ty1min,ty2min,ty1c,ty2c,TRUE,d,m1,
		y1c,y2c,rcirc,crit_points);
	    add_crit(rcmin,rcusp_l,ty1c,ty2c,ty1,ty2,TRUE,d,m1,
		y1c,y2c,rcirc,crit_points);
	  } else {
	    add_crit(rmin,rcusp_l,ty1min,ty2min,ty1,ty2,TRUE,d,m1,
		y1c,y2c,rcirc,crit_points);
	  }
	  /* there must be a absolute curvature minimum */
	  rcmin = get_rootdef(curvature_der,rcusp_u,rmax,TRUE,d,m1,-1,RTPREC);
	  if (perp_caust(rcmin,TRUE,d,m1,y1c,y2c) == 0.0) 
	    add_to_list(crit_points,rcmin,TRUE,d,m1,FALSE);
	  if (perp_caust(rcmin,TRUE,d,m1,y1c,-y2c) == 0.0) 
	    add_to_list(crit_points,rcmin,TRUE,d,m1,TRUE);
	  cval = get_cos(rcmin,TRUE,d,m1);
	  x1x2(rcmin,cval,d,m1,&x1,&x2);
	  tangent_y(x1,x2,d,m1,&ty1c,&ty2c);
	  cval = get_cos(rcusp_u,TRUE,d,m1);
	  x1x2(rcusp_u,cval,d,m1,&x1,&x2);
	  tangent_y(x1,x2,d,m1,&ty1,&ty2);
	  add_crit(rcusp_u,rcmin,ty1,ty2,ty1c,ty2c,TRUE,d,m1,
	    y1c,y2c,rcirc,crit_points);
	  /* at rmax, tangent vector is (-1.0,0.0) */
	  add_crit(rcmin,rmax,ty1c,ty2c,-1.0,0.0,TRUE,d,m1,
		y1c,y2c,rcirc,crit_points);
	  rmin = r_d1;
	  rmax = r_pm;
	  ra = get_rootdef(offdiag,rmin,rmax,FALSE,d,m1,-1,RTPREC);
	  get_rootdef_lu(cspf,rmin,ra,FALSE,d,m1,1,
	    &rcusp_l,&rcusp_u,CUSP_PREC);
	  cval = get_cos(rcusp_l,FALSE,d,m1);
	  x1x2(rcusp_l,cval,d,m1,&x1,&x2);
	  tangent_y(x1,x2,d,m1,&ty1,&ty2);
	  cval = get_cos(rmin,FALSE,d,m1);
	  x1x2(rmin,cval,d,m1,&x1,&x2);
	  tangent_y(x1,x2,d,m1,&ty1min,&ty2min);
	  curv_min = curvature_der(rmin,FALSE,d,m1);
	  if (curv_min > 0.0) {
	    rcmin = get_rootdef(curvature_der,rmin,rcusp_l,FALSE,d,m1,-1,RTPREC);
	    if (perp_caust(rcmin,FALSE,d,m1,y1c,y2c) == 0.0) 
	      add_to_list(crit_points,rcmin,FALSE,d,m1,FALSE);
	    if (perp_caust(rcmin,FALSE,d,m1,y1c,-y2c) == 0.0) 
	      add_to_list(crit_points,rcmin,FALSE,d,m1,TRUE);
	    cval = get_cos(rcmin,FALSE,d,m1);
	    x1x2(rcmin,cval,d,m1,&x1,&x2);
	    tangent_y(x1,x2,d,m1,&ty1c,&ty2c);
	    add_crit(rmin,rcmin,ty1min,ty2min,ty1c,ty2c,FALSE,d,m1,
		y1c,y2c,rcirc,crit_points);
	    add_crit(rcmin,rcusp_l,ty1c,ty2c,ty1,ty2,FALSE,d,m1,
		y1c,y2c,rcirc,crit_points);
	  } else {
	    add_crit(rmin,rcusp_l,ty1min,ty2min,ty1,ty2,FALSE,d,m1,
		y1c,y2c,rcirc,crit_points);
	  }
	  cval = get_cos(rcusp_u,FALSE,d,m1);
	  x1x2(rcusp_u,cval,d,m1,&x1,&x2);
	  tangent_y(x1,x2,d,m1,&ty1,&ty2);
	  /* there must be a absolute curvature minimum */
	  rcmin = get_rootdef(curvature_der,rcusp_u,rmax,FALSE,d,m1,-1,RTPREC);
	  if (perp_caust(rcmin,FALSE,d,m1,y1c,y2c) == 0.0) 
	    add_to_list(crit_points,rcmin,FALSE,d,m1,FALSE);
	  if (perp_caust(rcmin,FALSE,d,m1,y1c,-y2c) == 0.0) 
	    add_to_list(crit_points,rcmin,FALSE,d,m1,TRUE);
	  cval = get_cos(rcmin,FALSE,d,m1);
	  x1x2(rcmin,cval,d,m1,&x1,&x2);
	  tangent_y(x1,x2,d,m1,&ty1c,&ty2c);
	  cval = get_cos(rcusp_u,FALSE,d,m1);
	  x1x2(rcusp_u,cval,d,m1,&x1,&x2);
	  tangent_y(x1,x2,d,m1,&ty1,&ty2);
	  add_crit(rcusp_u,rcmin,ty1,ty2,ty1c,ty2c,FALSE,d,m1,
		y1c,y2c,rcirc,crit_points);
	  /* at rmax, tangent vector is (-1.0,0.0) */
	  add_crit(rcmin,rmax,ty1c,ty2c,-1.0,0.0,FALSE,d,m1,
		y1c,y2c,rcirc,crit_points);
 	}	
}

static FLOAT offdiag(FLOAT r, boolean pmflag, FLOAT d, FLOAT m1)
{
	FLOAT cval;
	FLOAT x1,x2,y11,y12,y22;

	cval = get_cos(r,pmflag,d,m1);
	x1x2(r,cval,d,m1,&x1,&x2);
	binlensd1(x1,x2,d,m1,&y11,&y12,&y22);
	return(y12);
}

static FLOAT cspf(FLOAT r, boolean pmflag, FLOAT d, FLOAT m1)
{
	FLOAT cval,x1,x2;
	FLOAT vz1, vz2;
	FLOAT y11,y12,y22,y111,y112,y122,y222;
	FLOAT gd1, gd2;

	cval = get_cos(r,pmflag,d,m1);
	x1x2(r,cval,d,m1,&x1,&x2);
	binlensd1(x1,x2,d,m1,&y11,&y12,&y22);
	binlensd2(x1,x2,d,m1,&y111,&y112,&y122,&y222);
	gd1 = y111*y22+y122*y11-2.0*y112*y12;
	gd2 = y112*y22+y222*y11-2.0*y122*y12;
	vz1 = y22;
	vz2 = -y12;
	return(vz1*gd1+vz2*gd2);
}



static void normal (FLOAT x1, FLOAT x2, FLOAT d, FLOAT m1,
	FLOAT *n1, FLOAT *n2)
{
	FLOAT y11,y12,y22,y111,y112,y122,y222;
	FLOAT gd1, gd2;
	FLOAT norm;

	binlensd1(x1,x2,d,m1,&y11,&y12,&y22);
	binlensd2(x1,x2,d,m1,&y111,&y112,&y122,&y222);
	gd1 = y111*y22+y122*y11-2.0*y112*y12;
	gd2 = y112*y22+y222*y11-2.0*y122*y12;
	norm = sqrt(gd1*gd1+gd2*gd2);
	norm = 1.0;
	*n1 = gd1/norm;
	*n2 = gd2/norm;
}

static void tangent_y(FLOAT x1, FLOAT x2, FLOAT d, FLOAT m1, 
	FLOAT *ty1, FLOAT *ty2)
{
	FLOAT n1, n2;	
	FLOAT y11,y12,y22;
// 	FLOAT norm;

        normal(x1,x2,d,m1,&n1,&n2);
	binlensd1(x1,x2,d,m1,&y11,&y12,&y22);
	*ty1 = n1*y12-n2*y11;
	*ty2 = n1*y22-n2*y12;
	/*
	norm = sqrt(SQ(*ty1)+SQ(*ty2));
	if (norm > 0) { 
	  *ty1 /= norm;
	  *ty2 /= norm;
	} 
	*/
}

static void binlens(FLOAT x1, FLOAT x2, FLOAT d, FLOAT m1, 
	FLOAT *y1, FLOAT *y2)
{
	FLOAT m2;
	FLOAT x2sq;
	FLOAT del11,del12,n1,n2;
	
	m2 = 1.0-m1;
	if (m2 < 0.0) {
	  m2 = 0.0;
	  m1 = 1.0;
	}
	del11 = x1+m2*d;
	del12 = x1-m1*d;
	x2sq = x2*x2;
	n1 = del11*del11+x2sq; 
	n2 = del12*del12+x2sq;
	if (n1 == 0.0 || n2 == 0.0) {
	  *y1 = -1.0e30;
	  *y2 = -1.0e30;
	  return;
	}
	*y1 = x1-m1*del11/n1-m2*del12/n2;
	*y2 = x2-m1*x2/n1-m2*x2/n2;
}

static void binlensd1(FLOAT x1, FLOAT x2, FLOAT d, FLOAT m1, 
	FLOAT *y11, FLOAT *y12, FLOAT *y22)
{
	FLOAT m2;
	FLOAT del11,del12,x2sq,n1,n2,n1sq,n2sq;
	FLOAT sum,sum11,sum12,sum22;

	m2 = 1.0-m1;
	if (m2 < 0.0) {
	  m2 = 0.0;
	  m1 = 1.0;
	}
	del11 = x1+m2*d;
	del12 = x1-m1*d;
	x2sq = x2*x2;
	n1 = del11*del11+x2sq; 
	n2 = del12*del12+x2sq;
	n1sq = n1*n1;
	n2sq = n2*n2;
	sum = m1/n1+m2/n2;
	sum11 = 2.0*m1*del11*del11/n1sq+2.0*m2*del12*del12/n2sq;
	sum12 = 2.0*m1*del11*x2/n1sq+2.0*m2*del12*x2/n2sq;
	sum22 = 2.0*m1*x2sq/n1sq+2.0*m2*x2sq/n2sq;
	*y11 = 1.0-sum+sum11;
	*y12 = sum12;
	*y22 = 1.0-sum+sum22;
}

static void binlensd2(FLOAT x1, FLOAT x2, FLOAT d, FLOAT m1, 
	FLOAT *y111, FLOAT *y112, FLOAT *y122, FLOAT *y222)
{
	FLOAT m2;
	FLOAT del11,del12,x2sq,n1,n2,n1sq,n2sq;
	FLOAT sum1,sum2,sum111,sum112,sum122,sum222;
	FLOAT del11sq, del12sq;

	m2 = 1.0-m1;
	if (m2 < 0.0) {
	  m2 = 0.0;
	  m1 = 1.0;
	}
	del11 = x1+m2*d;
	del12 = x1-m1*d;
	del11sq = del11*del11;
	del12sq = del12*del12;
	x2sq = x2*x2;
	n1 = del11*del11+x2sq; 
	n2 = del12*del12+x2sq;
	n1sq = n1*n1;
	n2sq = n2*n2;
	sum1 = m1/n1sq*del11+m2/n2sq*del12;
	sum2 = m1/n1sq*x2+m2/n2sq*x2;
	sum111 = m1*del11sq/n1sq*del11/n1+m2*del12sq/n2sq*del12/n2;
	sum112 = m1*del11sq/n1sq*x2/n1+m2*del12sq/n2sq*x2/n2;
	sum122 = m1*del11/n1*x2sq/n1sq+m2*del12/n2*x2sq/n2sq;
	sum222 = m1*x2/n1*x2sq/n1sq+m2*x2/n2*x2sq/n2sq;
	*y111 = 6.0*sum1-8.0*sum111;
	*y112 = 2.0*sum2-8.0*sum112;
	*y122 = 2.0*sum1-8.0*sum122;
	*y222 = 6.0*sum2-8.0*sum222;
}

static void binlensd3(FLOAT x1, FLOAT x2, FLOAT d, FLOAT m1, 
	FLOAT *y1111, FLOAT *y1112, FLOAT *y1122, FLOAT *y1222, FLOAT *y2222)
{
	FLOAT m2;
	FLOAT del11,del12,x2sq,n1,n2,n1sq,n2sq;
	FLOAT n1cb,n2cb,n1qt,n2qt;
	FLOAT sum,sum11,sum12,sum22,sum1111,sum1112,sum1122,sum1222,sum2222;
	FLOAT del11sq, del12sq;
	FLOAT del11cb,del11qt,del12cb,del12qt;
	FLOAT x2cb,x2qt;

	m2 = 1.0-m1;
	if (m2 < 0.0) {
	  m2 = 0.0;
	  m1 = 1.0;
	}
	del11 = x1+m2*d;
	del12 = x1-m1*d;
	del11sq = del11*del11;
	del12sq = del12*del12;
	del11cb = del11sq*del11;
	del12cb = del12sq*del12;
	del11qt = del11sq*del11sq;
	del12qt = del12sq*del12sq;
	x2sq = x2*x2;
	x2cb = x2sq*x2;
	x2qt = x2sq*x2sq;
	n1 = del11*del11+x2sq; 
	n2 = del12*del12+x2sq;
	n1sq = n1*n1;
	n2sq = n2*n2;
	n1cb = n1sq*n1;
	n2cb = n2sq*n2;
	n1qt = n1sq*n1sq;
	n2qt = n2sq*n2sq;
	sum = m1/n1sq+m2/n2sq;
	sum11 = m1/n1cb*del11sq+m2/n2cb*del12sq;
	sum12 = m1/n1cb*del11*x2+m2/n2cb*del12*x2;
	sum22 = m1/n1cb*x2sq+m2/n2cb*x2sq;
	sum1111 = m1/n1qt*del11qt+m2/n2qt*del12qt;
	sum1112 = m1/n1qt*del11cb*x2+m2/n2qt*del12cb*x2;
	sum1122 = m1/n1qt*del11sq*x2sq+m2/n2qt*del12sq*x2sq;
	sum1222 = m1/n1qt*del11*x2cb+m2/n2qt*del12*x2cb;
	sum2222 = m1/n1qt*x2qt+m2/n2qt*x2qt;
	*y1111 = 6.0*sum-48.0*sum11+48.0*sum1111;
	*y1112 = -24.0*sum12+48.0*sum1112;
	*y1122 = 2.0*sum-8.0*sum11-8.0*sum22+48.0*sum1122;
	*y1222 = -24.0*sum12+48.0*sum1222;
	*y2222 = 6.0*sum-48.0*sum22+48.0*sum2222;
}

static void binlensd4(FLOAT x1, FLOAT x2, FLOAT d, FLOAT m1, 
	FLOAT *y11111, FLOAT *y11112, FLOAT *y11122, 
	FLOAT *y11222, FLOAT *y12222, FLOAT *y22222)
{
	FLOAT m2;
	FLOAT del11,del12,x2sq,n1,n2,n1sq,n2sq;
	FLOAT n1cb,n2cb,n1q4,n2q4,n1q5,n2q5;
	FLOAT sum1,sum2,sum111,sum112,sum122,sum222;
	FLOAT sum11111,sum11112,sum11122,sum11222,sum12222,sum22222;
	FLOAT del11sq, del12sq;
	FLOAT del11cb,del11q4,del11q5,del12cb,del12q4,del12q5;
	FLOAT x2cb,x2q4,x2q5;

	m2 = 1.0-m1;
	if (m2 < 0.0) {
	  m2 = 0.0;
	  m1 = 1.0;
	}
	del11 = x1+m2*d;
	del12 = x1-m1*d;
	del11sq = del11*del11;
	del12sq = del12*del12;
	del11cb = del11sq*del11;
	del12cb = del12sq*del12;
	del11q4 = del11sq*del11sq;
	del12q4 = del12sq*del12sq;
        del11q5 = del11q4*del11;
	del12q5 = del12q4*del12;	
	x2sq = x2*x2;
	x2cb = x2sq*x2;
	x2q4 = x2sq*x2sq;
	x2q5 = x2q4*x2;
	n1 = del11*del11+x2sq; 
	n2 = del12*del12+x2sq;
	n1sq = n1*n1;
	n2sq = n2*n2;
	n1cb = n1sq*n1;
	n2cb = n2sq*n2;
	n1q4 = n1sq*n1sq;
	n2q4 = n2sq*n2sq;
	n1q5 = n1q4*n1;
	n2q5 = n2q4*n2;
	sum1 = m1/n1cb*del11+m2/n2cb*del12;
	sum2 = m1/n1cb*x2+m2/n2cb*x2;
	sum111 = m1/n1q4*del11cb+m2/n2q4*del12cb;
	sum112 = m1/n1q4*del11sq*x2+m2/n2q4*del12sq*x2;
	sum122 = m1/n1q4*del11*x2sq+m2/n2q4*del12*x2sq;
	sum222 = m1/n1q4*x2cb+m2/n2q4*x2cb;
	sum11111 = m1/n1q5*del11q5+m2/n2q5*del12q5;
	sum11112 = m1/n1q5*del11q4*x2+m2/n2q5*del12q4*x2;
	sum11122 = m1/n1q5*del11cb*x2sq+m2/n2q5*del12cb*x2sq;
	sum11222 = m1/n1q5*del11sq*x2cb+m2/n2q5*del12sq*x2cb;
	sum12222 = m1/n1q5*del11*x2q4+m2/n2q5*del12*x2q4;
	sum22222 = m1/n1q5*x2q5+m2/n2q5*x2q5;
	*y11111 = -120.0*sum1+480.0*sum111-384.0*sum11111;
	*y11112 = -24.0*sum2+288.0*sum112-384.0*sum11112;
	*y11122 = -24.0*sum1+48.0*sum111+144.0*sum122-384.0*sum11122;
	*y11222 = -24.0*sum2+48.0*sum222+144.0*sum112-384.0*sum11222;
	*y12222 = -24.0*sum1+288.0*sum122-384.0*sum12222;
	*y22222 = -120.0*sum2+480.0*sum222-384.0*sum22222;
}

static void binlensd5(FLOAT x1, FLOAT x2, FLOAT d, FLOAT m1, 
	FLOAT *y111111, FLOAT *y111112, FLOAT *y111122, 
	FLOAT *y111222, FLOAT *y112222, FLOAT *y122222, FLOAT *y222222)
{
	FLOAT m2;
	FLOAT del11,del12,x2sq,n1,n2,n1sq,n2sq;
	FLOAT n1cb,n2cb,n1q4,n2q4,n1q5,n2q5,n1sx,n2sx;
	FLOAT sum,sum11,sum12,sum22,sum1111,sum1112,sum1122,sum1222,sum2222;
	FLOAT sum111111,sum111112,sum111122,sum111222,
                sum112222,sum122222,sum222222;
	FLOAT del11sq, del12sq;
	FLOAT del11cb,del11q4,del11q5,del11sx,del12cb,del12q4,del12q5,del12sx;
	FLOAT x2cb,x2q4,x2q5,x2sx;

	m2 = 1.0-m1;
	if (m2 < 0.0) {
	  m2 = 0.0;
	  m1 = 1.0;
	}
	del11 = x1+m2*d;
	del12 = x1-m1*d;
	del11sq = del11*del11;
	del12sq = del12*del12;
	del11cb = del11sq*del11;
	del12cb = del12sq*del12;
	del11q4 = del11sq*del11sq;
	del12q4 = del12sq*del12sq;
        del11q5 = del11q4*del11;
	del12q5 = del12q4*del12;
	del11sx = del11cb*del11cb;
	del12sx = del12cb*del12cb;
	x2sq = x2*x2;
	x2cb = x2sq*x2;
	x2q4 = x2sq*x2sq;
	x2q5 = x2q4*x2;
	x2sx = x2cb*x2cb;
	n1 = del11*del11+x2sq; 
	n2 = del12*del12+x2sq;
	n1sq = n1*n1;
	n2sq = n2*n2;
	n1cb = n1sq*n1;
	n2cb = n2sq*n2;
	n1q4 = n1sq*n1sq;
	n2q4 = n2sq*n2sq;
	n1q5 = n1q4*n1;
	n2q5 = n2q4*n2;
	n1sx = n1cb*n1cb;
	n2sx = n2cb*n2cb;
	sum = m1/n1cb+m2/n2cb;
	sum11 = m1/n1q4*del11sq+m2/n2q4*del12sq;
	sum12 = m1/n1q4*del11*x2+m2/n2q4*del12*x2;
	sum22 = m1/n1q4*x2sq+m2/n2q4*x2sq;
	sum1111 = m1/n1q5*del11q4+m2/n2q5*del12q4;
	sum1112 = m1/n1q5*del11cb*x2+m2/n2q5*del12cb*x2;
	sum1122 = m1/n1q5*del11sq*x2sq+m2/n2q5*del12sq*x2sq;
	sum1222 = m1/n1q5*del11*x2cb+m2/n2q5*del12*x2cb;
	sum2222 = m1/n1q5*x2q4+m2/n2q5*x2q4;
	sum111111 = m1/n1sx*del11sx+m2/n2sx*del12sx;
	sum111112 = m1/n1sx*del11q5*x2+m2/n2sx*del12q5*x2;
	sum111122 = m1/n1sx*del11q4*x2sq+m2/n2sx*del12q4*x2sq;
	sum111222 = m1/n1sx*del11cb*x2cb+m2/n2sx*del12cb*x2cb;
	sum112222 = m1/n1sx*del11sq*x2q4+m2/n2sx*del12sq*x2q4;
	sum122222 = m1/n1sx*del11*x2q5+m2/n2sx*del12*x2q5;
	sum222222 = m1/n1sx*x2sx+m2/n2sx*x2sx;
	*y111111 = -120.0*sum+2160.0*sum11-5760.0*sum1111+3840.0*sum111111;
	*y111112 = 720.0*sum12-3840.0*sum1112+3840.0*sum111112;
	*y111122 = -24.0*sum+288.0*sum11+144.0*sum22
	  -384.0*sum1111-2304.0*sum1122+3840.0*sum111122;
	*y111222 = 432.0*sum12-1152.0*sum1112-1152.0*sum1222+3840.0*sum111222;
	*y112222 = -24.0*sum+288.0*sum22+144.0*sum11
	  -384.0*sum2222-2304.0*sum1122+3840.0*sum112222;
	*y122222 = 720.0*sum12-3840.0*sum1222+3840.0*sum122222;
	*y222222 = -120.0*sum+2160.0*sum22-5760.0*sum2222+3840.0*sum222222;
}


void test_der(void)
{
	FLOAT y1111,y1112,y1122,y1222,y2222;
	FLOAT z1111,z1112,z1122,z1222,z2222;
	FLOAT y11111,y11112,y11122,y11222,y12222,y22222;
	FLOAT z11111,z11112,z11122,z11222,z12222,z22222;
	FLOAT y111111,y111112,y111122,y111222,y112222,y122222,y222222;
	FLOAT x1,x2,d,m1;
	FLOAT delta;

	x1 = 0.7;
	x2 = 0.4;
	d = 0.3;
	m1 = 0.45;
	delta = 1e-10;
	binlensd3(x1,x2,d,m1,&y1111,&y1112,&y1122,&y1222,&y2222);
	binlensd3(x1+delta,x2,d,m1,
           &z1111,&z1112,&z1122,&z1222,&z2222);
	binlensd4(x1,x2,d,m1,
	   &y11111,&y11112,&y11122,&y11222,&y12222,&y22222);
	printf("%-8.6le %-8.6le\n",(z1111-y1111)/delta,y11111);
	printf("%-8.6le %-8.6le\n",(z1112-y1112)/delta,y11112);
	printf("%-8.6le %-8.6le\n",(z1122-y1122)/delta,y11122);
	printf("%-8.6le %-8.6le\n",(z1222-y1222)/delta,y11222);
	printf("%-8.6le %-8.6le\n",(z2222-y2222)/delta,y12222);
	printf("\nx2-der\n");
	binlensd3(x1,x2+delta,d,m1,
           &z1111,&z1112,&z1122,&z1222,&z2222);
	printf("%-8.6le %-8.6le\n",(z1111-y1111)/delta,y11112);
	printf("%-8.6le %-8.6le\n",(z1112-y1112)/delta,y11122);
	printf("%-8.6le %-8.6le\n",(z1122-y1122)/delta,y11222);
	printf("%-8.6le %-8.6le\n",(z1222-y1222)/delta,y12222);
	printf("%-8.6le %-8.6le\n",(z2222-y2222)/delta,y22222);
	printf("\n");
	binlensd4(x1,x2,d,m1,&y11111,&y11112,&y11122,&y11222,&y12222,&y22222);
	binlensd4(x1+delta,x2,d,m1,
           &z11111,&z11112,&z11122,&z11222,&z12222,&z22222);
	binlensd5(x1,x2,d,m1,
	   &y111111,&y111112,&y111122,&y111222,&y112222,&y122222,&y222222);
	printf("%-8.6le %-8.6le\n",(z11111-y11111)/delta,y111111);
	printf("%-8.6le %-8.6le\n",(z11112-y11112)/delta,y111112);
	printf("%-8.6le %-8.6le\n",(z11122-y11122)/delta,y111122);
	printf("%-8.6le %-8.6le\n",(z11222-y11222)/delta,y111222);
	printf("%-8.6le %-8.6le\n",(z12222-y12222)/delta,y112222);
	printf("%-8.6le %-8.6le\n",(z22222-y22222)/delta,y122222);
	printf("\nx2-der\n");
	binlensd4(x1,x2+delta,d,m1,
           &z11111,&z11112,&z11122,&z11222,&z12222,&z22222);
	printf("%-8.6le %-8.6le\n",(z11111-y11111)/delta,y111112);
	printf("%-8.6le %-8.6le\n",(z11112-y11112)/delta,y111122);
	printf("%-8.6le %-8.6le\n",(z11122-y11122)/delta,y111222);
	printf("%-8.6le %-8.6le\n",(z11222-y11222)/delta,y112222);
	printf("%-8.6le %-8.6le\n",(z12222-y12222)/delta,y122222);
	printf("%-8.6le %-8.6le\n",(z22222-y22222)/delta,y222222);
}

