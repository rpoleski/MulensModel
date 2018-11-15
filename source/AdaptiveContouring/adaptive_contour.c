/* adaptive_contour.c */
/* Implementation of adaptive contouring algorithm */
/* MD 24-2-06 */
/* Revision 3-3-06 */
/* Revision 7-3-06 */
/* Revision 9-3-06 */
/* Revision 14-7-06 */
/* Revision 14-9-06 */
/* Revision 29-9-06 */
/* Revision 3-10-06 */
/* Revision RP 19-12-17 */
/* Revision RP 02-02-18 */
/* Revision RP 24-02-18 */
/* Revision RP 15-11-18 */


#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

#ifndef FLOAT
#define FLOAT double
#endif

#ifndef boolean
typedef enum {FALSE,TRUE} boolean;
#endif

#define SQ(a) ((a)*(a))

#define FMAX(a,b) (fabs(a) > fabs(b) ? fabs(a) : fabs(b))

#define MAXLVL 128 
#define LVLOFF 63 

#define BIGG 1e37

typedef struct region {
	int level;
	FLOAT x1,x2;
	int sign1, sign2;
	unsigned long long factor1, factor2;
	unsigned long long fraction1, fraction2;
	struct region *sub[4];
	struct point *points[4];
	struct region *top;
	boolean do_divide;
} REGION, *REGION_PTR;

typedef struct point {
	FLOAT x1,x2;
	int sign1, sign2;
	unsigned long long factor1, factor2;
	unsigned long long fraction1, fraction2;
	struct region *square[4];
	struct point *left, *right, *top, *bottom;
	boolean is_inside;
} POINT, *POINT_PTR;

typedef struct regstack {
	REGION_PTR sq;
	struct regstack *next;
} REGNODE, *REGNODE_PTR;
	
typedef struct intstack {
	int val;
	struct intstack *next;
} INTNODE, *INTNODE_PTR;

typedef struct ptlist {
	FLOAT x1,x2;
	int min_lvl;
	struct ptlist *next;
} PTLIST, *PTLIST_PTR;
	
static int succ[4] = {1,3,0,2};
static int pred[4] = {2,0,3,1};
static int opp[4] = {3,2,1,0};
static int hmirror[4] = {1,0,3,2};
static int vmirror[4] = {2,3,0,1};

FLOAT adapsimpsonf1(FLOAT (*f)(FLOAT,FLOAT,
          FLOAT (*)(), FLOAT*, void (*)(), FLOAT*,
          FLOAT (*)(int,FLOAT*,FLOAT),int,FLOAT*),
        FLOAT a, FLOAT b, FLOAT p2,
        FLOAT (*rho_func)(), FLOAT *ps,
        void (*lenseq_func)(), FLOAT *pl,
        FLOAT (*ld_func)(int,FLOAT*,FLOAT),
        int n, FLOAT *gam, FLOAT tol);
FLOAT adapsimpsonf2(FLOAT (*f)(FLOAT,FLOAT,
          FLOAT (*)(), FLOAT*, void (*)(), FLOAT*,
          FLOAT (*)(int,FLOAT*,FLOAT),int,FLOAT*),
	FLOAT p1, FLOAT a, FLOAT b, 
        FLOAT (*rho_func)(), FLOAT *ps,
        void (*lenseq_func)(), FLOAT *pl,
	FLOAT (*ld_func)(int,FLOAT*,FLOAT),
	int n, FLOAT *gam, FLOAT tol);
static void push (REGNODE_PTR *regstack, REGION_PTR newreg);
static REGION_PTR pop(REGNODE_PTR *regstack);
static void delete_from_stack(REGION_PTR reg, REGNODE_PTR *stack); 
static void push_num (INTNODE_PTR *numstack, int num);
static int pop_num(INTNODE_PTR *numstack);
static REGION_PTR right_square(REGION_PTR dsquare);
static REGION_PTR left_square(REGION_PTR dsquare);
static REGION_PTR top_square(REGION_PTR dsquare);
static REGION_PTR bottom_square(REGION_PTR dsquare);
static boolean differ_between(int side, REGION_PTR reg);
static void get_cut(REGION_PTR reg, int side, FLOAT *x1, FLOAT *x2); 
static void get_cut2(REGION_PTR reg, int side, FLOAT *x1, FLOAT *x2); 
static void deallocate_struct(REGION_PTR sc_struct);
static void init_struct(REGION_PTR *sc_struct, FLOAT sc_unitsize,
	FLOAT zero_x1, FLOAT zero_x2, 
	REGNODE_PTR *ctstack, REGNODE_PTR *divstack,
	boolean (*inside_func)(), FLOAT *ps,
        void (*lenseq_func)(), FLOAT *pl);
static boolean subdivide_square(REGION_PTR sc_struct, 
  REGION_PTR dsquare, FLOAT sc_unitsize, 
  REGNODE_PTR *ctstack, REGNODE_PTR *divstack,
  boolean (*inside_func)(), FLOAT *ps,
  void (*lenseq_func)(), FLOAT *pl);
static void extend_struct(REGION_PTR sc_struct, FLOAT sc_unitsize,
	FLOAT zero_x1, FLOAT zero_x2,
	REGNODE_PTR *ctstack, REGNODE_PTR *divstack,
	boolean (*inside_func)(), FLOAT *ps,
        void (*lenseq_func)(), FLOAT *pl);
static void insert_point(REGION_PTR sc_struct, FLOAT sc_unitsize,
  FLOAT zero_x1, FLOAT zero_x2, REGION_PTR top, FLOAT x1, FLOAT x2,
  boolean status, int min_lvl,
  REGNODE_PTR *ctstack, REGNODE_PTR *divstack,
  boolean (*inside_func)(), FLOAT *ps,
  void (*lenseq_func)(), FLOAT *pl, boolean *extend);
static void set_minlevels_insertions(FLOAT sc_unitsize,
  PTLIST_PTR *pointlist, PTLIST_PTR *holelist, 
  boolean (*inside_func)(), FLOAT *ps,
  void (*lenseq_func)(), FLOAT *pl);
static void get_region_with_point(REGION_PTR top, FLOAT x1, FLOAT x2,
	REGION_PTR *reg1, REGION_PTR *reg2, POINT_PTR *pt);
static void get_area(FLOAT (*rho_func)(), FLOAT *ps,
	void (*lenseq_func)(), FLOAT *pl,
	FLOAT (*ld_func)(int n, FLOAT gam[], FLOAT rho),
	int n, FLOAT gam[], FLOAT ld_acc,
	FLOAT sc_unitsize, FLOAT *area, FLOAT *error,
	REGNODE_PTR *ctstack);
static boolean inside(FLOAT x1, FLOAT x2, boolean (*inside_func)(),
	FLOAT *ps, void (*lenseq_func)(), FLOAT *pl);
static FLOAT varprofile(FLOAT x1, FLOAT x2, 
	FLOAT (*rho_func)(), FLOAT *ps,
	void (*lenseq_func)(), FLOAT *pl,
	FLOAT (*ld_func)(int,FLOAT*,FLOAT), int n, FLOAT gam[]);


static void push (REGNODE_PTR *regstack, REGION_PTR newreg)
{
	REGNODE_PTR new_;

	new_ = (REGNODE_PTR) malloc(sizeof(REGNODE));
	if (new_ == NULL) {
	  fprintf(stderr,"Sorry, not enough memory available\n");
	  exit(1);
	}
	new_->sq = newreg;
	new_->next = *regstack;
	*regstack = new_;
}

static REGION_PTR pop(REGNODE_PTR *regstack)
{
	REGNODE_PTR top;
	REGION_PTR entry;

	top = *regstack;
	entry = top->sq;
	*regstack = top->next;	
	free(top);
	return(entry);
}

static void delete_from_stack(REGION_PTR reg, REGNODE_PTR *stack) 
{
	REGNODE_PTR curr, prev;

	curr = stack[reg->level+LVLOFF];
	prev = NULL;
	while (curr != NULL) {
	  if (curr->sq == reg) {
	    if (prev == NULL) 
	      stack[reg->level+LVLOFF] = curr->next;
	    else
	      prev->next = curr->next; 
	    free(curr);
	    break;
	  }
	  prev = curr;
	  curr = curr->next;
	}

}


static void push_num (INTNODE_PTR *numstack, int num)
{
	INTNODE_PTR new_;

	new_ = (INTNODE_PTR) malloc(sizeof(INTNODE));
	if (new_ == NULL) {
	  fprintf(stderr,"Sorry, not enough memory available\n");
	  exit(1);
	}
	new_->val = num;
	new_->next = *numstack;
	*numstack = new_;
}

static int pop_num(INTNODE_PTR *numstack)
{
	INTNODE_PTR top;
	int entry;

	top = *numstack;
	entry = top->val;
	*numstack = top->next;	
	free(top);
	return(entry);
}

static REGION_PTR right_square(REGION_PTR dsquare)
{
	REGION_PTR top;
	INTNODE_PTR num_stack;
	int num;

	top = dsquare->top;
	num_stack = NULL;
	while (top != NULL) {
	  num = 0;
	  while (top->sub[num] != dsquare) {
	    num++; 
	    if (num > 3) {
	      fprintf(stderr,"Error in right_square\n");
	      exit(1);
	    }
	  }
	  push_num(&num_stack,num);
	  if (num == 0 || num == 2) break; 
	  dsquare = top;
	  top = dsquare->top;
	}
	if (top == NULL) {
	  while (num_stack != NULL) 
	    pop_num(&num_stack);
	  return(NULL);
	}	  
	while (num_stack != NULL) {
	  dsquare = top->sub[hmirror[pop_num(&num_stack)]];
	  if (dsquare == NULL) break;
	  top = dsquare;
	} 
	while (num_stack != NULL) 
	  pop_num(&num_stack);
	return(top);
}

static REGION_PTR left_square(REGION_PTR dsquare)
{
	REGION_PTR top;
	INTNODE_PTR num_stack;
	int num;

	top = dsquare->top;
	num_stack = NULL;
	while (top != NULL) {
	  num = 0;
	  while (top->sub[num] != dsquare) {
	    num++; 
	    if (num > 3) {
              fprintf(stderr,"Error in left_square\n");
	      exit(1);
            }
	  }
	  push_num(&num_stack,num);
	  if (num == 1 || num == 3) break; 
	  dsquare = top;
	  top = dsquare->top;
	}
	if (top == NULL) {
	  while (num_stack != NULL) 
	    pop_num(&num_stack);
	  return(NULL);
	}	  
	while (num_stack != NULL) {
	  dsquare = top->sub[hmirror[pop_num(&num_stack)]];
	  if (dsquare == NULL) break;
	  top = dsquare;
	} 
	while (num_stack != NULL) 
	  pop_num(&num_stack);
	return(top);
}

static REGION_PTR top_square(REGION_PTR dsquare)
{
	REGION_PTR top;
	INTNODE_PTR num_stack;
	int num;

	top = dsquare->top;
	num_stack = NULL;
	while (top != NULL) {
	  num = 0;
	  while (top->sub[num] != dsquare) {
	    num++; 
	    if (num > 3) {
              fprintf(stderr,"Error in top_square\n");
	      exit(1);
	    }
	  }
	  push_num(&num_stack,num);
	  if (num == 0 || num == 1) break; 
	  dsquare = top;
	  top = dsquare->top;
	}
	if (top == NULL) {
	  while (num_stack != NULL) 
	    pop_num(&num_stack);
	  return(NULL);
	}	  
	while (num_stack != NULL) {
	  dsquare = top->sub[vmirror[pop_num(&num_stack)]];
	  if (dsquare == NULL) break;
	  top = dsquare;
	} 
	while (num_stack != NULL) 
	  pop_num(&num_stack);
	return(top);
}

static REGION_PTR bottom_square(REGION_PTR dsquare)
{
	REGION_PTR top;
	INTNODE_PTR num_stack;
	int num;

	top = dsquare->top;
	num_stack = NULL;
	while (top != NULL) {
	  num = 0;
	  while (top->sub[num] != dsquare) {
	    num++; 
	    if (num > 3) {
              fprintf(stderr,"Error in bottom_square\n");
              exit(1);
            }
	  }
	  push_num(&num_stack,num);
	  if (num == 2 || num == 3) break; 
	  dsquare = top;
	  top = dsquare->top;
	}
	if (top == NULL) {
	  while (num_stack != NULL) 
	    pop_num(&num_stack);
	  return(NULL);
	}	  
	while (num_stack != NULL) {
	  dsquare = top->sub[vmirror[pop_num(&num_stack)]];
	  if (dsquare == NULL) break;
	  top = dsquare;
	} 
	while (num_stack != NULL) 
	  pop_num(&num_stack);
	return(top);
}

static boolean differ_between(int side, REGION_PTR reg)
{
	/* check whether inside/outside differs for point in between
	    on 'side' for region 'reg' */

	if (reg == NULL) return(FALSE);
	if (reg->sub[0] == NULL) return(FALSE);
	if (reg->points[side]->is_inside == reg->points[succ[side]]->is_inside)
	  if (reg->sub[side]->points[succ[side]]->is_inside !=
              reg->points[side]->is_inside)
	      return(TRUE);
	if (differ_between(side,reg->sub[side])) return(TRUE);
	if (differ_between(side,reg->sub[succ[side]])) return(TRUE);
	return(FALSE);
}

static void get_cut(REGION_PTR reg, int side, FLOAT *x1, FLOAT *x2)
{
	/* get point between in/out change */
	
	REGION_PTR sq;

	if (reg->points[side]->square[pred[side]] == NULL) {
	  *x1 = 0.5*(reg->points[side]->x1 + reg->points[succ[side]]->x1);
	  *x2 = 0.5*(reg->points[side]->x2 + reg->points[succ[side]]->x2);
	  return;
	}
	/* otherwise: smallest adjoining square may be larger */
	if (reg->points[side]->square[pred[side]]->level <= 
	  reg->level) { 
	  *x1 = 0.5*(reg->points[side]->x1 + reg->points[succ[side]]->x1);
	  *x2 = 0.5*(reg->points[side]->x2 + reg->points[succ[side]]->x2);
	  return;
	}
	/* otherwise: determine adjoing square of same level */
	sq = reg->points[side]->square[pred[side]];
	while (sq->level > reg->level) {
	  sq = sq->top;
	  if (sq == NULL) {
	    fprintf(stderr,"Error: searching for top square unsuccessful\n");
	    exit(1);
	  }
	}
	get_cut2(sq,side,x1,x2);
}

static void get_cut2(REGION_PTR reg, int side, FLOAT *x1, FLOAT *x2)
{
	if (reg->sub[0] == NULL) {
	  *x1 = 0.5*(reg->points[pred[side]]->x1 + reg->points[opp[side]]->x1);
	  *x2 = 0.5*(reg->points[pred[side]]->x2 + reg->points[opp[side]]->x2);
	  return;
	}
	if (reg->sub[pred[side]]->points[opp[side]]->is_inside
             != reg->points[pred[side]]->is_inside) 
	  get_cut2(reg->sub[pred[side]],side,x1,x2);
        else
	  get_cut2(reg->sub[opp[side]],side,x1,x2);
}

static void deallocate_struct(REGION_PTR sc_struct)
{
	int i;
	boolean sub_exists;

	if (sc_struct == NULL) return;
	if (sc_struct->top == NULL) 
	  for (i=0; i<=3; i++)
	    free(sc_struct->points[i]);
	sub_exists = FALSE;
	for (i=0; i<=3; i++) 
	  if (sc_struct->sub[i] != NULL)
	    sub_exists = TRUE; 
	if (sub_exists) {
	  for (i=0; i<=3; i++) {
	    if (i == 3) {
	      if (sc_struct->sub[i]->points[succ[i]]->top == NULL)
	        free(sc_struct->sub[i]->points[succ[i]]);
	    } else if (i == 1) {
	      if (sc_struct->sub[i]->points[succ[i]]->right == NULL)
	        free(sc_struct->sub[i]->points[succ[i]]);
	    }
	    else
	      free(sc_struct->sub[i]->points[succ[i]]);
	  }
	  free(sc_struct->sub[0]->points[3]); 
	  for (i=0; i<=3; i++) 
	    deallocate_struct(sc_struct->sub[i]);
	}
	free(sc_struct);
}

static void init_struct(REGION_PTR *sc_struct, FLOAT sc_unitsize,
	FLOAT zero_x1, FLOAT zero_x2, 
	REGNODE_PTR *ctstack, REGNODE_PTR *divstack,
        boolean (*inside_func)(), FLOAT *ps,
        void (*lenseq_func)(), FLOAT *pl)
{
	int i,j;
	POINT_PTR pt;
	boolean extend;

	*sc_struct = (REGION_PTR) malloc(sizeof(REGION));	
	if (*sc_struct == NULL) {
	  fprintf(stderr,"Sorry, not enough memory available\n");
	  exit(1);
	}
	(*sc_struct)->level = -1;
	(*sc_struct)->x1 = zero_x1;
	(*sc_struct)->x2 = zero_x2;
	(*sc_struct)->sign1 = 1;
	(*sc_struct)->sign2 = 1;
	(*sc_struct)->factor1 = 0;
	(*sc_struct)->fraction1 = 0;
	(*sc_struct)->factor2 = 0;
	(*sc_struct)->fraction2 = 0;
	(*sc_struct)->top = NULL;
	(*sc_struct)->do_divide = FALSE;
	for (i=0; i<=3; i++) {
	  (*sc_struct)->sub[i] = NULL;	
	  (*sc_struct)->points[i] = (POINT_PTR) malloc(sizeof(POINT));	
	  if ((*sc_struct)->points[i] == NULL) {
	    fprintf(stderr,"Sorry, not enough memory available\n");
	    exit(1);
	  }
	  pt = (*sc_struct)->points[i];
	  if (i == 0 || i == 2)
	    pt->sign1 = -1;
	  else
	    pt->sign1 = 1;
	  if (i == 0 || i == 1)
	    pt->sign2 = -1;
	  else
	    pt->sign2 = 1;
	  pt->factor1 = 1;
	  pt->fraction1 = 0;
	  pt->factor2 = 1;
	  pt->fraction2 = 0;
	  pt->x1 = zero_x1+pt->sign1*sc_unitsize;
	  pt->x2 = zero_x2+pt->sign2*sc_unitsize;
	  pt->is_inside = inside(pt->x1,pt->x2,
           inside_func,ps,lenseq_func,pl);
	  for (j=0; j<=3; j++)
            if (i==j)
	      pt->square[j] = *sc_struct; 
	    else
	      pt->square[j] = NULL;
	  pt->top = NULL;
          pt->bottom = NULL;
	  pt->left = NULL;
	  pt->right = NULL;
	} 
	(*sc_struct)->points[0]->right = (*sc_struct)->points[1];
	(*sc_struct)->points[1]->left = (*sc_struct)->points[0];
	(*sc_struct)->points[0]->top = (*sc_struct)->points[2];
	(*sc_struct)->points[2]->bottom = (*sc_struct)->points[0];
	(*sc_struct)->points[1]->top = (*sc_struct)->points[3];
	(*sc_struct)->points[3]->bottom = (*sc_struct)->points[1];
	(*sc_struct)->points[2]->right = (*sc_struct)->points[3];
	(*sc_struct)->points[3]->left = (*sc_struct)->points[2];
	/* subdivide initial square */
	extend = 
          subdivide_square(*sc_struct,*sc_struct,
                    sc_unitsize,ctstack,divstack,
		    inside_func,ps,lenseq_func,pl); 
	for (i=0; i<=3; i++)
	  if ((*sc_struct)->points[0]->is_inside) extend = TRUE;
	if (extend) extend_struct(*sc_struct,sc_unitsize,
          zero_x1,zero_x2,ctstack,divstack,
	  inside_func,ps,lenseq_func,pl);
}


static boolean subdivide_square(REGION_PTR sc_struct,
  REGION_PTR dsquare, FLOAT sc_unitsize, 
  REGNODE_PTR *ctstack, REGNODE_PTR *divstack,
  boolean (*inside_func)(), FLOAT *ps,
  void (*lenseq_func)(), FLOAT *pl)
{
	/* subdivide square pointed to with "dsquare" */
	/* if adjoining squares require subdivision (and are not marked),
	     mark these and put them on the stack */

	REGION_PTR nsq,sq,adjsq,reg; 
	POINT_PTR /*newsq,pt,*/newpt[4],newcenter;
	int i,j,lvl;
	FLOAT delx1, delx2;
	unsigned long long dfac,dfrac;
	int sign1,sign2;
	FLOAT f;
	boolean inserted[4] = {FALSE, FALSE, FALSE, FALSE};
	boolean extend = FALSE;
	boolean ins;
	boolean subdivided;

	dsquare->do_divide = FALSE; /* mark subdivision as done */
	if (dsquare->sub[0] != NULL) return(FALSE);
 		/* nothing to do */ 
	/* Note: multiple requests on divstack
		can happen, no error, no problem */
	lvl = dsquare->level;

	if (lvl >= MAXLVL-1-LVLOFF) {
	  fprintf(stderr,"request to subdivide square of depth %d, which" \
		" would exceed the maximum %d\n",lvl,MAXLVL-1-LVLOFF);
	  exit(1);
	}
	for (i=0; i<=3; i++) {
	  dsquare->sub[i] = (REGION_PTR) malloc(sizeof(REGION));
	  if (dsquare->sub[i] == NULL) {
	    fprintf(stderr,"Sorry, not enough memory available\n");
	    exit(1);
	  }
	  nsq = dsquare->sub[i];
	  nsq->level = lvl + 1;
	  if (dsquare->factor1 == 0 && dsquare->fraction1 == 0
               && dsquare->factor2 == 0 && dsquare->fraction2 == 0) {
	    if (i==0 || i==2) 
              nsq->sign1 = -1;
            else 
              nsq->sign1 = 1;
            if (i==0 || i==1)
              nsq->sign2 = -1;
            else
              nsq->sign2 = 1;
	  } else {
	    nsq->sign1 = dsquare->sign1;
	    nsq->sign2 = dsquare->sign2;
	  }
	  if (i==0 || i==2) { 
	    delx1 = -sc_unitsize;
	    sign1 = -1;
	  }
	  else {
	    delx1 = sc_unitsize;
	    sign1 = 1;
	  }
	  if (i==0 || i==1) { 
	    delx2 = -sc_unitsize;
	    sign2 = -1;
	  }
	  else {
	    delx2 = sc_unitsize;
	    sign2 = 1;
	  }
	  if (nsq->level == -1) {
	    dfac = 1ULL;
	    dfrac = 0;
	  } else if (nsq->level < -1) {
            f = 1ULL<<(-nsq->level-1);
	    delx1 *= f;
	    delx2 *= f; 	
	    dfac = 1ULL<<(-nsq->level-1);
	    dfrac = 0;
	  } else {
	    f = 1ULL<<(nsq->level+1);
	    delx1 /= f;
            delx2 /= f;
	    dfac = 0;
	    dfrac = 1ULL<<(8*sizeof(unsigned long long)-nsq->level-1);
          }
	  nsq->x1 = dsquare->x1 + delx1;
	  nsq->x2 = dsquare->x2 + delx2; 
	  nsq->factor1 = dsquare->factor1 + sign1*nsq->sign1*dfac;
	  nsq->fraction1 = dsquare->fraction1 + sign1*nsq->sign1*dfrac;
	  nsq->factor2 = dsquare->factor2 + sign2*nsq->sign2*dfac;
	  nsq->fraction2 = dsquare->fraction2 + sign2*nsq->sign2*dfrac;
	  nsq->top = dsquare;
	  nsq->do_divide = FALSE;
	  for (j=0; j<=3; j++)
	    nsq->sub[j] = NULL;
	  nsq->points[i] = dsquare->points[i];
	  dsquare->points[i]->square[i] = nsq;
        }
	/* insert newcenter */
	newcenter = (POINT_PTR) malloc(sizeof(POINT));
	if (newcenter == NULL) {
	  fprintf(stderr,"Sorry, not enough memory available\n");
	  exit(1);
	}
        newcenter->x1 = dsquare->x1;
        newcenter->x2 = dsquare->x2;
        newcenter->sign1 = dsquare->sign1;
        newcenter->sign2 = dsquare->sign2;
        newcenter->factor1 = dsquare->factor1;
        newcenter->fraction1 = dsquare->fraction1;
        newcenter->factor2 = dsquare->factor2;
        newcenter->fraction2 = dsquare->fraction2;
        newcenter->is_inside = inside(newcenter->x1,newcenter->x2,
            inside_func,ps,lenseq_func,pl);
	for (i=0; i<=3; i++) { 
	  dsquare->sub[i]->points[opp[i]] = newcenter;
	  newcenter->square[i] = dsquare->sub[opp[i]]; 
	}
	for (i=0; i<=3; i++) {
	  /* insert new[i], or point to existing node */
	  newpt[i] = NULL;
	  if (dsquare->points[i]->square[pred[i]] != NULL) {
	    sq = dsquare->points[i]->square[pred[i]];
	    if (sq->level >= dsquare->level+1) {
	      while (sq->level > dsquare->level+1)
	        sq = sq->top;
	      /* found */
	      newpt[i] = sq->points[succ[succ[i]]];
            }
          }
	  if (newpt[i] == NULL) { /* insert */
	    inserted[i] = TRUE;
	    newpt[i] = (POINT_PTR) malloc(sizeof(POINT));
	    if (newpt[i] == NULL) {
	      fprintf(stderr,"Sorry, not enough memory available\n");
	      exit(1);
	    }
	    newpt[i]->left = NULL;
	    newpt[i]->right = NULL;
	    newpt[i]->top = NULL;
	    newpt[i]->bottom = NULL;
	    if (dsquare->factor1 == 0 && dsquare->fraction1 == 0
                 && dsquare->factor2 == 0 && dsquare->fraction2 == 0) {
	      if (i==2) 
                newpt[i]->sign1 = -1;
              else 
                newpt[i]->sign1 = 1;
              if (i==0)
                newpt[i]->sign2 = -1;
              else
                newpt[i]->sign2 = 1;
	    } else {
	      newpt[i]->sign1 = dsquare->sign1;
	      newpt[i]->sign2 = dsquare->sign2;
	    }
	    sign1 = 0;
	    sign2 = 0;
	    delx1 = 0.0;
	    delx2 = 0.0;
	    if (i==1) {
	      delx1 = sc_unitsize;
	      sign1 = 1;
	    }
	    else if (i==2) { 
	      delx1 = -sc_unitsize;
	      sign1 = -1;
	    }
	    else if (i==0) { 
	      delx2 = -sc_unitsize;
	      sign2 = -1;
	    }
	    else if (i==3) { 
	      delx2 = sc_unitsize;
	      sign2 = 1;
	    }
	    if (dsquare->sub[i]->level == 0) {
	      dfac = 1ULL;
	      dfrac = 0;
	    } else if (dsquare->sub[i]->level < 0) {
              f = 1ULL<<(-dsquare->sub[i]->level);
	      delx1 *= f;
	      delx2 *= f; 	
	      dfac = 1ULL<<(-dsquare->sub[i]->level);
	      dfrac = 0;
	    } else {
	      f = 1ULL<<(dsquare->sub[i]->level);
	      delx1 /= f;
              delx2 /= f;
	      dfac = 0;
	      dfrac = 
		1ULL<<(8*sizeof(unsigned long long)-dsquare->sub[i]->level);
            }
	    newpt[i]->x1 = dsquare->x1 + delx1;
	    newpt[i]->x2 = dsquare->x2 + delx2; 
	    newpt[i]->factor1 = dsquare->factor1 
	      + sign1*newpt[i]->sign1*dfac;
	    newpt[i]->fraction1 = dsquare->fraction1 
	      + sign1*newpt[i]->sign1*dfrac;
	    if (sign1*newpt[i]->sign1 > 0 && newpt[i]->fraction1 == 0ULL) 
	      newpt[i]->factor1 += 1ULL; /* carry bit */ 
	    newpt[i]->factor2 = dsquare->factor2 
	      + sign2*newpt[i]->sign2*dfac;
	    newpt[i]->fraction2 = dsquare->fraction2 
	      + sign2*newpt[i]->sign2*dfrac;
	    if (sign2*newpt[i]->sign2 > 0 && newpt[i]->fraction2 == 0ULL) 
	      newpt[i]->factor2 += 1ULL; /* carry bit */
	    newpt[i]->is_inside = inside(newpt[i]->x1,newpt[i]->x2,
		inside_func,ps,lenseq_func,pl);
	    newpt[i]->square[opp[i]] = NULL;
	    newpt[i]->square[pred[i]] = NULL;
	  }
	}
	for (i=0; i<=3; i++) {
	  newpt[i]->square[i] = dsquare->sub[succ[i]];
	  newpt[i]->square[succ[i]] = dsquare->sub[i]; 
	  dsquare->sub[i]->points[succ[i]] = newpt[i];
	  dsquare->sub[i]->points[pred[i]] = newpt[pred[i]]; 
	}
	dsquare->sub[3]->points[0]->right = 
         dsquare->sub[3]->points[1];
	dsquare->sub[3]->points[1]->left = 
         dsquare->sub[3]->points[0];
 	dsquare->sub[2]->points[3]->bottom =
	 dsquare->sub[2]->points[1];	
 	dsquare->sub[2]->points[1]->top =
	 dsquare->sub[2]->points[3];	
 	dsquare->sub[0]->points[2]->right =
	 dsquare->sub[0]->points[3];	
 	dsquare->sub[0]->points[3]->left =
	 dsquare->sub[0]->points[2];	
 	dsquare->sub[1]->points[2]->bottom =
	 dsquare->sub[1]->points[0];	
 	dsquare->sub[1]->points[0]->top =
	 dsquare->sub[1]->points[2];
	if (inserted[1]) {
	  dsquare->sub[1]->points[3]->top =	
	    dsquare->points[3];	
	  dsquare->points[3]->bottom =
	    dsquare->sub[1]->points[3];
	  dsquare->sub[3]->points[1]->bottom =
	    dsquare->points[1];
	  dsquare->points[1]->top =
 	    dsquare->sub[3]->points[1];
	  /* check for extension */
	  if (dsquare->sub[1]->points[3]->is_inside)
            if ((dsquare->sub[1]->points[3]->sign1 == 1)
                && (sc_struct->points[1]->factor1 == 
                       dsquare->sub[1]->points[3]->factor1)
                && (sc_struct->points[1]->fraction1 == 
                       dsquare->sub[1]->points[3]->fraction1))
	      extend = TRUE;
	  /* check for neighbouring square requiring subdivision
	      (new point has different inside/outside status than surrounding 
		 points)) */  
	  if ((dsquare->points[1]->is_inside == dsquare->points[3]->is_inside)
            && (dsquare->sub[1]->points[3]->is_inside 
             != dsquare->points[1]->is_inside)) {
	    adjsq = right_square(dsquare);
	    if (adjsq != NULL) {
	      if (adjsq->sub[0] != NULL) {
	        fprintf(stderr,"Error: Request for subdividing divided square (right)\n");
	        exit(1);
	      }
	        adjsq->do_divide = TRUE;
	        push(&divstack[adjsq->level+LVLOFF],adjsq);
	    }
	  } 
	}
	if (inserted[3]) {
	  dsquare->sub[3]->points[2]->left =	
	    dsquare->points[2];	
	  dsquare->points[2]->bottom =
	    dsquare->sub[3]->points[2];
	  dsquare->sub[2]->points[3]->right =
	    dsquare->points[3];
	  dsquare->points[3]->left =
 	    dsquare->sub[2]->points[3];
	  /* check for extension */ /* rule 2 */
	  if (dsquare->sub[3]->points[2]->is_inside)
            if ((dsquare->sub[3]->points[2]->sign2 == 1)
                && (sc_struct->points[3]->factor2 == 
                       dsquare->sub[3]->points[2]->factor2)
                && (sc_struct->points[3]->fraction2 == 
                       dsquare->sub[3]->points[2]->fraction2))
	      extend = TRUE;
	  /* check for neighbouring square requiring subdivision
	      (new point has different inside/outside status than surrounding 
		 points)) */ /* rule 3a */  
	  if ((dsquare->points[3]->is_inside == dsquare->points[2]->is_inside)
            && (dsquare->sub[3]->points[2]->is_inside 
             != dsquare->points[3]->is_inside)) {
	    adjsq = top_square(dsquare);
	    if (adjsq != NULL) {
	      if (adjsq->sub[0] != NULL) {
	        fprintf(stderr,"Error: Request for subdividing divided square (top)\n");
	        exit(1);
	      }
	        adjsq->do_divide = TRUE;
	        push(&divstack[adjsq->level+LVLOFF],adjsq);
	    }
	  }
	}
	if (inserted[2]) {
	  dsquare->sub[2]->points[0]->bottom =	
	    dsquare->points[0];	
	  dsquare->points[0]->top =
	    dsquare->sub[2]->points[0];
	  dsquare->sub[0]->points[2]->top =
	    dsquare->points[2];
	  dsquare->points[2]->bottom =
 	    dsquare->sub[0]->points[2];
	  /* check for extension */ /* rule 2 */
	  if (dsquare->sub[2]->points[0]->is_inside)
            if ((dsquare->sub[2]->points[0]->sign1 == -1)
                && (sc_struct->points[2]->factor1 == 
                       dsquare->sub[2]->points[0]->factor1)
                && (sc_struct->points[2]->fraction1 == 
                       dsquare->sub[2]->points[0]->fraction1))
	      extend = TRUE;
	  /* check for neighbouring square requiring subdivision
	      (new point has different inside/outside status than surrounding 
		 points)) */ /* rule 3a */  
	  if ((dsquare->points[2]->is_inside == dsquare->points[0]->is_inside)
            && (dsquare->sub[2]->points[0]->is_inside 
             != dsquare->points[2]->is_inside)) {
	    adjsq = left_square(dsquare);
	    if (adjsq != NULL) {
	      if (adjsq->sub[0] != NULL) {
	        fprintf(stderr,"Error: Request for subdividing divided square (left)\n");
	        exit(1);
	      }
	        adjsq->do_divide = TRUE;
	        push(&divstack[adjsq->level+LVLOFF],adjsq);
	    }
	  } 
	}
	if (inserted[0]) {
	  dsquare->sub[0]->points[1]->right =	
	    dsquare->points[1];	
	  dsquare->points[1]->left =
	    dsquare->sub[0]->points[1];
	  dsquare->sub[1]->points[0]->left =
	    dsquare->points[0];
	  dsquare->points[0]->right =
 	    dsquare->sub[1]->points[0];
	  /* check for extension */ /* rule 2 */
	  if (dsquare->sub[0]->points[1]->is_inside)
            if ((dsquare->sub[0]->points[1]->sign2 == -1)
                && (sc_struct->points[0]->factor2 == 
                       dsquare->sub[0]->points[1]->factor2)
                && (sc_struct->points[0]->fraction2 == 
                       dsquare->sub[0]->points[1]->fraction2))
	      extend = TRUE;
	  /* check for neighbouring square requiring subdivision
	      (new point has different inside/outside status than surrounding 
		 points)) */ /* rule 3a */ 
	  if ((dsquare->points[0]->is_inside == dsquare->points[1]->is_inside)
            && (dsquare->sub[0]->points[1]->is_inside 
             != dsquare->points[0]->is_inside)) {
	    adjsq = bottom_square(dsquare);
	    if (adjsq != NULL) {
	      if (adjsq->sub[0] != NULL) {
	        fprintf(stderr,"Error: Request for subdividing divided square (bottom)\n");
	        exit(1);
	      }
	        adjsq->do_divide = TRUE;
	        push(&divstack[adjsq->level+LVLOFF],adjsq);
	    }
	  } 
	}
	for (i=0; i<=3; i++) {
	  subdivided = FALSE;
	  /* if inside/outside/inside/outside sequence, subdivide
		recursively */ /* rule 1 */
	  if ((dsquare->sub[i]->points[i]->is_inside
               == dsquare->sub[i]->points[opp[i]]->is_inside)	  
	   && (dsquare->sub[i]->points[pred[i]]->is_inside
               == dsquare->sub[i]->points[succ[i]]->is_inside)	  
	   && (dsquare->sub[i]->points[i]->is_inside
               != dsquare->sub[i]->points[succ[i]]->is_inside)) {	  
	    extend |= subdivide_square(sc_struct,dsquare->sub[i],
	      sc_unitsize,ctstack,divstack,
              inside_func,ps,lenseq_func,pl);
	    subdivided = TRUE;
	  }
	  /* subdivide of inner point is on different side that surrounding
		points */  /* rule 3b */
          else {
            if (!inserted[i]) { 
	      /* largest test square is on same level as dsquare->sub[i] */
	      adjsq = dsquare->sub[i]->points[i]->square[pred[i]];
	      if (adjsq == NULL) {
	        fprintf(stderr,"Error: Adjoining square must exist, but does not\n");
	        exit(1);
	      }
	      while (adjsq->level > dsquare->sub[i]->level) {
	        adjsq = adjsq->top;
	        if (adjsq == NULL) {
		  fprintf(stderr,"Error: Missing top square\n");
	          exit(1);
                }
	      }
	      if (differ_between(opp[i],adjsq)) {
	        extend |= subdivide_square(sc_struct,dsquare->sub[i],
                  sc_unitsize,ctstack,divstack,
		  inside_func,ps,lenseq_func,pl);
	        subdivided = TRUE;
              }
	    }
            if (!subdivided && !inserted[pred[i]]) { 
	      /* largest test square is on same level as dsquare->sub[i] */
	      adjsq = dsquare->sub[i]->points[i]->square[succ[i]];
	      if (adjsq == NULL) {
	        fprintf(stderr,"Error: Adjoining square must exist, but does not\n");
	        exit(1);
	      }
	      while (adjsq->level > dsquare->sub[i]->level) {
	        adjsq = adjsq->top;
	        if (adjsq == NULL) {
		  fprintf(stderr,"Error: Missing top square\n");
	          exit(1);
                }
	      }
	      if (differ_between(succ[i],adjsq)) {
	        extend |= subdivide_square(sc_struct,dsquare->sub[i],
		  sc_unitsize,ctstack,divstack,
		  inside_func,ps,lenseq_func,pl);
	        subdivided = TRUE;
              }
            }
	  }
	  if (!subdivided) {  /* check whether new region is to be inserted into
                       ctstack */
	    reg = dsquare->sub[i];
	    ins = reg->points[0]->is_inside;
	    j = 1;
	    while (j<=3) {
	      if (reg->points[j]->is_inside != ins) {
		if (!reg->do_divide) {
	          reg->do_divide = TRUE;
	          push(&ctstack[reg->level+LVLOFF],reg);
		}
	        break;
	      }
	      j++;
	    }
	  }
	}
	return(extend);
}


static void extend_struct(REGION_PTR sc_struct, FLOAT sc_unitsize,
	FLOAT zero_x1, FLOAT zero_x2,
	REGNODE_PTR *ctstack, REGNODE_PTR *divstack,
	boolean (*inside_func)(), FLOAT *ps,
        void (*lenseq_func)(), FLOAT *pl)
{
	REGION_PTR new_reg[4],subreg,sq;
	int i,j,k;
	FLOAT del;
	FLOAT f;
	unsigned long long sdfac,sdfrac;
	unsigned long long dfac,dfrac;
	unsigned long long ldfac,ldfrac;
	POINT_PTR newcorner[4];
	POINT_PTR newmid[4];
	POINT_PTR short_succ[4];
	POINT_PTR short_pred[4];
	int lvl;
	boolean repeat;
	boolean subdivided;
	boolean ins;

      do {
        repeat = FALSE;
	lvl = sc_struct->level;
	if (lvl <= -LVLOFF) {
	  fprintf(stderr,"request to extend square of depth %d, which" \
		" would exceed the minimum %d\n",lvl,-LVLOFF);
	  exit(1);
	}
	del = sc_unitsize;
	if (lvl == -1) {
	  dfac = 1ULL;
	  dfrac = 0;
	} else if (lvl < -1) {
          f = 1ULL<<(-lvl-1);
	  del *= f;
	  dfac = 1ULL<<(-lvl-1);
	  dfrac = 0;
	} else {
	  f = 1ULL<<(lvl+1);
          del /= f;
	  dfac = 0;
	  dfrac = 1ULL<<(8*sizeof(unsigned long long)-lvl-1);
        }
	if (lvl == 0) {
	  ldfac = 1ULL;
	  ldfrac = 0;
	} else if (lvl < 0) {
	  ldfac = 1ULL<<(-lvl);
	  ldfrac = 0;
	} else {
	  ldfac = 0;
	  ldfrac = 1ULL<<(8*sizeof(unsigned long long)-lvl);
        }
	if (lvl == -2) {
	  sdfac = 1ULL;
	  sdfrac = 0;
	} else if (lvl < -2) {
	  sdfac = 1ULL<<(-lvl-2);
	  sdfrac = 0;
	} else {
	  sdfac = 0;
	  sdfrac = 1ULL<<(8*sizeof(unsigned long long)-lvl-2);
        }
	for (i=0; i<=3; i++) {
	  new_reg[i] = (REGION_PTR) malloc(sizeof(REGION));
	  if (new_reg[i] == NULL) {
	    fprintf(stderr,"Sorry, not enough memory available\n");
	    exit(1);
	  }
	  new_reg[i]->level = lvl;
	  if (i==0 || i==2) 
            new_reg[i]->sign1 = -1;
	  else
            new_reg[i]->sign1 = 1;
          if (i==0 || i==1) 
            new_reg[i]->sign2 = -1;
	  else
            new_reg[i]->sign2 = 1;
	  new_reg[i]->x1 = zero_x1 + new_reg[i]->sign1*del;
	  new_reg[i]->x2 = zero_x2 + new_reg[i]->sign2*del; 
	  new_reg[i]->factor1 = dfac;
	  new_reg[i]->fraction1 = dfrac;
	  new_reg[i]->factor2 = dfac;
	  new_reg[i]->fraction2 = dfrac;
	  new_reg[i]->do_divide = FALSE;
	  new_reg[i]->top = sc_struct;
	  new_reg[i]->sub[opp[i]] = sc_struct->sub[i];
	  new_reg[i]->sub[opp[i]]->top = new_reg[i];
	  /* create new points */
	  /* corner point */
	  newcorner[i] = (POINT_PTR) malloc(sizeof(POINT));
	  if (newcorner[i] == NULL) {
	    fprintf(stderr,"Sorry, not enough memory available\n");
	    exit(1);
	  }
	  newcorner[i]->sign1 = new_reg[i]->sign1;
	  newcorner[i]->sign2 = new_reg[i]->sign2;
	  newcorner[i]->x1 = zero_x1 + 2.0*new_reg[i]->sign1*del;
	  newcorner[i]->x2 = zero_x2 + 2.0*new_reg[i]->sign2*del; 
	  newcorner[i]->factor1 = ldfac;
	  newcorner[i]->fraction1 = ldfrac;
	  newcorner[i]->factor2 = ldfac;
	  newcorner[i]->fraction2 = ldfrac;
	  newcorner[i]->is_inside = inside(newcorner[i]->x1,newcorner[i]->x2,
		inside_func,ps,lenseq_func,pl);
	  if (newcorner[i]->is_inside)
            repeat = TRUE;
	  if (i==0 || i==1) 
	    newcorner[i]->bottom = NULL;
	  if (i==1 || i==3)
	    newcorner[i]->right = NULL;
	  if (i==3 || i==2)
	    newcorner[i]->top = NULL;
	  if (i==2 || i==0)
	    newcorner[i]->left = NULL; 
	  /* succeeding far point */
	  newmid[i] = (POINT_PTR) malloc(sizeof(POINT));
	  if (newmid[i] == NULL) {
	    fprintf(stderr,"Sorry, not enough memory available\n");
	    exit(1);
	  }
	  if (i == 2)
	    newmid[i]->sign1 = -1;
	  else
	    newmid[i]->sign1 = 1;
	  if (i == 0)
	    newmid[i]->sign2 = -1;
	  else
	    newmid[i]->sign2 = 1;
	  if (i==0 || i==3) {
	    newmid[i]->x1 = zero_x1;
	    newmid[i]->x2 = zero_x2 + 2.0*newmid[i]->sign2*del; 
	    newmid[i]->factor1 = 0;
	    newmid[i]->fraction1 = 0;
	    newmid[i]->factor2 = ldfac;
	    newmid[i]->fraction2 = ldfrac;
	  } else {
	    newmid[i]->x1 = zero_x1 + 2.0*newmid[i]->sign1*del;
	    newmid[i]->x2 = zero_x2; 
	    newmid[i]->factor1 = ldfac;
	    newmid[i]->fraction1 = ldfrac;
	    newmid[i]->factor2 = 0;
	    newmid[i]->fraction2 = 0;
	  }
	  newmid[i]->is_inside = inside(newmid[i]->x1,newmid[i]->x2,
		inside_func,ps,lenseq_func,pl);
	  if (newmid[i]->is_inside)
            repeat = TRUE;
	  if (i==0) 
	    newmid[i]->bottom = NULL;
	  if (i==1)
	    newmid[i]->right = NULL;
	  if (i==3)
	    newmid[i]->top = NULL;
	  if (i==2)
	    newmid[i]->left = NULL; 
	  /* succeeding and preceeding close points */
	  short_succ[i] = (POINT_PTR) malloc(sizeof(POINT));
	  if (short_succ[i] == NULL) {
	    fprintf(stderr,"Sorry, not enough memory available\n");
	    exit(1);
	  }
	  short_pred[i] = (POINT_PTR) malloc(sizeof(POINT));
	  if (short_pred[i] == NULL) {
	    fprintf(stderr,"Sorry, not enough memory available\n");
	    exit(1);
	  }
	  if (i==0 || i==2) {
	    short_succ[i]->sign1 = -1;
	    short_pred[i]->sign1 = -1;
	  } else {
	    short_succ[i]->sign1 = 1;
	    short_pred[i]->sign1 = 1;
	  }
	  if (i==0 || i==1) {
	    short_succ[i]->sign2 = -1;
	    short_pred[i]->sign2 = -1;
	  } else {
	    short_succ[i]->sign2 = 1;
	    short_pred[i]->sign2 = 1;
	  }
	  if (i==0 || i==3) {
	    short_succ[i]->x1 = zero_x1 + short_succ[i]->sign1*del;
	    short_succ[i]->x2 = zero_x2 + short_succ[i]->sign2*2.0*del; 
	    short_succ[i]->factor1 = dfac;
	    short_succ[i]->fraction1 = dfrac;
	    short_succ[i]->factor2 = ldfac;
	    short_succ[i]->fraction2 = ldfrac;
	    short_pred[i]->x1 = zero_x1 + short_pred[i]->sign1*2.0*del;
	    short_pred[i]->x2 = zero_x2 + short_pred[i]->sign2*del; 
	    short_pred[i]->factor1 = ldfac;
	    short_pred[i]->fraction1 = ldfrac;
	    short_pred[i]->factor2 = dfac;
	    short_pred[i]->fraction2 = dfrac;
	  } else {
	    short_succ[i]->x1 = zero_x1 + short_succ[i]->sign1*2.0*del;
	    short_succ[i]->x2 = zero_x2 + short_succ[i]->sign2*del; 
	    short_succ[i]->factor1 = ldfac;
	    short_succ[i]->fraction1 = ldfrac;
	    short_succ[i]->factor2 = dfac;
	    short_succ[i]->fraction2 = dfrac;
	    short_pred[i]->x1 = zero_x1 + short_pred[i]->sign1*del;
	    short_pred[i]->x2 = zero_x2 + short_pred[i]->sign2*2.0*del; 
	    short_pred[i]->factor1 = dfac;
	    short_pred[i]->fraction1 = dfrac;
	    short_pred[i]->factor2 = ldfac;
	    short_pred[i]->fraction2 = ldfrac;
	  }
	  short_succ[i]->is_inside = 
	     inside(short_succ[i]->x1,short_succ[i]->x2,
		inside_func,ps,lenseq_func,pl);
	  short_pred[i]->is_inside = 
	     inside(short_pred[i]->x1,short_pred[i]->x2,
		inside_func,ps,lenseq_func,pl);
	  if (short_succ[i]->is_inside)
	    repeat = TRUE;
	  if (short_pred[i]->is_inside)
	    repeat = TRUE;
	  if (i==0) { 
	    short_succ[i]->bottom = NULL;
	    short_pred[i]->left = NULL;
	  } else if (i==1) {
	    short_succ[i]->right = NULL;
	    short_pred[i]->bottom = NULL;
	  } else if (i==3) {
	    short_succ[i]->top = NULL;
	    short_pred[i]->right = NULL;
	  } else if (i==2) {
	    short_succ[i]->left = NULL;
	    short_pred[i]->top = NULL;
	  }
	  /* make new subregions */
	  for (j=0; j<=3; j++) {
	    if (j != opp[i]) {
	      new_reg[i]->sub[j] = (REGION_PTR) malloc(sizeof(REGION));
	      if (new_reg[i] == NULL) {
	        fprintf(stderr,"Sorry, not enough memory available\n");
	        exit(1);
	      }
	      subreg = new_reg[i]->sub[j];
	      subreg->top = new_reg[i];
	      subreg->level = new_reg[i]->level + 1;
	      subreg->do_divide = FALSE;
	      for (k=0; k<=3; k++)
                subreg->sub[k] = NULL; 
	      if (i==0 || i==2) 
		subreg->sign1 = -1;
	      else
		subreg->sign1 = 1;
	      if (i==0 || i==1)
	        subreg->sign2 = -1;
	      else
		subreg->sign2 = 1;
	      if (i == j) {
		subreg->x1 = zero_x1 + subreg->sign1*1.5*del;
	        subreg->x2 = zero_x2 + subreg->sign2*1.5*del; 
	        subreg->factor1 = dfac+sdfac;
		subreg->fraction1 = dfrac+sdfrac;
		subreg->factor2 = dfac+sdfac;
	        subreg->fraction2 = dfrac+sdfrac;
	      } else if ((i==0 && j==1) || (i==1 && j==0)
                           || (i==2 && j==3) || (i==3 && j==2)) {
		subreg->x1 = zero_x1 + subreg->sign1*0.5*del;
	        subreg->x2 = zero_x2 + subreg->sign2*1.5*del; 
	        subreg->factor1 = sdfac;
		subreg->fraction1 = sdfrac;
		subreg->factor2 = dfac+sdfac;
	        subreg->fraction2 = dfrac+sdfrac;
	      } else {
		subreg->x1 = zero_x1 + subreg->sign1*1.5*del;
	        subreg->x2 = zero_x2 + subreg->sign2*0.5*del; 
	        subreg->factor1 = dfac+sdfac;
		subreg->fraction1 = dfrac+sdfrac;
		subreg->factor2 = sdfac;
	        subreg->fraction2 = sdfrac;
	      }
	    }
	  }
	}
	for (i=0; i<=3; i++) {
	  new_reg[i]->points[i] = newcorner[i];
	  new_reg[i]->points[succ[i]] = newmid[i];
	  new_reg[i]->points[pred[i]] = newmid[pred[i]]; 
	  new_reg[i]->points[opp[i]] = sc_struct->sub[i]->points[opp[i]];
	  new_reg[i]->sub[i]->points[i] = newcorner[i];
	  new_reg[i]->sub[i]->points[succ[i]] = short_succ[i];
	  new_reg[i]->sub[i]->points[pred[i]] = short_pred[i];
	  new_reg[i]->sub[i]->points[opp[i]] = sc_struct->points[i];
	  newcorner[i]->square[i] = new_reg[i]->sub[i];
	  short_succ[i]->square[succ[i]] = new_reg[i]->sub[i];
	  short_pred[i]->square[pred[i]] = new_reg[i]->sub[i];
	  sc_struct->points[i]->square[opp[i]] = new_reg[i]->sub[i];
	  new_reg[i]->sub[succ[i]]->points[i] = short_succ[i];
	  new_reg[i]->sub[succ[i]]->points[opp[i]] = 
	    sc_struct->sub[i]->points[succ[i]];
	  new_reg[i]->sub[succ[i]]->points[succ[i]] = newmid[i];
	  new_reg[i]->sub[succ[i]]->points[pred[i]] = sc_struct->points[i];
	  short_succ[i]->square[i] = new_reg[i]->sub[succ[i]];
	  sc_struct->sub[i]->points[succ[i]]->square[opp[i]]
            = new_reg[i]->sub[succ[i]];
	  newmid[i]->square[succ[i]] = new_reg[i]->sub[succ[i]];    
	  sc_struct->points[i]->square[pred[i]] = new_reg[i]->sub[succ[i]];
	  new_reg[i]->sub[pred[i]]->points[i] = short_pred[i];
	  new_reg[i]->sub[pred[i]]->points[opp[i]] = 
	   sc_struct->sub[i]->points[pred[i]];
	  new_reg[i]->sub[pred[i]]->points[succ[i]] = sc_struct->points[i];
	  new_reg[i]->sub[pred[i]]->points[pred[i]] = newmid[pred[i]];
	  short_pred[i]->square[i] = new_reg[i]->sub[pred[i]];
	  sc_struct->sub[i]->points[pred[i]]->square[opp[i]]
            =  new_reg[i]->sub[pred[i]];
	  sc_struct->points[i]->square[succ[i]] =  new_reg[i]->sub[pred[i]];
	  newmid[pred[i]]->square[pred[i]] =  new_reg[i]->sub[pred[i]];
	  newmid[i]->square[opp[i]] = NULL;
	  newmid[i]->square[pred[i]] = NULL;
	  short_succ[i]->square[opp[i]] = NULL;
	  short_succ[i]->square[pred[i]] = NULL;
	  short_pred[i]->square[opp[i]] = NULL;
	  short_pred[i]->square[succ[i]] = NULL;
	  for (j=0; j<=3; j++) {
            if (i != j) newcorner[i]->square[j] = NULL;
	    /* already done by new_reg[i]->sub[opp[i]] = sc_struct->sub[i] */
	    /*
	    new_reg[i]->sub[opp[i]]->points[j] = sc_struct->sub[i]->points[j];
	    if (sc_struct->sub[i]->points[j]->square[j]
                   == sc_struct->sub[i])
	      sc_struct->sub[i]->points[j]->square[j] = new_reg[i]->sub[opp[i]];
	    */
	  }
	}
	newcorner[0]->right = short_succ[0];
	short_succ[0]->right = newmid[0];
	newmid[0]->right = short_pred[1];
	short_pred[1]->right = newcorner[1];
	newcorner[1]->top = short_succ[1];
	short_succ[1]->top = newmid[1];
	newmid[1]->top = short_pred[3];
	short_pred[3]->top = newcorner[3];
	newcorner[3]->left = short_succ[3];
	short_succ[3]->left = newmid[3];
	newmid[3]->left = short_pred[2];
	short_pred[2]->left = newcorner[2];
	newcorner[2]->bottom = short_succ[2];
	short_succ[2]->bottom = newmid[2];
	newmid[2]->bottom = short_pred[0];
	short_pred[0]->bottom = newcorner[0];
	newcorner[0]->top = short_pred[0];
	short_pred[0]->top = newmid[2];
	newmid[2]->top = short_succ[2];
	short_succ[2]->top = newcorner[2];
	newcorner[2]->right = short_pred[2];
	short_pred[2]->right = newmid[3];
	newmid[3]->right = short_succ[3];
	short_succ[3]->right = newcorner[3];
	newcorner[3]->bottom = short_pred[3];
	short_pred[3]->bottom = newmid[1];
	newmid[1]->bottom = short_succ[1];
	short_succ[1]->bottom = newcorner[1];
	newcorner[1]->left = short_pred[1];
	short_pred[1]->left = newmid[0];
	newmid[0]->left = short_succ[0];
	short_succ[0]->left = newcorner[0];
	sc_struct->sub[0]->points[0]->left = short_pred[0];
	short_pred[0]->right = sc_struct->sub[0]->points[0];
	sc_struct->sub[0]->points[0]->bottom = short_succ[0];	
	short_succ[0]->top = sc_struct->sub[0]->points[0];
	sc_struct->sub[0]->points[1]->bottom = newmid[0];
	newmid[0]->top = sc_struct->sub[0]->points[1];
	sc_struct->sub[1]->points[1]->bottom = short_pred[1];
	short_pred[1]->top = sc_struct->sub[1]->points[1];
	sc_struct->sub[1]->points[1]->right = short_succ[1];	
	short_succ[1]->left = sc_struct->sub[1]->points[1];
	sc_struct->sub[1]->points[3]->right = newmid[1];
	newmid[1]->left = sc_struct->sub[1]->points[3];
	sc_struct->sub[3]->points[3]->right = short_pred[3];
	short_pred[3]->left = sc_struct->sub[3]->points[3];
	sc_struct->sub[3]->points[3]->top = short_succ[3];	
	short_succ[3]->bottom = sc_struct->sub[3]->points[3];
	sc_struct->sub[3]->points[2]->top = newmid[3];
	newmid[3]->bottom = sc_struct->sub[3]->points[2];
	sc_struct->sub[2]->points[2]->top = short_pred[2];
	short_pred[2]->bottom = sc_struct->sub[2]->points[2];
	sc_struct->sub[2]->points[2]->left = short_succ[2];	
	short_succ[2]->right = sc_struct->sub[2]->points[2];
	sc_struct->sub[2]->points[0]->left = newmid[2];
	newmid[2]->right = sc_struct->sub[2]->points[0];
	/* set new subregions */
	for (i=0; i<=3; i++) { 
	  sc_struct->sub[i]->top = new_reg[i];
	  sc_struct->sub[i] = new_reg[i];
	  sc_struct->points[i] = new_reg[i]->points[i];
	}
	sc_struct->level--;
	/* sc_struct needs different points */
	for (i=0; i<=3; i++) {
	  for (j=0; j<=3; j++) {
	    if (j==opp[i]) continue;
	    sq = sc_struct->sub[i]->sub[j]; 
	    subdivided = FALSE;
	    /* if inside/outside/inside/outside sequence, subdivide
		recursively */ /* rule 1 */
	    if ((sq->points[j]->is_inside
               == sq->points[opp[j]]->is_inside)	  
	     && (sq->points[pred[j]]->is_inside
               == sq->points[succ[j]]->is_inside)	  
	     && (sq->points[j]->is_inside
               != sq->points[succ[j]]->is_inside)) {	  
	      repeat |= subdivide_square(sc_struct,sq,
	        sc_unitsize,ctstack,divstack,
		inside_func,ps,lenseq_func,pl);
	      subdivided = TRUE;
	    }
	    /* subdivide of inner point is on different side that surrounding
	 	points */  /* rule 3b */
            else if (j==succ[i]) { 
	      if (differ_between(i,
		   sc_struct->sub[i]->sub[succ[j]])) {
	        repeat |= subdivide_square(sc_struct,sq,
                  sc_unitsize,ctstack,divstack,
		  inside_func,ps,lenseq_func,pl);
	        subdivided = TRUE;
              }
	    }
            else if (j==pred[i]) {
	      if (differ_between(j,
	           sc_struct->sub[i]->sub[pred[j]])) {
	        repeat |= subdivide_square(sc_struct,sq,
                  sc_unitsize,ctstack,divstack,
		  inside_func,ps,lenseq_func,pl);
	        subdivided = TRUE;
              }
            }
	    if (!subdivided) {  
              /* check whether new region is to be inserted into
                         ctstack */
	      sq = sc_struct->sub[i]->sub[j];
	      ins = sq->points[0]->is_inside;
	      k = 1;
	      while (k<=3) {
	        if (sq->points[k]->is_inside != ins) {
	          if (!sq->do_divide) {
	            sq->do_divide = TRUE;
	            push(&ctstack[sq->level+LVLOFF],sq);
	          }
	          break;
	        }
	        k++;
	      }
	    }
	  }
	}
      } while (repeat);
}


static void get_area(FLOAT (*rho_func)(), FLOAT *ps,
	void (*lenseq_func)(), FLOAT *pl,
	FLOAT (*ld_func)(int n, FLOAT gam[], FLOAT rho),
	int n, FLOAT gam[], FLOAT ld_acc,
	FLOAT sc_unitsize, FLOAT *area, FLOAT *error,
	REGNODE_PTR *ctstack)
{
	int i,k;
	FLOAT magarea;
	FLOAT J11,J12,J21,J22;
	FLOAT rterror;
	FLOAT x11,x12,x21,x22;
	REGION_PTR reg;
	REGNODE_PTR regnode, curr, prev;

	/* new version of get_area */
	/* simply get through all regions on stack */
	/* stack should only contain regions with segment of contour */
	/* 3-10-06: omit squares that have been subdivided and delete
		them from stack */
	
	magarea = 0.0;
	*area = 0.0;
	*error = 0.0;
	for (k=0; k<=MAXLVL-1; k++) {
	  regnode = ctstack[k];
	  prev = NULL;
	  while (regnode != NULL) {
	    reg = regnode->sq;
	    if (reg->sub[0] != NULL) { /* 3-10-06 */
	      /* subsquares exist; omit and delete from stack */ 
	      curr = regnode;
	      free(curr);  
	      regnode = regnode->next;
	      if (prev == NULL)
	        ctstack[k] = regnode;
	      else
		prev->next = regnode;
	    }
	    else {
	      for (i=0; i<=3; i++) {
	        if ((reg->points[0]->is_inside != reg->points[1]->is_inside)
	           && (reg->points[1]->is_inside != reg->points[3]->is_inside)
	           && (reg->points[3]->is_inside != reg->points[2]->is_inside)
	           && (reg->points[2]->is_inside != reg->points[0]->is_inside)) {
	          fprintf(stderr,"getarea: in/out/in/out (problem!)\n");
	          exit(1);
	        }
	        if (reg->points[i]->is_inside 
	           && !reg->points[succ[i]]->is_inside)
	          get_cut(reg,i,&x11,&x12);
	        if (!reg->points[i]->is_inside
	           && reg->points[succ[i]]->is_inside)
                  get_cut(reg,i,&x21,&x22);
	      }
	      *area += x11*x22-x12*x21;
	      if (n != -1) {
	        J11 = adapsimpsonf1(varprofile,0,x11,0.5*(x12+x22),
			rho_func,ps,lenseq_func,pl,ld_func,n,gam,ld_acc);
	        J12 = adapsimpsonf2(varprofile,0.5*(x11+x21),0,x12,
		        rho_func,ps,lenseq_func,pl,ld_func,n,gam,ld_acc);
	        J21 = adapsimpsonf1(varprofile,0,x21,0.5*(x12+x22),
			rho_func,ps,lenseq_func,pl,ld_func,n,gam,ld_acc);
	        J22 = adapsimpsonf2(varprofile,0.5*(x11+x21),0,x22,
			rho_func,ps,lenseq_func,pl,ld_func,n,gam,ld_acc);
	        magarea += 0.5*((J11+J21)*(x22-x12)-(J12+J22)*(x21-x11));
	      }	
	      rterror = sc_unitsize;
	      if (reg->level < 0) 
                rterror *= 1ULL<<(-reg->level);
	      else if (reg->level > 0) 
	        rterror /= 1ULL<<(reg->level);
	      *error += SQ(rterror); 
	      prev = regnode;
	      regnode = regnode->next; 
	    }
	  }
	}
	if (n != -1) {
	 *area *= 2.0*(*ld_func)(n,gam,1.0)-(*ld_func)(n,gam,0.0);
	 *area += magarea; 
	}
}


static void get_region_with_point(REGION_PTR top, FLOAT x1, FLOAT x2,
	REGION_PTR *reg1, REGION_PTR *reg2, POINT_PTR *pt)
{
	int i;
	FLOAT x1div, x2div;
	REGION_PTR dum;

	*reg1 = NULL;
	*reg2 = NULL;
	*pt = NULL;
	if (top == NULL) return;
	if (!(top->points[0]->x1 <= x1 && top->points[1]->x1 >= x1
           && top->points[0]->x2 <= x2 && top->points[2]->x2 >= x2))
	  return; /* outside top region */
	for (i=0; i<=3; i++)
          if (top->points[i]->x1 == x1 && top->points[i]->x2 == x2) {
	    *pt = top->points[i];
	    return;
	  }
	if (top->sub[0] == NULL) {
	  *reg1 = top;
	  return;
	}
	x1div = top->sub[0]->points[1]->x1;
	x2div = top->sub[0]->points[3]->x2;
	if (x1 < x1div && x2 < x2div) 
	  get_region_with_point(top->sub[0],x1,x2,reg1,reg2,pt); 
	if (x1 > x1div && x2 < x2div) 
	  get_region_with_point(top->sub[1],x1,x2,reg1,reg2,pt); 
	if (x1 > x1div && x2 > x2div) 
	  get_region_with_point(top->sub[3],x1,x2,reg1,reg2,pt); 
	if (x1 < x1div && x2 > x2div) 
	  get_region_with_point(top->sub[2],x1,x2,reg1,reg2,pt); 
	if (x1 == x1div && x2 == x2div) {
	  *pt = top->sub[0]->points[3];
	  return;
	}
	/* Note: splitting can occur only once */
        if (x1 == x1div) {
	  if (x2 < x2div) {
	    get_region_with_point(top->sub[0],x1,x2,reg1,&dum,pt);
	    if (*pt == NULL)
	      get_region_with_point(top->sub[1],x1,x2,reg2,&dum,pt);	
	  }
	  else {
	    get_region_with_point(top->sub[3],x1,x2,reg1,&dum,pt);
	    if (*pt == NULL)
	      get_region_with_point(top->sub[2],x1,x2,reg2,&dum,pt);	
	  }
        }
	else if (x2 == x2div) {
	  if (x1 < x1div) {
	    get_region_with_point(top->sub[0],x1,x2,reg1,&dum,pt);
	    if (*pt == NULL)
	      get_region_with_point(top->sub[2],x1,x2,reg2,&dum,pt);	
	  }
	  else {
	    get_region_with_point(top->sub[3],x1,x2,reg1,&dum,pt);
	    if (*pt == NULL)
	      get_region_with_point(top->sub[1],x1,x2,reg2,&dum,pt);	
	  }
        }
}


static void insert_point(REGION_PTR sc_struct, FLOAT sc_unitsize,
  FLOAT zero_x1, FLOAT zero_x2,  
  REGION_PTR top, FLOAT x1, FLOAT x2,
  boolean status, int min_lvl,
  REGNODE_PTR *ctstack, REGNODE_PTR *divstack,
  boolean (*inside_func)(), FLOAT *ps,
  void (*lenseq_func)(), FLOAT *pl, boolean *extend)
{
	REGION_PTR reg1,reg2;
	POINT_PTR pt;

	/* require neighbourhood to have "status" with respect to 
	     is_inside */
	/* for ptlist: status = TRUE */
	/* for holelist: status = FALSE */

	if (inside(x1,x2,inside_func,ps,lenseq_func,pl)!=status) return; /* modified 24-05-06 */
	while (!(sc_struct->points[0]->x1 < x1
                 && sc_struct->points[1]->x1 > x1
                 && sc_struct->points[0]->x2 < x2
                 && sc_struct->points[2]->x2 > x2))
	    extend_struct(sc_struct,sc_unitsize,zero_x1,zero_x2,
	      ctstack,divstack,inside_func,ps,lenseq_func,pl);	
	get_region_with_point(top,x1,x2,&reg1,&reg2,&pt);
	if (pt != NULL) return;
	if (!(reg1->points[0]->is_inside == status &&
                reg1->points[1]->is_inside == status &&
                reg1->points[2]->is_inside == status &&
		reg1->points[3]->is_inside == status &&
		reg1->level >= min_lvl)) {
	  /* region may have been already marked for subdivision */
	  delete_from_stack(reg1,ctstack);
	  delete_from_stack(reg1,divstack);
          *extend |= 
	    subdivide_square(sc_struct,reg1,sc_unitsize,ctstack,divstack,
	      inside_func,ps,lenseq_func,pl);
	  insert_point(sc_struct,sc_unitsize,zero_x1,zero_x2,
	    reg1,x1,x2,status,min_lvl,ctstack,divstack,
	    inside_func,ps,lenseq_func,pl,extend);
	}
	if (reg2 != NULL)
	  if (!(reg2->points[0]->is_inside &&
                reg2->points[1]->is_inside &&
                reg2->points[2]->is_inside &&
		reg2->points[3]->is_inside &&
		reg2->level >= min_lvl)) {
	  /* region may have been already marked for subdivision */
	    delete_from_stack(reg2,ctstack);
	    delete_from_stack(reg2,divstack);
            *extend |= 
	      subdivide_square(sc_struct,reg2,sc_unitsize,ctstack,divstack,
		inside_func,ps,lenseq_func,pl);
	    insert_point(sc_struct,sc_unitsize,zero_x1,zero_x2,
	      reg2,x1,x2,status,min_lvl,ctstack,divstack,
	      inside_func,ps,lenseq_func,pl,extend);
	  }
}


static boolean inside(FLOAT x1, FLOAT x2, boolean (*inside_func)(FLOAT, FLOAT, FLOAT*),
	FLOAT *ps, void (*lenseq_func)(FLOAT, FLOAT, FLOAT*, FLOAT*, FLOAT*), FLOAT *pl)
{
	FLOAT y1,y2;
	
	(*lenseq_func)(x1,x2,&y1,&y2,pl);
	return((*inside_func)(y1,y2,ps));
}
	
static FLOAT varprofile(FLOAT x1, FLOAT x2, 
	FLOAT (*rho_func)(FLOAT, FLOAT, FLOAT*), 
        FLOAT *ps, void (*lenseq_func)(FLOAT, FLOAT, FLOAT*, FLOAT*, FLOAT*), FLOAT *pl,
	FLOAT (*ld_func)(int,FLOAT*,FLOAT), int n, FLOAT gam[])
{
	FLOAT y1,y2;
	FLOAT rho;
	FLOAT I0,I1;

	/* modified variable part of the ld-profile, continued for
	    rho > 1, so that it is > 0 and vanishes as rho -> infty */

	(*lenseq_func)(x1,x2,&y1,&y2,pl);
	rho = rho_func(y1,y2,ps);
	I0 = (*ld_func)(n,gam,0.0);
	
	if (rho == 1.0) { 
	  I1 = (*ld_func)(n,gam,1.0);
	  return(I0-I1);
	}
	else if (rho < 1.0) {
	  I1 = (*ld_func)(n,gam,1.0);
	  return((*ld_func)(n,gam,rho)+I0-2.0*I1);
	}
	else  
	  return(I0-(*ld_func)(n,gam,1.0/rho));
// RP: This function can be faster if I0 and I1 are static 
// variables that are changed only when n or gam[] changes. 
}



static void set_minlevels_insertions(FLOAT sc_unitsize,
  PTLIST_PTR *pointlist, PTLIST_PTR *holelist, 
  boolean (*inside_func)(), FLOAT *ps, 
  void (*lenseq_func)(), FLOAT *pl)
{
	PTLIST_PTR curr, prev, next;
	PTLIST_PTR curr2;
	FLOAT dist;
	int min_lvl;

	curr = *pointlist;
	prev = NULL;
	while (curr != NULL) {
	  if (!inside(curr->x1,curr->x2,inside_func,ps,lenseq_func,pl)) {
	    if (prev == NULL) 
	       *pointlist = curr->next;
	    else
	       prev->next = curr->next;
	    next = curr->next;
	    free(curr);
	    curr = next;
	  }
	  else {
	    prev = curr;
	    curr = curr->next;
	  }
	}
	curr = *holelist;
	prev = NULL;
	while (curr != NULL) {
	  if (inside(curr->x1,curr->x2,inside_func,ps,lenseq_func,pl)) {
	    if (prev == NULL) 
	       *holelist = curr->next;
	    else
	       prev->next = curr->next;
	    next = curr->next;
	    free(curr);
	    curr = next;
	  }
	  else {
	    prev = curr;
	    curr = curr->next;
	  }
	}
	
	curr = *pointlist;
	curr->min_lvl = -10000;
	while (curr != NULL) {	
	  curr2 = *pointlist;
	  while (curr2 != NULL) {
	    if (curr2 != curr) {
	      dist = FMAX(curr->x1-curr2->x1,curr->x2-curr2->x2);
	      min_lvl = ceil(-log(dist/(2.0*sc_unitsize))/M_LN2);
	      if (min_lvl > curr->min_lvl)
	        curr->min_lvl = min_lvl;
	    }
	    curr2 = curr2->next;	
	  }
	  curr2 = *holelist;
	  while (curr2 != NULL) {
	    if (curr2 != curr) {
	      dist = FMAX(curr->x1-curr2->x1,curr->x2-curr2->x2);
	      min_lvl = ceil(-log(dist/(2.0*sc_unitsize))/M_LN2);
	      if (min_lvl > curr->min_lvl)
	        curr->min_lvl = min_lvl;
	    }
	    curr2 = curr2->next;	
	  }
	  curr = curr->next;
	}
	curr = *holelist;
	curr->min_lvl = -10000;
	while (curr != NULL) {	
	  curr2 = *pointlist;
	  while (curr2 != NULL) {
	    if (curr2 != curr) {
	      dist = FMAX(curr->x1-curr2->x1,curr->x2-curr2->x2);
	      min_lvl = ceil(-log(dist/(2.0*sc_unitsize))/M_LN2);
	      if (min_lvl > curr->min_lvl)
	        curr->min_lvl = min_lvl;
	    }
	    curr2 = curr2->next;	
	  }
	  curr2 = *holelist;
	  while (curr2 != NULL) {
	    if (curr2 != curr) {
	      dist = FMAX(curr->x1-curr2->x1,curr->x2-curr2->x2);
	      min_lvl = ceil(-log(dist/(2.0*sc_unitsize))/M_LN2);
	      if (min_lvl > curr->min_lvl)
	        curr->min_lvl = min_lvl;
	    }
	    curr2 = curr2->next;	
	  }
	  curr = curr->next;
	}
}



FLOAT adaptive_contour(FLOAT acc, FLOAT ld_acc,
        PTLIST_PTR *pointlist, PTLIST_PTR *holelist,
        FLOAT (*ld_func)(int n, FLOAT gam[], FLOAT rho),
        int n, FLOAT gam[],
	FLOAT (*rho_func)(), boolean (*inside_func)(), void (*lenseq_func)(),
        int ns, int nl, ...)
{
	REGION_PTR scan_struct;
	FLOAT zero_x1, zero_x2;
	FLOAT scan_unitsize;
	va_list ap;
	FLOAT *ps;
	FLOAT *pl;

	boolean extend;
	REGNODE_PTR ctstack[MAXLVL];
	REGNODE_PTR divstack[MAXLVL];
	FLOAT area,error;
	int i,k,l;
	boolean repeat;
	PTLIST_PTR succ;
	FLOAT /*errfrac,*/oldarea;
// 	boolean unchanged;
	
	ps = (FLOAT *) calloc(ns,sizeof(FLOAT));
	pl = (FLOAT *) calloc(nl,sizeof(FLOAT));
	va_start(ap,nl);
	for (i=0; i<ns; i++) 
	  ps[i] = va_arg(ap,FLOAT);
	for (i=0; i<nl; i++)
	  pl[i] = va_arg(ap,FLOAT);
	va_end(ap);

	/* divstack is a collection of stacks for all levels 1..MAXLVL
	    where divstack[i] corresponds to level i+LVLOFF (i = 0 .. MAXLVL-1) */  
	

	scan_struct = NULL;
	zero_x1 = 0.0;
	zero_x2 = 0.0;
	scan_unitsize = 1.0;

	for (l=0; l<=MAXLVL-1; l++) {
	  ctstack[l] = NULL;
	  divstack[l] = NULL;
	}
	init_struct(&scan_struct,scan_unitsize,
	   zero_x1,zero_x2,ctstack,divstack,inside_func,ps,lenseq_func,pl);
	extend = FALSE;
	set_minlevels_insertions(scan_unitsize,pointlist,holelist,
		inside_func,ps,lenseq_func,pl);
	while (*pointlist != NULL) {
  	  insert_point(scan_struct, scan_unitsize, zero_x1, zero_x2,
	    scan_struct,(*pointlist)->x1,(*pointlist)->x2,TRUE,
	    (*pointlist)->min_lvl,
	    ctstack,divstack,
	    inside_func,ps,lenseq_func,pl,&extend);  
	  succ = (*pointlist)->next;
	  free(*pointlist);
	  *pointlist = succ;
	}
	while (*holelist != NULL) {
  	  insert_point(scan_struct, scan_unitsize, zero_x1, zero_x2,
	    scan_struct,(*holelist)->x1,(*holelist)->x2,FALSE,
	    (*holelist)->min_lvl,
	    ctstack,divstack,
	    inside_func,ps,lenseq_func,pl,&extend);  
	  succ = (*holelist)->next;
	  free(*holelist);
	  *holelist = succ;
	}
        if (extend) {
          extend_struct(scan_struct,scan_unitsize,
		zero_x1,zero_x2,ctstack,divstack,
		inside_func,ps,lenseq_func,pl);
          extend = FALSE;
        }
	/* output_struct(scan_struct); */
	k = 0;
	oldarea = 0.0;
	do {
	  extend = FALSE;
	  /* change: find squares to be divided at next level
	      already with subdivide_square */
	  /* work through stacks from highest level 
            -> inserting into next one, until current one reached */
	  do {
	    repeat = FALSE;
	    for (l=0; l<=k+LVLOFF; l++) {
	      while (ctstack[l] != NULL) { 
	        extend |= 
                  subdivide_square(scan_struct,pop(&ctstack[l]),
		    scan_unitsize,ctstack,divstack,
		    inside_func,ps,lenseq_func,pl);
	        if (extend) {
		  l = -1; /* 29-9-06 */
	          extend_struct(scan_struct,scan_unitsize,
			zero_x1,zero_x2,ctstack,divstack,
			inside_func,ps,lenseq_func,pl);
	          extend = FALSE;
	        }
	      }
	    } 
	    for (l=0; l<=MAXLVL-1; l++) {
	      while (divstack[l] != NULL) { 
		repeat = TRUE;
	        extend |= 
                  subdivide_square(scan_struct,pop(&divstack[l]),
		    scan_unitsize,ctstack,divstack,
		    inside_func,ps,lenseq_func,pl);
	        if (extend) {
	          extend_struct(scan_struct,scan_unitsize,
		       zero_x1,zero_x2,ctstack,divstack,
		       inside_func,ps,lenseq_func,pl);
	          extend = FALSE;
	        }
	      }
            }
          } while (repeat);
	  /* for get_area: just go through ALL divstacks (all orders) */
	  /*
	  if (k==4) output_struct(scan_struct); 
	  */
	  get_area(rho_func,ps,lenseq_func,pl,ld_func,n,gam,ld_acc,
		scan_unitsize,&area,&error,ctstack);
	  /*
	  printf("%2d %10.8lf %10.8lf %10.8lf \n",k,area/2.0,error/2.0,
	    fabs(error/area));  
	  */
	  oldarea = area;
	  k++;
	} while (k <= MAXLVL-2-LVLOFF && (area <= 0.0 || fabs(error/area) > acc));
	/* output_struct(scan_struct); */
	for (l=0; l<MAXLVL; l++) {
	  while (ctstack[l] != NULL)
	    pop(&ctstack[l]);
	  while (divstack[l] != NULL)
	    pop(&divstack[l]);
	}
	deallocate_struct(scan_struct);
	scan_struct = NULL;
	free(ps);
	free(pl);
	if (k > MAXLVL-2-LVLOFF) {
	  fprintf(stderr,"Depth exceeded in adaptive_contour\n");
	  exit(1);
	}
	/*
	  printf("%2d %10.8lf %10.8lf %10.8lf \n",k,area/2.0,error/2.0,
	    fabs(error/area));  
	*/
	return(area/2.0);
}

