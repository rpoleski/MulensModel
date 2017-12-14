#include <stdio.h>
#include <stdlib.h>

#define FLOAT double

FLOAT *vector(int nl,int nh)
{
	FLOAT *v;

	v=(FLOAT *)calloc((unsigned) (nh-nl+1),sizeof(FLOAT));
	if (!v) fprintf(stderr,"allocation failure in vector()");
	return v-nl;
}

int *ivector(int nl,int nh)
{
	int *v;

	v=(int *)calloc((unsigned) (nh-nl+1),sizeof(int));
	if (!v) fprintf(stderr,"allocation failure in ivector()");
	return v-nl;
}

FLOAT **matrix(int nrl,int nrh,int ncl,int nch)
{
	int i;
	FLOAT **m;

	m=(FLOAT **) calloc((unsigned) (nrh-nrl+1),sizeof(FLOAT*));
	if (!m) fprintf(stderr,"allocation failure 1 in matrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(FLOAT *) calloc((unsigned) (nch-ncl+1),sizeof(FLOAT));
		if (!m[i]) fprintf(stderr,"allocation failure 2 in matrix()");
		m[i] -= ncl;
	}
	return m;
}


void free_vector(FLOAT *v,int nl,int nh)
{
	free((v+nl));
}

void free_ivector(int *v,int nl,int nh)
{
	free((v+nl));
}


void free_matrix(FLOAT **m,int nrl,int nrh,int ncl,int nch)
{
	int i;

	for(i=nrh;i>=nrl;i--) 
			free((m[i]+ncl));
	free((m+nrl));
}

