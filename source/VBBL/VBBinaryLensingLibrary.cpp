// VBBinaryLensing v1.2.1
// This code has been developed by Valerio Bozza, University of Salerno.
// Any use of this code for scientific publications should be acknowledged by a citation to
// V. Bozza, MNRAS 408 (2010) 2188

#ifdef _WIN32
char systemslash='\\';
#else
char systemslash='/';
#endif

#include "VBBinaryLensingLibrary.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef __unmanaged
using namespace VBBinaryLensingLibrary;
#endif

double abs(complex);
complex conj(complex);
complex sqrt(complex);
double real(complex);
double imag(complex);
complex operator+(complex,complex);
complex operator-(complex,complex);
complex operator*(complex,complex);
complex operator/(complex,complex);
complex operator+(complex,double);
complex operator-(complex,double);
complex operator*(complex,double);
complex operator/(complex,double);
complex operator+(double,complex);
complex operator-(double,complex);
complex operator*(double,complex);
complex operator/(double,complex);
complex operator+(int,complex);
complex operator-(int,complex);
complex operator*(int,complex);
complex operator/(int,complex);
complex operator+(complex,int);
complex operator-(complex,int);
complex operator*(complex,int);
complex operator/(complex,int);
complex operator-(complex);


VBBinaryLensing::VBBinaryLensing(){
	Obj[0]=-0.0397317;
	Obj[1]=0.998164;
	Obj[2]=-0.045714;
	// reference is ecliptic with x-axis toward the equinox.
	// axial tilt at J2000 is 23:26:21.406 from JPL fundamental ephemeris
	Eq2000[0]=1;
	Eq2000[1]=Eq2000[2]=Quad2000[0]=North2000[0]=0;
	Quad2000[1]=0.9174820003578725;
	Quad2000[2]=-0.3977772982704228;
	North2000[1]=0.3977772982704228;
	North2000[2]=0.9174820003578725; 
	Tol=1.e-2;
	tsat=0;
	possat=0;
	nsat=0;
	ndatasat=0;
	satellite=0;
	parallaxsystem=0;
	minannuli=1;
}

VBBinaryLensing::~VBBinaryLensing(){
	if(nsat){
		for(int i=0;i<nsat;i++){
			for(int j=0;j<ndatasat[i];j++) free(possat[i][j]);
			free(tsat[i]);
			free(possat[i]);
		}
		free(tsat);
		free(possat);
		free(ndatasat);
	}
}

_sols *VBBinaryLensing::PlotCrit(double a1,double q1){
	complex  a,q,ej,zr[4],x1,x2;
	int NPS=200;
	_sols *CriticalCurves;
	_curve *Prov,*Prov2,*isso;
	double SD,MD,CD,centeroffset;

	a=complex(a1,0.0);
	q=complex(q1,0.0);
	centeroffset=a1/2.0*(1.0-q1)/(1.0+q1);

	CriticalCurves=new _sols;
	for(int i=0;i<4;i++){
		Prov=new _curve;
		CriticalCurves->append(Prov);
	}

	for(int j=0;j<NPS;j++){
		ej=complex (cos(2*j*M_PI/NPS),-sin(2*j*M_PI/NPS));
		complex  coefs[5]={a*a/16.0*(4.0-a*a*ej)*(1.0+q),a*(q-1.0),(q+1.0)*(1.0+a*a*ej/2.0),0.0,-(1.0+q)*ej};
		zroots(coefs,4,zr,1,1.0e-10);
		Prov=CriticalCurves->first;
		for(int i=0;i<4;i++){
			Prov->append(zr[i].re+centeroffset,zr[i].im);
			Prov=Prov->next;
		}
	}

	Prov=CriticalCurves->first;
	while(Prov->next){
		SD=*(Prov->first) - *(Prov->last);
		MD=100.0;
		for(Prov2=Prov->next;Prov2;Prov2=Prov2->next){
			CD=*(Prov2->first) - *(Prov->last);
			if(CD<MD){
				MD=CD;
				isso=Prov2;
			}
		}
		if(MD<SD){
			CriticalCurves->drop(isso);
			Prov->join(isso);
		}else{
			Prov=Prov->next;
		}
	}

// Caustics

	for(Prov=CriticalCurves->last; Prov; Prov=Prov->prev){
		Prov2=new _curve;
		for(_point *scanpoint=Prov->first; scanpoint; scanpoint=scanpoint->next){
			x1=complex (scanpoint->x1-centeroffset,0.0);
			x2=complex (scanpoint->x2,0.0);
			Prov2->append(real(_L1)+centeroffset,real(_L2));
		}
		CriticalCurves->append(Prov2);
	}		
	return CriticalCurves;
}

void VBBinaryLensing::PrintCau(double a, double q){
	_sols *CriticalCurves;
	_curve *scancurve;
	_point *scanpoint;
	FILE *f;
	int ncc;

	CriticalCurves=PlotCrit(a,q);
	f=fopen("outcurves.causticdata","w");
	ncc=CriticalCurves->length/2;
	scancurve=CriticalCurves->first;
	for(int i=0;i<ncc;i++){
		scancurve=scancurve->next;
	}

	for(int i=0;i<ncc;i++){
		fprintf(f,"Curve: %d\n",i+1);
		for(scanpoint=scancurve->first; scanpoint; scanpoint=scanpoint->next){
			fprintf(f,"%lf %lf\n",scanpoint->x1,scanpoint->x2);
		}
		scancurve=scancurve->next;
	}
	fclose(f);
	delete CriticalCurves;
}

void VBBinaryLensing::SetObjectCoordinates(char *modelfile,char *sateltabledir){
	double RA,Dec,dis,hr,mn,sc,phiprec;
	char filename[256];
	int ic;
	FILE *f;

	if(nsat){
		for(int i=0;i<nsat;i++){
			for(int j=0;j<ndatasat[i];j++) free(possat[i][j]);
			free(tsat[i]);
			free(possat[i]);
		}
		free(tsat);
		free(possat);
		free(ndatasat);
	}

	f=fopen(modelfile,"r");
	fscanf(f,"%lf:%lf:%lf",&hr,&mn,&sc);
	RA=(hr+mn/60+sc/3600)*M_PI/12,
	fscanf(f,"%lf:%lf:%lf",&hr,&mn,&sc);
	Dec=(fabs(hr)+mn/60+sc/3600)*M_PI/180;
	if(hr<0) Dec=-Dec;

	for(int i=0;i<3;i++){
		Obj[i]=(cos(RA)*cos(Dec)*Eq2000[i]+sin(RA)*cos(Dec)*Quad2000[i]+sin(Dec)*North2000[i]);
		rad[i]=Eq2000[i];
		tang[i]=North2000[i];
	}
	fclose(f);


// Looking for satellite table files in the specified directory
	sprintf(filename,"%s%csatellite*.txt",sateltabledir,systemslash);
	nsat=0;
	for(unsigned char c=32;c<255;c++){
		filename[strlen(filename)-5]=c;
		f=fopen(filename,"r");
		if(f!=0){
			nsat++;
			fclose(f);
		}
	}


	tsat=(double **)malloc(sizeof(double *)*nsat);
	possat=(double ***)malloc(sizeof(double **)*nsat);
	ndatasat=(int *)malloc(sizeof(int)*nsat);

// Reading satellite table files
	ic=0;
	for(unsigned char c=32;c<255;c++){
		filename[strlen(filename)-5]=c;
		f=fopen(filename,"r");
		if(f!=0){
			int flag2=0;
			long startpos=0;
			char teststring[1000];
			ndatasat[ic]=1;

			// Finding start of data
			while(!feof(f)){
				fscanf(f,"%s",teststring);
				if(!feof(f)){
					fseek(f,1,SEEK_CUR);
					teststring[5]=0;
					if(strcmp(teststring,"$$SOE")==0){
						flag2=1;
						break;
					}
				}
			}
			// Finding end of data
			if(flag2){
				flag2=0;
				startpos=ftell(f);
				while(!feof(f)){
					fscanf(f,"%[^\n]s",teststring);
					if(!feof(f)){
						fseek(f,1,SEEK_CUR);
						teststring[5]=0;
						if(strcmp(teststring,"$$EOE")==0){
							flag2=1;
							break;
						}else{
							ndatasat[ic]++;
						}
					}
				}					
			}
			
			// Allocating memory according to the length of the table
			tsat[ic]=(double *)malloc(sizeof(double)*ndatasat[ic]);
			possat[ic]=(double **)malloc(sizeof(double *)*ndatasat[ic]);
			for(int j=0;j<ndatasat[ic];j++){
				possat[ic][j]=(double *)malloc(sizeof(double)*3);
			}
			ndatasat[ic]--;

			// Reading data
			if(f){
				fseek(f,startpos,SEEK_SET);
				for(int id=0;id<ndatasat[ic];id++){

					if(fscanf(f,"%lf %lf %lf %lf %lf",&(tsat[ic][id]),&RA,&Dec,&dis,&phiprec)==5){
						tsat[ic][id]-=2450000;
						RA*=M_PI/180;
						Dec*=M_PI/180;
						for(int i=0;i<3;i++){
							possat[ic][id][i]=dis*(cos(RA)*cos(Dec)*Eq2000[i]+sin(RA)*cos(Dec)*Quad2000[i]+sin(Dec)*North2000[i]);
						}
					}else{
						ndatasat[ic]=id;
						break;
					}
				}
				fclose(f);
			}

			ic++;
		}
	}

}

void VBBinaryLensing::ComputeParallax(double t,double t0,double *Et){
	static double a0=1.00000261, adot=0.00000562; // Ephemeris from JPL website 
	static double e0=0.01671123, edot=-0.00004392;
	static double inc0=-0.00001531, incdot=-0.01294668;
	static double L0=100.46457166, Ldot=35999.37244981;
	static double om0=102.93768193, omdot=0.32327364;
	static double deg=M_PI/180;
	double a,e,inc,L,om,M,EE,dE,dM;
	double x1,y1,vx,vy,Ear[3],vEar[3];
	double Et0[2],vt0[2],t0old=0.,r,sp,ty,Spit;
	int c=0,ic;
	

	if(t0!=t0old){
		t0old=t0;
		ty=(t0-1545)/36525.0;

		a=a0+adot*ty;
		e=e0+edot*ty;
		inc=(inc0+incdot*ty)*deg;
		L=(L0+Ldot*ty)*deg;
		om=(om0+omdot*ty)*deg;

		M=L-om;
		M-=floor((M+M_PI)/(2*M_PI))*2*M_PI;

		EE=M+e*sin(M);
		dE=1;
		while(fabs(dE)>1.e-8){
			dM=M-(EE-e*sin(EE));
			dE=dM/(1-e*cos(EE));
			EE+=dE;
		}
		x1=a*(cos(EE)-e);
		y1=a*sqrt(1-e*e)*sin(EE);
//		r=a*(1-e*cos(EE));
		vx=-a/(1-e*cos(EE))*sin(EE)*Ldot*deg/36525;
		vy=a/(1-e*cos(EE))*cos(EE)*sqrt(1-e*e)*Ldot*deg/36525;

		Ear[0]=x1*cos(om)-y1*sin(om);
		Ear[1]=x1*sin(om)*cos(inc)+y1*cos(om)*cos(inc);
		Ear[2]=x1*sin(om)*sin(inc)+y1*cos(om)*sin(inc);
		vEar[0]=vx*cos(om)-vy*sin(om);
		vEar[1]=vx*sin(om)*cos(inc)+vy*cos(om)*cos(inc);
		vEar[2]=vx*sin(om)*sin(inc)+vy*cos(om)*sin(inc);

		sp=0;
		switch(parallaxsystem){
			case 1:
				for(int i=0;i<3;i++) sp+=North2000[i]*Obj[i];
				for(int i=0;i<3;i++) rad[i]=-North2000[i]+sp*Obj[i];
				break;
			default:
				for(int i=0;i<3;i++) sp+=Ear[i]*Obj[i];
				for(int i=0;i<3;i++) rad[i]=Ear[i]-sp*Obj[i];
				break;
		}

		r=sqrt(rad[0]*rad[0]+rad[1]*rad[1]+rad[2]*rad[2]);
		rad[0]/=r;
		rad[1]/=r;
		rad[2]/=r;
		tang[0]=rad[1]*Obj[2]-rad[2]*Obj[1];
		tang[1]=rad[2]*Obj[0]-rad[0]*Obj[2];
		tang[2]=rad[0]*Obj[1]-rad[1]*Obj[0];

		Et0[0]=Et0[1]=vt0[0]=vt0[1]=0;
		for(int i=0;i<3;i++){
			Et0[0]+=Ear[i]*rad[i];
			Et0[1]+=Ear[i]*tang[i];
			vt0[0]+=vEar[i]*rad[i];
			vt0[1]+=vEar[i]*tang[i];
		}
	}

	ty=(t-1545)/36525.0;

	a=a0+adot*ty;
	e=e0+edot*ty;
	inc=(inc0+incdot*ty)*deg;
	L=(L0+Ldot*ty)*deg;
	om=(om0+omdot*ty)*deg;

	M=L-om;
	M-=floor((M+M_PI)/(2*M_PI))*2*M_PI;

	EE=M+e*sin(M);
	dE=1;
	while(dE>1.e-8){
		dM=M-(EE-e*sin(EE));
		dE=dM/(1-e*cos(EE));
		EE+=dE;
	}
	x1=a*(cos(EE)-e);
	y1=a*sqrt(1-e*e)*sin(EE);
//	r=a*(1-e*cos(EE));

	Ear[0]=x1*cos(om)-y1*sin(om);
	Ear[1]=x1*sin(om)*cos(inc)+y1*cos(om)*cos(inc);
	Ear[2]=x1*sin(om)*sin(inc)+y1*cos(om)*sin(inc);
	Et[0]=Et[1]=0;
	for(int i=0;i<3;i++){
		Et[0]+=Ear[i]*rad[i];
		Et[1]+=Ear[i]*tang[i];
	}
	Et[0]+=-Et0[0]-vt0[0]*(t-t0);
	Et[1]+=-Et0[1]-vt0[1]*(t-t0);

	if(satellite>0 && satellite <=nsat){
		if(ndatasat[satellite-1]>2){
			int left,right;
			if(t<tsat[satellite-1][0]){
				ic=0;
			}else{
				if(t>tsat[satellite-1][ndatasat[satellite-1]-1]){
					ic=ndatasat[satellite-1]-2;
				}else{
					left=0;
					right=ndatasat[satellite-1]-1;
					while(right-left>1){
						ic=(right+left)/2;
						if(tsat[satellite-1][ic]>t){
							right=ic;
						}else{
							left=ic;
						}
					}
					ic=left;
				}
			}
			ty=t-tsat[satellite-1][ic];
			for(int i=0;i<3;i++){
				Spit=possat[satellite-1][ic][i]*(1-ty)+possat[satellite-1][ic+1][i]*ty;
				Et[0]+=Spit*rad[i];
				Et[1]+=Spit*tang[i];
			}
		}
	}
}

double VBBinaryLensing::PSPLCurve(double *pr,double t){
	double u,u1;
	u1=(t-pr[2])/exp(pr[1]);
	u=exp(pr[0]);
	u=u1*u1+u*u;
	return (u+2)/sqrt(u*(u+4));
}

double VBBinaryLensing::PSPLParallaxCurve(double *pr,double t){
	double Mag=-1.0;
	double Et[2];
	double u0,tE,t0,pai1,pai2,th;
	
	u0=pr[0];
	tE=exp(pr[1]);
	t0=pr[2];
	pai1=pr[3];
	pai2=pr[4];
	th=0;

	ComputeParallax(t,t0,Et);

	// Et is the displacement of the Sun from the Earth relative to the position 
	//it would have with a uniform rectilinear motion. 
	// Its components are Et[0] in the Earth->Sun direction projected on the sky plane
	// Et[1] in the perpendicular direction on the sky plane (Earth->Sun \vec Object)

	// We have (D t, D u) = (pai . Et, pai \vec Et)
	// Then pai1 is the component of the parallax vector parallel to the Earth->Sun direction
	// pai2 is the component orthogonal to this and parallel to Earth->Sun \vec Obj.

	y_1=(u0+pai1*Et[1]-pai2*Et[0])*sin(th)-((t-t0)/tE+pai1*Et[0]+pai2*Et[1])*cos(th);
	y_2=(-u0-pai1*Et[1]+pai2*Et[0])*cos(th)-((t-t0)/tE+pai1*Et[0]+pai2*Et[1])*sin(th);

	//y_1=(u0+pai1*Et[0]+pai2*Et[1])*sin(th)-((t-t0)/tE+pai2*Et[0]-pai1*Et[1])*cos(th);
	//y_2=(-u0-pai1*Et[0]-pai2*Et[1])*cos(th)-((t-t0)/tE+pai2*Et[0]-pai1*Et[1])*sin(th);

    u0=y_1*y_1+y_2*y_2;
	return (u0+2)/sqrt(u0*(u0+4));
}


double VBBinaryLensing::ESPLCurve(double *pr,double t){
	double u,r;
	double n,k,b1,b2,b3,iF,iE,iP;
	r=(t-pr[2])/exp(pr[1]);
	u=exp(pr[0]);
	u=sqrt(r*r+u*u);
	r=exp(pr[3]);
	n=4*u*r/((r+u)*(r+u));
	k=sqrt(4*n/(4+(u-r)*(u-r)));
	b1=(r-u)*(8-r*r+u*u);
	b2=(4+(r-u)*(r-u))*(r+u);
	b3=4*(1+r*r)*(r-u)*(r-u)/(r+u);
	iF=ellf(M_PI_2,k);
	iE=elle(M_PI_2,k);
	iP=ellpi(M_PI_2,n,k);
	return (b1*iF+b2*iE+b3*iP)/M_PI/2/(r*r)/sqrt(4+(r-u)*(r-u));
}


double VBBinaryLensing::ESPLParallaxCurve(double *pr,double t){
	double Mag=-1.0;
	double Et[2];
	double u0,tE,t0,pai1,pai2,th;
	double u,rr;
	double n,k,b1,b2,b3,iF,iE,iP;
	
	u0=pr[0];
	tE=exp(pr[1]);
	t0=pr[2];
	pai1=pr[4];
	pai2=pr[5];
	th=0;

	ComputeParallax(t,t0,Et);

	y_1=(u0+pai1*Et[1]-pai2*Et[0])*sin(th)-((t-t0)/tE+pai1*Et[0]+pai2*Et[1])*cos(th);
	y_2=(-u0-pai1*Et[1]+pai2*Et[0])*cos(th)-((t-t0)/tE+pai1*Et[0]+pai2*Et[1])*sin(th);

	//y_1=(u0+pai1*Et[0]+pai2*Et[1])*sin(th)-((t-t0)/tE+pai2*Et[0]-pai1*Et[1])*cos(th);
	//y_2=(-u0-pai1*Et[0]-pai2*Et[1])*cos(th)-((t-t0)/tE+pai2*Et[0]-pai1*Et[1])*sin(th);

    u=sqrt(y_1*y_1+y_2*y_2);
	rr=exp(pr[3]);
	n=4*u*rr/((rr+u)*(rr+u));
	k=sqrt(4*n/(4+(u-rr)*(u-rr)));
	b1=(rr-u)*(8-rr*rr+u*u);
	b2=(4+(rr-u)*(rr-u))*(rr+u);
	b3=4*(1+rr*rr)*(rr-u)*(rr-u)/(rr+u);
	iF=ellf(M_PI_2,k);
	iE=elle(M_PI_2,k);
	iP=ellpi(M_PI_2,n,k);
	return (b1*iF+b2*iE+b3*iP)/(2*M_PI*(rr*rr)*sqrt(4+(rr-u)*(rr-u)));
}


double VBBinaryLensing::BinaryLightCurve(double *pr,double t){
	double Mag=-1.0,Tolv=Tol;
	int c=0;
	_sols *Images;
	double tn=(t-pr[6])/exp(pr[5]);
	//if((b>-.13)&&(t>0.179)){
	//	c=c;
	//}

	y_1=pr[2]*sin(pr[3])-tn*cos(pr[3]);
	y_2=-pr[2]*cos(pr[3])-tn*sin(pr[3]);
	if(fabs(tn)<10.){
		while((Mag<0.9)&&(c<3)){
			Mag= BinaryMag(exp(pr[0]),exp(pr[1]),y_1,y_2,exp(pr[4]),Tolv,&Images);
			Tolv/=10;
			c++;
			delete Images;
		}
	}else{
		Mag=-Mag;
	}

	return Mag;
}

double VBBinaryLensing::BinaryLightCurveW(double *pr,double t){
	double Mag=-1.0,Tolv=Tol,a,q,xc,t0,u0,tn,tE;
	int c=0;
	_sols *Images;
	a=exp(pr[0]);
	q=exp(pr[1]);
	tE=exp(pr[5]);
	xc=(a-1/a)/(1+exp(pr[1]));
	if(xc<0) xc=0.;
	t0=pr[6]+xc*cos(pr[3])*tE;
	u0=pr[2]+xc*sin(pr[3]);
	tn=(t-t0)/tE;

	y_1=u0*sin(pr[3])-tn*cos(pr[3]);
	y_2=-u0*cos(pr[3])-tn*sin(pr[3]);
	if(fabs(tn)<10.){
		while((Mag<0.9)&&(c<3)){
			Mag= BinaryMag(a,q,y_1,y_2,exp(pr[4]),Tolv,&Images);
			Tolv/=10;
			c++;
			delete Images;
		}
	}else{
		Mag=-Mag;
	}

	return Mag;
}

// ds/dt = a w1. dth/dt = w2



double VBBinaryLensing::BinaryLightCurveParallax(double *pr,double t){
//(double a,double q,double u0,double th,double RSv,double tE,double t0,double t,double pai1,double pai2,double Tol){
	double Mag=-1.0,Tolv=Tol;
	double a=exp(pr[0]),q=exp(pr[1]),u0=pr[2],th=pr[3],RSv=exp(pr[4]),tE=exp(pr[5]),t0=pr[6],pai1=pr[7],pai2=pr[8];
	double Et[2];
	int c=0;
	_sols *Images;
	
	ComputeParallax(t,t0,Et);

	y_1=(u0+pai1*Et[1]-pai2*Et[0])*sin(th)-((t-t0)/tE+pai1*Et[0]+pai2*Et[1])*cos(th);
	y_2=(-u0-pai1*Et[1]+pai2*Et[0])*cos(th)-((t-t0)/tE+pai1*Et[0]+pai2*Et[1])*sin(th);

	//y_1=(u0+pai1*Et[0]+pai2*Et[1])*sin(th)-((t-t0)/tE+pai2*Et[0]-pai1*Et[1])*cos(th);
	//y_2=(-u0-pai1*Et[0]-pai2*Et[1])*cos(th)-((t-t0)/tE+pai2*Et[0]-pai1*Et[1])*sin(th);

	while((Mag<0.9)&&(c<3)){
		Mag= BinaryMag(a,q,y_1,y_2,RSv,Tolv,&Images);
		Tolv/=10;
		c++;
		delete Images;
	}
	return Mag;
}

double VBBinaryLensing::BinaryLightCurveOrbital(double *pr,double t){
//(double a,double q,double u0,double th,double RSv,double tE,double t0,double t,double w1f,double w2,double w3v,double pai1,double pai2,double Tol){
	double Mag=-1.0,Tolv=Tol;
	double a=exp(pr[0]),q=exp(pr[1]),u0=pr[2],th=pr[3],RSv=exp(pr[4]),tE=exp(pr[5]),t0=pr[6],pai1=pr[7],pai2=pr[8],w1=pr[9],w2=pr[10],w3=pr[11];
	static double tEv,w,phi0,inc,phi,Cinc,Sinc,Cth,Sth,Cphi,Sphi,Cphi0,Sphi0,COm,SOm;
	static double w13,w123,den,den0,u,ty;
	static double Et[2];
	int c=0;
	_sols *Images;


	w13=w1*w1+w3*w3;
	w123=sqrt(w13+w2*w2);
	w13=sqrt(w13);
	if(w13>1.e-8){
		w3=(w3>1.e-8)? w3 : 1.e-8;
		w=w3*w123/w13;
		inc=acos(w2*w3/w13/w123);
		phi0=atan2(-w1*w123,w3*w13);
	}else{
		w=w2;
		inc=0.;
	    phi0=0.;
	}
	Cphi0=cos(phi0);
	Sphi0=sin(phi0);
	Cinc=cos(inc);
	Sinc=sin(inc);
	den0=sqrt(Cphi0*Cphi0+Cinc*Cinc*Sphi0*Sphi0);
	phi=(t-t0)*w+phi0;
	Cphi=cos(phi);
	Sphi=sin(phi);
	Cth=cos(th);
	Sth=sin(th);
	COm=(Cphi0*Cth-Cinc*Sth*Sphi0)/den0;
	SOm=(Cphi0*Sth+Cinc*Cth*Sphi0)/den0;
	tEv=tE/(1-0*tE*u0*w2);	
	den=sqrt(Cphi*Cphi+Cinc*Cinc*Sphi*Sphi);
	av=a*den/den0;

	ComputeParallax(t,t0,Et);

	u=u0+pai1*Et[1]-pai2*Et[0];
	ty=(t-t0)/tEv+pai1*Et[0]+pai2*Et[1];
	y_1=(Cphi*(u*SOm-ty*COm)-Cinc*Sphi*(u*COm+ty*SOm))/den;
	y_2=(-Cphi*(u*COm+ty*SOm)+Cinc*Sphi*(ty*COm-u*SOm))/den;


	if(y_1*y_1+y_2*y_2<100.){
		while((Mag<0.9)&&(c<3)){
			Mag= BinaryMag(av,q,y_1,y_2,RSv,Tolv,&Images);
			Tolv/=10;
			c++;
			delete Images;
		}
	}else{
		Mag=-Mag;
	}

	return Mag;
}

double VBBinaryLensing::BinaryMagDark(double a,double q,double y1,double y2,double RSv,double a1, double Tol){
	double Mag=-1.0,Magold=0.,Tolv=Tol;
//	static double a1=0.51;
	double tc,lb,rb,lc,rc,cb,cc,r2,cr2,scr2;
	int c=0,flag;
	double currerr,maxerr;
	annulus *first,*scan,*scan2;
	int nannold;
//		static double Mags[256];	
	_sols *Images;

	y_1=y1;
	y_2=y2;
//	if(fabs(t)<10.){
		while((Mag<0.9)&&(c<3)){
			
			first = new annulus;
			first->bin=0.;
			first->cum=0.;
			first->Mag=BinaryMag0(a,q,y_1,y_2,&Images);
			first->nim=Images->length;
			delete Images;
			first->f=3/(3-a1);
			first->err=0;
			first->prev=0;


			first->next= new annulus;
			scan=first->next;
			scan->prev=first;
			scan->next=0;
			scan->bin=1.;
			scan->cum=1.;
			scan->Mag=BinaryMag(a,q,y_1,y_2,RSv,Tolv,&Images);
			phi0=phip+1;
			scan->nim=Images->length;
			delete Images;
			scan->f=3/(3-a1)*(1-a1);
			if(scan->nim==scan->prev->nim){
				scan->err=fabs((scan->Mag-scan->prev->Mag)*(scan->prev->f-scan->f)/4);
			}else{
				scan->err=fabs((scan->Mag)*(scan->prev->f-scan->f)/4);
			}

			Magold=Mag=scan->Mag;
//			scan->err+=scan->Mag*Tolv*0.25; //Impose calculation of intermediate annulus at mag>4. Why?
			currerr=scan->err;
			flag=0;
			nannuli=nannold=1;
			while(((flag<nannold+5)&&(currerr>Tolv))|| (nannuli<minannuli)){
				maxerr=0;
				for(scan2=first->next;scan2;scan2=scan2->next){
#ifdef _PRINT_ERRORS_DARK
						printf("\n%d %lf %le | %lf %le",nannuli,scan2->Mag,scan2->err,Mag,currerr);
#endif
					if(scan2->err>maxerr){
						maxerr=scan2->err;
						scan=scan2;
					}
				}

				nannuli++;
				Magold=Mag;
				Mag-=(scan->Mag*scan->bin*scan->bin-scan->prev->Mag*scan->prev->bin*scan->prev->bin)*(scan->cum-scan->prev->cum)/(scan->bin*scan->bin-scan->prev->bin*scan->prev->bin);
				currerr-=scan->err;
				rb=scan->bin;
				rc=scan->cum;
				lb=scan->prev->bin;
				lc=scan->prev->cum;
				tc=(lc+rc)/2;
				do{
					cb=rb+(tc-rc)*(rb-lb)/(rc-lc);
					r2=cb*cb;
					cr2=1-r2;
					scr2=sqrt(cr2);
					cc=(3*r2*(1-a1)-2*a1*(scr2*cr2-1))/(3-a1);
					if(cc>tc){
						rb=cb;
						rc=cc;
					}else{
						lb=cb;
						lc=cc;
					}
				}while(fabs(cc-tc)>1.e-5);
				scan->prev->next=new annulus;
				scan->prev->next->prev=scan->prev;
				scan->prev=scan->prev->next;
				scan->prev->next=scan;
				scan->prev->bin=cb;
				scan->prev->cum=cc;
				scan->prev->f=3/(3-a1)*(1-a1*(1-scr2));
				scan->prev->Mag=BinaryMag(a,q,y_1,y_2,RSv*cb,Tolv,&Images);
				phi0+=phip;
				scan->prev->nim=Images->length;
				if(scan->prev->prev->nim==scan->prev->nim){				
					scan->prev->err=fabs((scan->prev->Mag-scan->prev->prev->Mag)*(scan->prev->prev->f-scan->prev->f)*(scan->prev->bin*scan->prev->bin-scan->prev->prev->bin*scan->prev->prev->bin)/4);
				}else{
					scan->prev->err=fabs((scan->prev->bin*scan->prev->bin*scan->prev->Mag-scan->prev->prev->bin*scan->prev->prev->bin*scan->prev->prev->Mag)*(scan->prev->prev->f-scan->prev->f)/4);
				}
				if(scan->nim==scan->prev->nim){
					scan->err=fabs((scan->Mag-scan->prev->Mag)*(scan->prev->f-scan->f)*(scan->bin*scan->bin-scan->prev->bin*scan->prev->bin)/4);
				}else{
					scan->err=fabs((scan->bin*scan->bin*scan->Mag-scan->prev->bin*scan->prev->bin*scan->prev->Mag)*(scan->prev->f-scan->f)/4);
				}
				rb=(scan->Mag+scan->prev->prev->Mag-2*scan->prev->Mag);
				scan->prev->err+=fabs(rb*(scan->prev->prev->f-scan->prev->f)*(scan->prev->bin*scan->prev->bin-scan->prev->prev->bin*scan->prev->prev->bin));
				scan->err+=fabs(rb*(scan->prev->f-scan->f)*(scan->bin*scan->bin-scan->prev->bin*scan->prev->bin));
#ifdef _PRINT_ERRORS_DARK
						printf("\n%d",Images->length);
#endif
				delete Images;
				
				Mag+=(scan->bin*scan->bin*scan->Mag-cb*cb*scan->prev->Mag)*(scan->cum-scan->prev->cum)/(scan->bin*scan->bin-scan->prev->bin*scan->prev->bin);
				Mag+=(cb*cb*scan->prev->Mag-scan->prev->prev->bin*scan->prev->prev->bin*scan->prev->prev->Mag)*(scan->prev->cum-scan->prev->prev->cum)/(scan->prev->bin*scan->prev->bin-scan->prev->prev->bin*scan->prev->prev->bin);
				currerr+=scan->err + scan->prev->err;

				if(fabs(Magold-Mag)*2<Tolv){
					flag++;
				}else{
					flag=0;
					nannold=nannuli;
				}
				
			}

			while(first){
				scan=first->next;
				delete first;
				first=scan;
			}

			//for(int i=0; i<nb;i++){
			//	Magp=Bins[i]*Bins[i]*BinaryMag(a,q,y_1,y_2,RSv*Bins[i],Tolv,&Images);
			//		Mags[i]=Magp;
			//	Mag+=Magp*Fl[i];
			//	if(i<nb-1) Mag-=Magp*Fl[i+1];
			//	delete Images;
			//}
			Tolv/=10;
			c++;
		}
	//}else{
	//	Mag=-Mag;
	//}

	return Mag;
}


double VBBinaryLensing::BinSourceMag(double *pr,double t){
//(double tE, double FR,double u1,double u2,double t01,double t02,double t){

	double tE=exp(pr[0]),FR=exp(pr[1]),u1=pr[2],u2=pr[3],t01=pr[4],t02=pr[5];

	y_1=u1*u1+(t-t01)*(t-t01)/tE/tE;
	y_2=u2*u2+(t-t02)*(t-t02)/tE/tE;

	return ((y_1+2)/sqrt(y_1*(y_1+4))+FR*(y_2+2)/sqrt(y_2*(y_2+4)))/(1+FR);
}

double VBBinaryLensing::BinSourceParallaxMag(double *pr, double t){
//(double tE, double FR,double u1,double u2,double t01,double t02,double t,double pai1,double pai2){
	double Mag=-1.0;
	double tE=exp(pr[0]),FR=exp(pr[1]),u1=pr[2],u2=pr[3],t01=pr[4],t02=pr[5],pai1=pr[6],pai2=pr[7];
	static double EtE[2];

//	static double r,r0,ty,fphi,Dfphi,phir,phil,Tp,tp,t0;
	static double y1,y2;
	
	t0=t01;

	ComputeParallax(t,t0,EtE);


	y_1=((t-t02)/tE+pai1*EtE[0]+pai2*EtE[1]);
	y_2=(u2+pai1*EtE[1]-pai2*EtE[0]);
	y2=(y_1)*(y_1)+(y_2)*(y_2);

	y_1=((t-t01)/tE+pai1*EtE[0]+pai2*EtE[1]);
	y_2=(u1+pai1*EtE[1]-pai2*EtE[0]);
	y1=(y_1)*(y_1)+(y_2)*(y_2);

	return ((y1+2)/sqrt(y1*(y1+4))+FR*(y2+2)/sqrt(y2*(y2+4)))/(1+FR);
}

double VBBinaryLensing::BinSourceXallarapMag(double *pr,double t){
//(double tE, double FR,double u1,double u2,double t01,double t02,double t,double pai1,double pai2,double w1f,double w2,double w3v){
	double tE=exp(pr[0]),FR=exp(pr[1]),u1=pr[2],u2=pr[3],t01=pr[4],t02=pr[5],pai1=pr[6],pai2=pr[7],w1=pr[8],w2=pr[9],w3=pr[10];

	double Mag=-1.0;
	static double tEv,w,phi0,inc,phi,Cinc,Sinc,Cth,Sth,Cphi,Sphi,Cphi0,Sphi0,COm,SOm;
	static double w13,w123,den,den0,u,ty;
	static double Et[2];
	static double y1,y2,q,a,th,u0;
	int c=0;

	q=sqrt(sqrt(FR)); // Attenzione! q dovrebbe essere in realtà un parametro libero ma non troppo. Da rivedere
	a=sqrt((u1-u2)*(u1-u2)+(t01-t02)*(t01-t02)/(tE*tE));
	th=atan((u1-u2)*tE/(t01-t02));
	u0=(u1+u2*q)/(1+q);
	t0=(t01+t02*q)/(1+q);

	w13=w1*w1+w3*w3;
	w123=sqrt(w13+w2*w2);
	w13=sqrt(w13);
	if(w13>1.e-8){
		w3=(w3>1.e-8)? w3 : 1.e-8;
		w=w3*w123/w13;
		inc=acos(w2*w3/w13/w123);
		phi0=atan2(-w1*w123,w3*w13);
	}else{
		w=w2;
		inc=0.;
	    phi0=0.;
	}
	Cphi0=cos(phi0);
	Sphi0=sin(phi0);
	Cinc=cos(inc);
	Sinc=sin(inc);
	den0=sqrt(Cphi0*Cphi0+Cinc*Cinc*Sphi0*Sphi0);
	phi=(t-t0)*w+phi0;
	Cphi=cos(phi);
	Sphi=sin(phi);
	Cth=cos(th);
	Sth=sin(th);
	COm=(Cphi0*Cth-Cinc*Sth*Sphi0)/den0;
	SOm=(Cphi0*Sth+Cinc*Cth*Sphi0)/den0;
	tEv=tE;	
	den=sqrt(Cphi*Cphi+Cinc*Cinc*Sphi*Sphi);
	av=a*den/den0;

	ComputeParallax(t,t0,Et);

	u=u0+pai1*Et[1]-pai2*Et[0];
	ty=(t-t0)/tEv+pai1*Et[0]+pai2*Et[1];
	y_1=(Cphi*(u*SOm-ty*COm)-Cinc*Sphi*(u*COm+ty*SOm))/den;
	y_2=(-Cphi*(u*COm+ty*SOm)+Cinc*Sphi*(ty*COm-u*SOm))/den;

	y1=(y_1+av*q/(1+q));
	y1*=y1;
	y1+=y_2*y_2;

	y2=(y_1-av/(1+q));
	y2*=y2;
	y2+=y_2*y_2;

	return ((y1+2)/sqrt(y1*(y1+4))+FR*(y2+2)/sqrt(y2*(y2+4)))/(1+FR);
}


///////////////////////////////////////////////
///////////////////////////////////////////////
///////////////////////////////////////////////
///////////////////////////////////////////////
///////////////////////////////////////////////

_curve *VBBinaryLensing::NewImages(complex yi,complex  *coefs,_theta *theta){
	static complex  y,yc,z,zc,zr[5]={0.,0.,0.,0.,0.},dy,dz,dJ;
	static double dzmax,dlmax=1.0e-6,good[5];
	static int worst1,worst2,worst3,bad,f1;
	static double av=0.0,m1v=0.0,disim,disisso;
	static _curve *Prov;
	static _point *scan,*prin,*fifth,*left,*right,*center;

#ifdef _PRINT_TIMES
	static double tim0,tim1;
#endif

	y=yi+coefs[11];
	yc=conj(y);

	/* coefs[6]=a*a; coefs[7]=a*a*a; coefs[8]=m2*m2; coefs[9]=a*a*m2*m2; coefs[10]=a*m2; coefs[11]=a*m1; coefs[20]=a; coefs[21]=m1; coefs[22]=m2;*/

	coefs[0]=coefs[9]*y;
	coefs[1]=coefs[10]*(coefs[20]*(coefs[21]+y*(2*yc-coefs[20]))-2*y);
	coefs[2]=y*(1-coefs[7]*yc)-coefs[20]*(coefs[21]+2*y*yc*(1+coefs[22]))+coefs[6]*(yc*(coefs[21]-coefs[22])+y*(1+coefs[22]+yc*yc));
	coefs[3]=2*y*yc+coefs[7]*yc+coefs[6]*(yc*(2*y-yc)-coefs[21])-coefs[20]*(y+2*yc*(yc*y-coefs[22]));
	coefs[4]=yc*(2*coefs[20]+y);
	coefs[4]=yc*(coefs[4]-1)-coefs[20]*(coefs[4]-coefs[21]);
	coefs[5]=yc*(coefs[20]-yc);

	bad=1;
	dzmax=1.0e-12;
	disim=-1.;
	f1=0;
	while(bad){

#ifdef _PRINT_TIMES
		tim0=Environment::TickCount;
#endif
		zroots(coefs,5,zr,1,dzmax);

#ifdef _PRINT_TIMES
		tim1=Environment::TickCount;
		inc+=tim1-tim0;
#endif
		// apply lens equation to check if it is really solved
		for(int i=0;i<5;i++){
			z=zr[i];
			zc=conj(z);
			good[i]=abs(_LL); // Lens equation check
			switch(i){
				case 0:
					worst1=i;
					break;
				case 1:
					if(good[i]>good[worst1]){
						worst2=worst1;
						worst1=i;
					}else worst2=i;
					break;
				case 2:
					if(good[i]>good[worst1]){
						worst3=worst2;
						worst2=worst1;
						worst1=i;
					}else if(good[i]>good[worst2]){
						worst3=worst2;
						worst2=i;
					}else worst3=i;
					break;
				default:
					if(good[i]>good[worst1]){
						worst3=worst2;
						worst2=worst1;
						worst1=i;
					}else if(good[i]>good[worst2]){
						worst3=worst2;
						worst2=i;
					}else if(good[i]>good[worst3]){
						worst3=i;
					}
			}
			//if(good[i]>good[worst1]){
			//	worst2=worst1;
			//	worst1=i;
			//}else if((good[i]>good[worst2])&&(i>0)){
			//	worst2=i;
			//}
		}
		if((good[worst3]<dlmax)&&((good[worst1]<dlmax)||(good[worst2]>1.e2*good[worst3]))){
			bad=0;
		} else{
			if((disim>0)&&(good[worst3]/disim>0.99)){
				if(f1>1){
					bad=0;
				}else{
					dzmax/=10;
					f1++;
				}
			}else{
				disim=good[worst3];
				dzmax/=10;
			}
		}
	}
	Prov=new _curve;
	if(good[worst1]>dlmax){
		for(int i=0;i<5;i++){
			if((i!=worst1)&&(i!=worst2)){
				//if((i==worst3)&&(good[i]>dlmax)&&(good[worst2]>1.e2*good[worst3])){
				//	zr[i]=(coefs[21].re<coefs[22].re)? 0.5*coefs[20]+coefs[21]/(0.5*coefs[20]-yc-coefs[22]/coefs[20]) : -0.5*coefs[20]+coefs[22]/(-0.5*coefs[20]-yc+coefs[21]/coefs[20]);
				//}
				Prov->append(zr[i].re,zr[i].im);

				zc=complex(Prov->last->x1, - Prov->last->x2);
				z=_J1c;
				dJ=1-z*conj(z);
				dy=complex(-sin(theta->th),cos(theta->th))*coefs[23];
				dz=(dy-z*conj(dy))/dJ;
				z=conj(zc);
				Prov->last->x1-=coefs[11].re;
				Prov->last->d=dz;
				Prov->last->dJ=dJ.re;
				Prov->last->ds=(imag(dy*dz*dz*_J2)+coefs[23].re*coefs[23].re)/dJ.re;

				Prov->last->theta=theta;

			}
		}
		theta->errworst=abs(zr[worst1]-zr[worst2]);
//		printf("\nworst: %le",abs(zr[worst1]-zr[worst2]));
	}else{
		f1=0;
		for(int i=0;i<5;i++){
			Prov->append(zr[i].re,zr[i].im);

			zc=complex(Prov->last->x1, - Prov->last->x2);
			z=_J1c;
			dJ=1-z*conj(z);
			dy=complex(-sin(theta->th),cos(theta->th))*coefs[23];
			dz=(dy-z*conj(dy))/dJ;
			z=conj(zc);
			Prov->last->x1-=coefs[11].re;
			Prov->last->d=dz;
			Prov->last->dJ=dJ.re;
			Prov->last->ds=(imag(dy*dz*dz*_J2)+coefs[23].re*coefs[23].re)/dJ.re;

			Prov->last->theta=theta;

			if(fabs(dJ.re)<1.e-5) f1=1;
		}
		theta->errworst=-1.e100;
// check degli Jacobiani in casi dubbi
		if(f1){
			left=right=center=fifth=0;
			dJ.re=0;
			for(scan=Prov->first;scan;scan=scan->next){
				if(_sign(scan->x2)==_sign(y.im)){
					prin=scan;
				}else{
					dz.re=fabs(scan->dJ);
					if(dz.re>dJ.re){ 
						fifth=scan;
						dJ.re=dz.re;
					}
				}
			}
			for(scan=Prov->first;scan;scan=scan->next){
				if((scan!=prin)&&(scan!=fifth)){
					if(left){
						if(scan->x1<left->x1){
							if(left!=right){
								center=left;
							}
							left=scan;
						}else{
							if(scan->x1>right->x1){
								if(left!=right){
									center=right;
								}
								right=scan;
							}else{
								center=scan;
							}
						}
					}else{
						left=right=center=scan;
					}
				}
			}
			if(left->dJ>0) left->dJ=-left->dJ;
			if(center->dJ<0) center->dJ=-center->dJ;
			if(right->dJ>0) right->dJ=-right->dJ;
		}
	}
	return Prov;
}

double VBBinaryLensing::BinaryMag0(double a1,double q1,double y1v,double y2v){
	_sols *images;
	double mag;
	mag=BinaryMag0(a1,q1,y1v,y2v,&images);
	delete images;
	return mag;
}

double VBBinaryLensing::BinaryMag0(double a1,double q1,double y1v,double y2v,_sols **Images){
	static complex a,q,m1,m2,y;
	static double av=-1.0,qv=-1.0;
	static complex  coefs[24],d1,d2,dy,dJ,dz;
	double Mag=-1.0;
	_theta *stheta;
	_curve *Prov, *Prov2;
	_point *scan1,*scan2;

	stheta=new _theta(0.);
	if((a1!=av)||(q1!=qv)){
		av=a1;
		qv=q1;
		if(q1<1){
			a=complex(-a1,0);
			q=complex(q1,0);
		}else{
			a=complex(a1,0);
			q=complex(1/q1,0);
		}
		m1=1.0/(1.0+q);
		m2=q*m1;

		coefs[20]=a;
		coefs[21]=m1;
		coefs[22]=m2;
		coefs[6]=a*a; 
		coefs[7]=coefs[6]*a; 
		coefs[8]=m2*m2; 
		coefs[9]=coefs[6]*coefs[8]; 
		coefs[10]=a*m2;
		coefs[11]=a*m1;
		coefs[23]=0;

	}
	y=complex (y1v,y2v);
	(*Images)=new _sols;
	Prov=NewImages(y,coefs,stheta);
	Mag=0.;
	for(scan1=Prov->first;scan1;scan1=scan2){
		scan2=scan1->next;
		Prov2= new _curve(scan1);
		Prov2->append(scan1->x1,scan1->x2);
		Prov2->last->theta=stheta;
		Prov2->last->d=Prov2->first->d;
		Prov2->last->dJ=Prov2->first->dJ;
		Prov2->last->ds=Prov2->first->ds;
		(*Images)->append(Prov2);
		Mag+=fabs(1/Prov2->last->dJ);
	}
	Prov->length=0;
	delete Prov;
	delete stheta;

	phip=1;
	return Mag;

}

double VBBinaryLensing::BinaryMag(double a1,double q1,double y1v,double y2v,double RSv,double Tol){
	_sols *images;
	double mag;
	mag=BinaryMag(a1,q1,y1v,y2v,RSv,Tol,&images);
	delete images;
	return mag;
}


double VBBinaryLensing::BinaryMag(double a1,double q1,double y1v,double y2v,double RSv,double Tol,_sols **Images){
	static complex a,q,m1,m2,y0,y,yc,z,zc;
	static double av=-1.0,qv=-1.0;
	static complex coefs[24],d1,d2,dy,dJ,dz;
	static double thoff=0.01020304;
	static double Mag=-1.0,th;
	static double errimage,maxerr,currerr,Magold;
	static int NPS,NPSmax,flag,NPSold;
	static _curve *Prov,*Prov2;
	static _point *scan1,*scan2;
	static _thetas *Thetas;
	static _theta *stheta,*itheta;

#ifdef _PRINT_TIMES
	static double tim0,tim1;
#endif


// Initialization of the equation coefficients

	if((a1!=av)||(q1!=qv)){
		av=a1;
		qv=q1;
		if(q1<1){
			a=complex(-a1,0);
			q=complex(q1,0);
		}else{
			a=complex(a1,0);
			q=complex(1/q1,0);
		}
		m1=1.0/(1.0+q);
		m2=q*m1;

		coefs[20]=a;
		coefs[21]=m1;
		coefs[22]=m2;
		coefs[6]=a*a; 
		coefs[7]=coefs[6]*a; 
		coefs[8]=m2*m2; 
		coefs[9]=coefs[6]*coefs[8]; 
		coefs[10]=a*m2;
		coefs[11]=a*m1;

	}
	coefs[23]=RSv;

	y0=complex (y1v,y2v);
	NPS=1;
	if(Tol>1.){
		errimage=0.;
		NPSmax=(int) (Tol);			
	}else{
		errimage=Tol*M_PI*RSv*RSv;
		NPSmax=32000;
	}

// Calculation of the images

	(*Images)=new _sols;
	Thetas=new _thetas;
	th=thoff;
	stheta=Thetas->insert(th);
	stheta->maxerr=1.e100;
	y=y0+complex(RSv*cos(thoff),RSv*sin(thoff));


#ifdef _PRINT_TIMES
		tim0=Environment::TickCount;
#endif
  	Prov=NewImages(y,coefs,stheta);
#ifdef _PRINT_TIMES
		tim1=Environment::TickCount;
		GM+=tim1-tim0;
#endif
	stheta=Thetas->insert(2.0*M_PI+thoff);
	stheta->maxerr=0.;
	stheta->Mag=0.;
	stheta->errworst=Thetas->first->errworst;
	for(scan1=Prov->first;scan1;scan1=scan2){
		scan2=scan1->next;
		Prov2= new _curve(scan1);
		Prov2->append(scan1->x1,scan1->x2);
		Prov2->last->theta=stheta;
		Prov2->last->d=Prov2->first->d;
		Prov2->last->dJ=Prov2->first->dJ;
		Prov2->last->ds=Prov2->first->ds;
		(*Images)->append(Prov2);
	}
	Prov->length=0;
	delete Prov;

	th=M_PI+thoff;
	flag=0;
	Magold=-1.;
	NPSold=2;
	do{
		//if(NPS==486){
		//	NPS=NPS;
		//}
		stheta=Thetas->insert(th);
		y=y0+complex(RSv*cos(th),RSv*sin(th));
#ifdef _PRINT_TIMES
		tim0=Environment::TickCount;
#endif
		Prov=NewImages(y,coefs,stheta);
#ifdef _PRINT_TIMES
		tim1=Environment::TickCount;
		GM+=tim1-tim0;
#endif
		OrderImages((*Images),Prov);
		if(stheta->th-stheta->prev->th<1.e-8){
			stheta->maxerr=0;
			stheta->prev->maxerr=0;
		}
		maxerr=currerr=Mag=0.;
		stheta=Thetas->first;
		while(stheta->next){
			currerr+=stheta->maxerr;
			Mag+=stheta->Mag;
#ifndef _uniform
			if(stheta->maxerr>maxerr){
				maxerr=stheta->maxerr;
#else
						if(stheta->next->th*0.99999-stheta->th>maxerr){
							maxerr=stheta->next->th-stheta->th;
#endif
				itheta=stheta;							
			}
			th=(itheta->th+itheta->next->th)/2;
			stheta=stheta->next;
#ifdef _selectimage
						if((NPS==NPSmax-1)&&(fabs(floor(stheta->th/M_PI*_npoints/2+0.5)-stheta->th/M_PI*_npoints/2)<1.e-8)){
							printf("%d %.15le\n",(int) floor(stheta->th/M_PI*_npoints/2+0.5),Mag);
						}
#endif
		}
		NPS++;
#ifndef _uniform
		if(fabs(Magold-Mag)*2<errimage){
			flag++;
		}else{
			flag=0;
			Magold=Mag;
			NPSold=NPS+1;
		}
#else
		currerr=2*errimage;
		if(NPS==2*NPSold){
			if(fabs(Magold-Mag)*2<errimage){
				flag=NPSold;
			}else{
				flag=0;
				NPSold=NPS;
				Magold=Mag;
			}
		}
#endif
#ifdef _PRINT_ERRORS2
		printf("\nNPS= %d Mag = %lf maxerr= %lg currerr =%lg th = %lf",NPS,Mag/(M_PI*RSv*RSv),maxerr/(M_PI*RSv*RSv),currerr/(M_PI*RSv*RSv),th);
#endif
	}while((currerr>errimage)&&(NPS<NPSmax)&&((flag<NPSold)/* || NPS<8 ||(currerr>10*errimage)*/)/*&&(flagits)*/);
	Mag/=(M_PI*RSv*RSv);
	phi=currerr/(M_PI*RSv*RSv);
	phip=NPS;

	delete Thetas;


	return Mag;
}

void VBBinaryLensing::OrderImages(_sols *Sols,_curve *Newpts){
	static double A[5][5];
	static _curve *cprec[5];
	static _curve *cpres[5];
	static _curve *cfoll[5];
	static _point *scan,*scan2,*scan3,*isso[2];
	static _curve *scurve,*scurve2;

	_theta *theta;
	double th,mi,cmp,cmp2;
	int nprec=0,npres,nfoll=0,issoc[2],ij;
	
	theta=Newpts->first->theta;
	th=theta->th;
	theta->Mag=theta->prev->Mag=theta->maxerr=theta->prev->maxerr=0;
	if(Newpts->length==3){
		mi=theta->next->errworst-theta->errworst;
		if((mi>theta->errworst)&&(theta->prev->errworst>0.)){
			theta->prev->maxerr=mi*mi;
		}
		mi=theta->prev->errworst-theta->errworst;
		if((mi>theta->errworst)&&(theta->next->errworst>0.)){
			theta->maxerr=mi*mi;
		}
	}

// Per ciascuna immagine troviamo il punto in cui inserire i nuovi punti
	scurve=Sols->first;
	for(int i=0;i<Sols->length;i++){
		if(th<scurve->first->theta->th){
			if(th>scurve->first->theta->prev->prev->th){
				cfoll[nfoll]=scurve; // immagine coinvolta all'inizio
				nfoll++;
				scurve2=scurve->next;
				Sols->drop(scurve);
				i--;
				scurve=scurve2;
			}else{
				scurve=scurve->next;
			}
		}else{
			if(th>scurve->last->theta->th){
				if(th<scurve->last->theta->next->next->th){
					cprec[nprec]=scurve; // immagine coinvolta alla fine
					nprec++;
				}
			}else{
                 // immagine coinvolta al centro
				scan=scurve->last;
				while(scan->theta->th>th){
					scan=scan->prev;
				}
				cfoll[nfoll]=scurve->divide(scan);
				nfoll++;
				cprec[nprec]=scurve;
				nprec++;
			}
			scurve=scurve->next;
		}
	}
	npres=Newpts->length;

			//if((theta->th>4.7116917419)&&(theta->th<4.711691759)){
			//	theta->th=theta->th;
			//}
			//if((theta->prev->th>4.7116917419)&&(theta->prev->th<4.711691759)){
			//	theta->th=theta->th;
			//}


// Caso di creazione nuove immagini
	if(nprec<npres){
		mi=100.;
		scan=Newpts->first;
		for(int i=0;i<Newpts->length-1;i++){
			scan2=scan->next;
			for(int j=i+1;j<Newpts->length;j++){
				cmp=(*scan2)- (*scan);
				if(cmp<mi){
					mi=cmp;
					isso[0]=scan;
					isso[1]=scan2;
				}
				scan2=scan2->next;
			}
			scan=scan->next;
		}
		Newpts->drop(isso[0]);
		Newpts->drop(isso[1]);
		scurve=new _curve(isso[0]);
		isso[0]->prev=isso[0]->next=0;
		scurve2=new _curve(isso[1]);
		isso[1]->prev=isso[1]->next=0;
		scurve->partneratstart=scurve2;
		scurve2->partneratstart=scurve;
		Sols->append(scurve);
		Sols->append(scurve2);
		cpres[3]=scurve;
		cpres[4]=scurve2;
		scan=isso[0];
		scan2=isso[1];

		cmp2=fabs(scan->d.re*scan2->d.re+scan->d.im*scan2->d.im);
		cmp=sqrt(mi/cmp2);
		mi=cmp*cmp*cmp/24;
		scurve->parabstart=-(-scan->ds+scan2->ds)*mi;

#ifdef _PRINT_ERRORS
		printf("\n%le %le %le %le %le %le %le %le",scan->x1,scan->x2,scan->dJ,(scan->x2+scan2->x2)*(scan2->x1-scan->x1)/2,scurve->parabstart,(scan->ds+scan2->ds)*mi/2,fabs(scurve->parabstart)*(cmp*cmp)/10,1.5*fabs(((scan->d.re-scan2->d.re)*(scan->x1-scan2->x1)+(scan->d.im-scan2->d.im)*(scan->x2-scan2->x2))-2*cmp*cmp2)*cmp);
#endif

		mi=fabs((scan->ds+scan2->ds)*mi/2)+fabs(scurve->parabstart)*(cmp*cmp/10)+1.5*fabs(((scan->d.re-scan2->d.re)*(scan->x1-scan2->x1)+(scan->d.im-scan2->d.im)*(scan->x2-scan2->x2))-2*cmp*cmp2)*cmp;
#ifdef _noparab
		mi=fabs(scurve->parabstart)*2+0*fabs((scan->x2+scan2->x2)*(scan2->x1-scan->x1)/2*cmp*cmp/6);
		scurve->parabstart=0.;
#endif
#ifdef _selectimage
			if(_selectionimage)
#endif
		theta->prev->Mag-=((scan->dJ>0)? -1 : 1)*((scan->x2+scan2->x2)*(scan2->x1-scan->x1)/2+scurve->parabstart);
		theta->prev->maxerr+=mi;
		scurve2->parabstart=-scurve->parabstart;

	}
	
// Caso di distruzione immagini
	if(nprec>npres){
		mi=100.;
		for(int i=0;i<nprec-1;i++){
			for(int j=i+1;j<nprec;j++){
				cmp=*(cprec[i]->last) - *(cprec[j]->last);
				if(cmp<mi){
					mi=cmp;
					issoc[0]=i;
					issoc[1]=j;
				}
			}
		}
		cprec[issoc[0]]->partneratend=cprec[issoc[1]];
		cprec[issoc[1]]->partneratend=cprec[issoc[0]];

		scan=cprec[issoc[0]]->last;
		scan2=cprec[issoc[1]]->last;

		cmp2=fabs(scan->d.re*scan2->d.re+scan->d.im*scan2->d.im);
		cmp=sqrt(mi/cmp2);
		mi=cmp*cmp*cmp/24;
		scan->parab=-(scan->ds-scan2->ds)*mi;

#ifdef _PRINT_ERRORS
		printf("\n%le %le %le %le %le %le %le %le",scan->x1,scan->x2,scan->dJ,(scan->x2+scan2->x2)*(scan2->x1-scan->x1)/2,scan->parab,(scan->ds+scan2->ds)*mi/2,fabs(scan->parab)*(cmp*cmp)/10,1.5*fabs(((scan->d.re-scan2->d.re)*(scan->x1-scan2->x1)+(scan->d.im-scan2->d.im)*(scan->x2-scan2->x2))+2*cmp*cmp2)*cmp);
#endif

		mi=fabs((scan->ds+scan2->ds)*mi/2)+fabs(scan->parab)*(cmp*cmp/10)+1.5*fabs(((scan->d.re-scan2->d.re)*(scan->x1-scan2->x1)+(scan->d.im-scan2->d.im)*(scan->x2-scan2->x2))+2*cmp*cmp2)*cmp;
#ifdef _noparab
		mi=fabs(scan->parab)*2+0*fabs((scan->x2+scan2->x2)*(scan2->x1-scan->x1)/2*cmp*cmp/6);
		scan->parab=0.;
#endif
#ifdef _selectimage
			if(_selectionimage)
#endif
		theta->prev->Mag+=((scan->dJ>0)? -1 : 1)*((scan->x2+scan2->x2)*(scan2->x1-scan->x1)/2+scan->parab);
		theta->prev->maxerr+=mi;
		scan2->parab=-scan->parab;

		nprec-=2;
		ij=0;
		for(int i=0;i<nprec;i++){
			if(i==issoc[0]) ij++;
			if(i==issoc[1]-1) ij++;
			cprec[i]=cprec[i+ij];
		}
	}

// Costruzione matrice distanze con immagini precedenti
	mi=100.;
	for(int i=0;i<nprec;i++){
		cpres[i]=cprec[i];
		scan=Newpts->first;
		for(int j=0;j<nprec;j++){
			A[i][j]=*(cprec[i]->last) - *scan;
			if(A[i][j]<mi){
				mi=A[i][j];
				issoc[0]=i;
				issoc[1]=j;
				isso[1]=scan;
			}
			scan=scan->next;
		}
	}

//  Associazione con le immagini che precedono
	while(nprec){
		scan=cprec[issoc[0]]->last;
		scan2=isso[1];

		cmp2=mi/fabs(scan->d.re*scan2->d.re+scan->d.im*scan2->d.im);
		cmp=(scan->theta->th-scan2->theta->th);
		mi=cmp*cmp*cmp/24;
		scan->parab=(scan->ds+scan2->ds)*mi;

#ifdef _PRINT_ERRORS
		printf("\n%le %le %le %le %le %le %le %le",scan->x1,scan->x2,scan->dJ,(scan->x2+scan2->x2)*(scan2->x1-scan->x1)/2,scan->parab,(scan->ds-scan2->ds)*mi/2,fabs(scan->parab)*(cmp2)/10,fabs(scan->parab)*(1.5*fabs(cmp2/(cmp*cmp)-1)));
#endif

		mi=fabs((scan->ds-scan2->ds)*mi/2)+fabs(scan->parab*(cmp2/10+1.5*fabs(cmp2/(cmp*cmp)-1)));
#ifdef _noparab
		mi=fabs(scan->parab)*2+0*fabs((scan->x2+scan2->x2)*(scan2->x1-scan->x1)/2*cmp*cmp/6);
		scan->parab=0.;
#endif
#ifdef _selectimage
			if(_selectionimage)
#endif
		theta->prev->Mag+=((scan->dJ>0)? -1 : 1)*((scan->x2+scan2->x2)*(scan2->x1-scan->x1)/2+scan->parab);
		theta->prev->maxerr+=mi;

		Newpts->drop(isso[1]);
		cprec[issoc[0]]->append(isso[1]);
		cprec[issoc[0]]->partneratend=0;

		nprec--;
		for(int i=issoc[0];i<nprec;i++){
			cprec[i]=cprec[i+1];
			for(int j=0;j<nprec+1;j++){
				A[i][j]=A[i+1][j];
			}
		}
		for(int j=issoc[1];j<nprec;j++){
			for(int i=0;i<nprec;i++){
				A[i][j]=A[i][j+1];
			}
		}
		mi=100.;
		for(int i=0;i<nprec;i++){
			scan=Newpts->first;
			for(int j=0;j<nprec;j++){
				if(A[i][j]<mi){
					mi=A[i][j];
					issoc[0]=i;
					issoc[1]=j;
					isso[1]=scan;
				}
				scan=scan->next;
			}
		}
	}
	delete Newpts;

#ifdef _PRINT_ERRORS
		printf("\nN");
#endif

// immagini seguenti
	if(nfoll){
	// Caso di creazione nuove immagini
		if(npres<nfoll){
			mi=100.;
			for(int i=0;i<nfoll-1;i++){
				for(int j=i+1;j<nfoll;j++){
					cmp=*(cfoll[i]->first) - *(cfoll[j]->first);
					if(cmp<mi){
						mi=cmp;
						issoc[0]=i;
						issoc[1]=j;
					}
				}
			}
			cfoll[issoc[0]]->partneratstart=cfoll[issoc[1]];
			cfoll[issoc[1]]->partneratstart=cfoll[issoc[0]];

			scan=cfoll[issoc[0]]->first;
			scan2=cfoll[issoc[1]]->first;

			cmp2=fabs(scan->d.re*scan2->d.re+scan->d.im*scan2->d.im);
			cmp=sqrt(mi/cmp2);
			mi=cmp*cmp*cmp/24;
			cfoll[issoc[0]]->parabstart=(scan->ds-scan2->ds)*mi;

#ifdef _PRINT_ERRORS
			printf("\n%le %le %le %le %le %le %le %le",scan->x1,scan->x2,scan->dJ,(scan->x2+scan2->x2)*(scan2->x1-scan->x1)/2,cfoll[issoc[0]]->parabstart,(scan->ds+scan2->ds)*mi/2,fabs(cfoll[issoc[0]]->parabstart)*(cmp*cmp)/10,1.5*fabs(((scan->d.re-scan2->d.re)*(scan->x1-scan2->x1)+(scan->d.im-scan2->d.im)*(scan->x2-scan2->x2))-2*cmp*cmp2)*cmp);
#endif
			mi=fabs((scan->ds+scan2->ds)*mi/2)+fabs(cfoll[issoc[0]]->parabstart)*(cmp*cmp/10)+1.5*fabs(((scan->d.re-scan2->d.re)*(scan->x1-scan2->x1)+(scan->d.im-scan2->d.im)*(scan->x2-scan2->x2))-2*cmp*cmp2)*cmp;
#ifdef _noparab
			mi=fabs(cfoll[issoc[0]]->parabstart)*2+0*fabs((scan->x2+scan2->x2)*(scan2->x1-scan->x1)/2*cmp*cmp/6);
			cfoll[issoc[0]]->parabstart=0.;
#endif
#ifdef _selectimage
			if(_selectionimage)
#endif
			theta->Mag-=((scan->dJ>0)? -1 : 1)*((scan->x2+scan2->x2)*(scan2->x1-scan->x1)/2+cfoll[issoc[0]]->parabstart);
			theta->maxerr+=mi;
			cfoll[issoc[1]]->parabstart=-cfoll[issoc[0]]->parabstart;

			Sols->append(cfoll[issoc[0]]);
			Sols->append(cfoll[issoc[1]]);
			nfoll-=2;
			ij=0;
			for(int i=0;i<nfoll;i++){
				if(i==issoc[0]) ij++;
				if(i==issoc[1]-1) ij++;
				cfoll[i]=cfoll[i+ij];
			}
		}
		
		
	// Caso di distruzione immagini
		if(npres>nfoll){
			mi=100.;
			for(int i=0;i<npres-1;i++){
				for(int j=i+1;j<npres;j++){
					cmp=*(cpres[i]->last) - *(cpres[j]->last);
					if(cmp<mi){
						mi=cmp;
						issoc[0]=i;
						issoc[1]=j;
					}
				}
			}
			cpres[issoc[0]]->partneratend=cpres[issoc[1]];
			cpres[issoc[1]]->partneratend=cpres[issoc[0]];
			
			scan=cpres[issoc[0]]->last;
			scan2=cpres[issoc[1]]->last;

			cmp2=fabs(scan->d.re*scan2->d.re+scan->d.im*scan2->d.im);
			cmp=sqrt(mi/cmp2);
			mi=cmp*cmp*cmp/24;
			scan->parab=-(scan->ds-scan2->ds)*mi;
		
#ifdef _PRINT_ERRORS
			printf("\n%le %le %le %le %le %le %le %le",scan->x1,scan->x2,scan->dJ,(scan->x2+scan2->x2)*(scan2->x1-scan->x1)/2,scan->parab,(scan->ds+scan2->ds)*mi/2,fabs(scan->parab)*(cmp*cmp)/10,1.5*fabs(((scan->d.re-scan2->d.re)*(scan->x1-scan2->x1)+(scan->d.im-scan2->d.im)*(scan->x2-scan2->x2))+2*cmp*cmp2)*cmp);
#endif

			mi=fabs((scan->ds+scan2->ds)*mi/2)+fabs(scan->parab)*(cmp*cmp/10)+1.5*fabs(((scan->d.re-scan2->d.re)*(scan->x1-scan2->x1)+(scan->d.im-scan2->d.im)*(scan->x2-scan2->x2))+2*cmp*cmp2)*cmp;
#ifdef _noparab
			mi=fabs(scan->parab)*2+0*fabs((scan->x2+scan2->x2)*(scan2->x1-scan->x1)/2*cmp*cmp/6);
			scan->parab=0.;
#endif
#ifdef _selectimage
			if(_selectionimage)
#endif
			theta->Mag+=((scan->dJ>0)? -1 : 1)*((scan->x2+scan2->x2)*(scan2->x1-scan->x1)/2+scan->parab);
			theta->maxerr+=mi;
			scan2->parab=-scan->parab;

			npres-=2;
			ij=0;
			for(int i=0;i<npres;i++){
				if(i==issoc[0]) ij++;
				if(i==issoc[1]-1) ij++;
				cpres[i]=cpres[i+ij];
			}
		}

	// Costruzione matrice distanze con immagini seguenti
		mi=100.;
		for(int i=0;i<npres;i++){
			for(int j=0;j<npres;j++){
				A[i][j]=*(cpres[i]->last) - *(cfoll[j]->first);
				if(A[i][j]<mi){
					mi=A[i][j];
					issoc[0]=i;
					issoc[1]=j;
				}
			}
		}


	//  Associazione con le immagini che seguono
		while(npres){
			scan=cpres[issoc[0]]->last;
			scan2=cfoll[issoc[1]]->first;

			cmp2=mi/fabs(scan->d.re*scan2->d.re+scan->d.im*scan2->d.im);
			cmp=(scan->theta->th-scan2->theta->th);
			mi=cmp*cmp*cmp/24;
			scan->parab=(scan->ds+scan2->ds)*mi;

#ifdef _PRINT_ERRORS
			printf("\n%le %le %le %le %le %le %le %le",scan->x1,scan->x2,scan->dJ,(scan->x2+scan2->x2)*(scan2->x1-scan->x1)/2,scan->parab,(scan->ds-scan2->ds)*mi/2,fabs(scan->parab)*(cmp2)/10,fabs(scan->parab)*(1.5*fabs(cmp2/(cmp*cmp)-1)));
#endif

			mi=fabs((scan->ds-scan2->ds)*mi/2)+fabs(scan->parab*(cmp2/10+1.5*fabs(cmp2/(cmp*cmp)-1)));
#ifdef _noparab
			mi=fabs(scan->parab)*2+0*fabs((scan->x2+scan2->x2)*(scan2->x1-scan->x1)/2*cmp*cmp/6);
			scan->parab=0.;
#endif
#ifdef _selectimage
			if(_selectionimage)
#endif
			theta->Mag+=((scan->dJ>0)? -1 : 1)*((scan->x2+scan2->x2)*(scan2->x1-scan->x1)/2+scan->parab);

			theta->maxerr+=mi;

			cpres[issoc[0]]->join(cfoll[issoc[1]]);
	
			npres--;
			for(int i=issoc[0];i<npres;i++){
				cpres[i]=cpres[i+1];
				for(int j=0;j<npres+1;j++){
					A[i][j]=A[i+1][j];
				}
			}
			for(int j=issoc[1];j<npres;j++){
				cfoll[j]=cfoll[j+1];
				for(int i=0;i<npres;i++){
					A[i][j]=A[i][j+1];
				}
			}
			mi=100.;
			for(int i=0;i<npres;i++){
				for(int j=0;j<npres;j++){
					if(A[i][j]<mi){
						mi=A[i][j];
						issoc[0]=i;
						issoc[1]=j;
					}
				}
			}
		}
	}

}


_point::_point(double x,double y,_theta *theta1){
	x1=x;
	x2=y;
	theta=theta1;
}

double _point::operator-(_point p2){
	return (x1-p2.x1)*(x1-p2.x1)+(x2-p2.x2)*(x2-p2.x2);
}

_curve::_curve(void){
	length=0;
	first=last=0;
	partneratstart=partneratend=0;
}

_curve::_curve(_point *p1){
	length=1;
	first=last=p1;
	p1->prev=p1->next=0;
	partneratstart=partneratend=0;
}

_curve::~_curve(void){
	_point *scan1,*scan2;
	scan1=first;
	for(int i=0;i<length;i++){
		scan2=scan1->next;
		delete scan1;
		scan1=scan2;
	}
}

_curve *_curve::divide(_point *ref){
	_point *scan;
	_curve *nc;
	int l1;

	l1=1;
	for(scan=first;scan!=ref;scan=scan->next) l1++;
	nc=new _curve();
	nc->first=ref->next;
	nc->first->prev=0;
	nc->last=last;
	nc->length=length-l1;
	nc->partneratend=partneratend;
	if(partneratend) partneratend->partneratend=nc;

	length=l1;
	last=ref;
	ref->next=0;
	partneratend=0;
	return nc;
}


void _curve::append(double x1,double x2){
	_point *pp;
	pp=new _point(x1,x2,0);
	if(length==0){
		first=pp;
		last=pp;
		pp->prev=0;
	}else{
		last->next=pp;
		pp->prev=last;
		last=pp;
	}
	pp->next=0;
	length++;
}

void _curve::append(_point *pp){

	pp->next=last->next;
	pp->prev=last;
	last->next=pp;
	last=pp;
	length++;
}

void _curve::prepend(double x1,double x2){
	_point *pp;
	pp=new _point(x1,x2,0);
	if(length==0){
		first=pp;
		last=pp;
		pp->next=0;
	}else{
		first->prev=pp;
		pp->next=first;
		first=pp;
	}
	pp->prev=0;
	length++;
}

_curve *_curve::join(_curve *nc){
	if(length>0){
		last->next=nc->first;
	}else{
		first=nc->first;
	};
	if(nc->length>0){
		nc->first->prev=last;
		last=nc->last;
	}
	length+=nc->length;
	partneratend=nc->partneratend;
	if(partneratend) partneratend->partneratend=this;
	nc->first=0;
	nc->last=0;
	nc->length=0;
	delete nc;
	return this;
}

_curve *_curve::joinbefore(_curve *nc){
	if(length>0){
		first->prev=nc->last;
	}else{
		last=nc->last;
	};
	if(nc->length>0){
		nc->last->next=first;
		first=nc->first;
	}
	length+=nc->length;
	nc->first=0;
	nc->last=0;
	nc->length=0;
	delete nc;
	return this;
}

_curve *_curve::reverse(void){
	_point *scan1,*scan2,*scambio;
	if(length>1){
		scan1=first;
		while(scan1){
			scan2=scan1->next;
			scambio=scan1->next;
			scan1->next=scan1->prev;
			scan1->prev=scambio;
			scan1=scan2;
		}
		scambio=first;
		first=last;
		last=scambio;
	}
	return this;
}

void _curve::drop(_point *ref){
	_point *scan;
	if(length){
		for(scan=last;scan&&(scan!=ref);scan=scan->prev);
		if(scan){
			if(length==1){
				first=last=0;
			}else{
				if(ref->prev){
					ref->prev->next=ref->next;
					if(ref==last){
						last=ref->prev;
					}
				}
				if(ref->next){
					ref->next->prev=ref->prev;
					if(ref==first){
						first=ref->next;
					}
				}
			}
			length--;
		}
	}
}

double _curve::closest2(_point *ref,_point **clos2){
	double mi=100.0,mi2=100.0,FP;
	_point *scan,*clos;
	if(length>1){
		clos=*clos2=first;
		for(scan=first; scan!=0; scan=scan->next){
			FP=*scan - *ref;
			if(FP<mi){
				mi2=mi;
				mi=FP;
				*clos2=clos;
				clos=scan;
			}else if(FP<mi2){
				mi2=FP;
				*clos2=scan;
			}
		}
	}else{
		*clos2=0;
	}
	return (**clos2 - *ref);
}

double _curve::closest(_point *ref,_point **clos){
	double mi=100.0,FP;
	_point *scan;
	for(scan=first; scan!=0; scan=scan->next){
		FP=*scan - *ref;
		if(FP<mi){
			mi=FP;
			*clos=scan;
		}
	}
	return mi;
}

void _curve::complement(_point **sott,int lensott, _point **res, int lenres){
	int flag,i;
	_point *scan;
	i=0;
	for(scan=first; scan!=0;scan=scan->next){
		flag=0;
		for(int j=0;(j<lensott)&&(!flag);j++){
			if(scan==sott[j]){
				flag=1;
			}
		}
		if((!flag)&&(i<lenres)){
			res[i]=scan;
			i++;
		}
	}
}
//
//void _curve::sort(){
//	_point *scan,*pres,*futu;
//	double maxmag=-1.;
//	for(pres=first; pres!=0;pres=futu){
//		if(pres->x2>maxmag) maxmag=pres->x2;
//		futu=pres->next;
//		for(scan=first;(scan!=pres)&&(scan->prev!=pres);scan=scan->next){
//			if(scan->x1>pres->x1){
//				scan->prev->next=pres;
//				scan->prev=pres;
//				pres->prev->next=futu;
//				if(futu) futu->prev=pres->prev;
//				pres->prev=scan->prev;
//				pres->next=scan;
//			}
//		}
//	}
//	return maxmag;
//}


_sols::_sols(void){
	length=0;
	first=last=0;
}

_sols::~_sols(void){
	_curve *scan1,*scan2;
	scan1=first;
	while(scan1){
		scan2=scan1->next;
		delete scan1;
		scan1=scan2;
	}
}

void _sols::append(_curve *cc){
	if(length==0){
		first=cc;
		last=cc;
		cc->prev=0;
	}else{
		last->next=cc;
		cc->prev=last;
		last=cc;
	}
	cc->next=0;
	length++;
}

void _sols::prepend(_curve *cc){
	if(length==0){
		first=cc;
		last=cc;
		cc->next=0;
	}else{
		first->prev=cc;
		cc->next=first;
		first=cc;
	}
	cc->prev=0;
	length++;
}

void _sols::drop(_curve *ref){
	_curve *scan;
	if(length){
		for(scan=last;scan&&(scan!=ref);scan=scan->prev);
		if(scan){
			if(length==1){
				first=last=0;
			}else{
				if(ref->prev){
					ref->prev->next=ref->next;
					if(ref==last){
						last=ref->prev;
					}
				}
				if(ref->next){
					ref->next->prev=ref->prev;
					if(ref==first){
						first=ref->next;
					}
				}
			}
			length--;
		}
	}
}

void _sols::join(_sols *nc){
	if(length>0){
		last->next=nc->first;
	}else{
		first=nc->first;
	};
	if(nc->length>0){
		nc->first->prev=last;
		last=nc->last;
	}
	length+=nc->length;
	nc->first=0;
	nc->last=0;
	nc->length=0;
	delete nc;
}

_theta::_theta(double th1){
	th=th1;
}
_thetas::_thetas(void){
	length=0;
}

_thetas::~_thetas(void){
	_theta *scan,*scan2;
	scan=first;
	while(scan){
		scan2=scan->next;
		delete scan;
		scan=scan2;
	}
}

_theta *_thetas::insert(double th){
	_theta *scan,*scan2;

	scan2=new _theta(th);
	if(length){
		if(th<first->th){
			first->prev=scan2;
			scan2->next=first;
			scan2->prev=0;
			first=scan2;
		}else{
			if(th>last->th){
				last->next=scan2;
				scan2->prev=last;
				scan2->next=0;
				last=scan2;
			}else{
				scan=first;
				while(scan->th<th) scan=scan->next;
				scan2->next=scan;
				scan2->prev=scan->prev;
				scan->prev->next=scan2;
				scan->prev=scan2;
			}
		}
	}else{
		first=scan2;
		last=scan2;
		scan2->next=0;
		scan2->prev=0;
	}
	length++;
//	scan2->maxerr=0.;
	return scan2;
}
complex::complex(double a,double b){
	re=a;
	im=b;
}

complex::complex(double a){
	re=a;
	im=0;
}

complex::complex(void){
	re=0;
	im=0;
}

double abs(complex z){
	return sqrt(z.re*z.re+z.im*z.im);
}

complex conj(complex z){
	return complex(z.re,-z.im);
}

complex sqrt(complex z){
	double md=sqrt(z.re*z.re+z.im*z.im);	
	return (md>0)? sqrt(md)*complex((sqrt((1+z.re/md)/2)*((z.im>0) ? 1 : -1)),sqrt((1-z.re/md)/2)) : 0.0;
}

double real(complex z){
	return z.re;
}

double imag(complex z){
	return z.im;
}

complex operator+(complex p1,complex p2){
	return complex(p1.re+p2.re,p1.im+p2.im);
}

complex operator-(complex p1,complex p2){
	return complex(p1.re-p2.re,p1.im-p2.im);
}

complex operator*(complex p1,complex p2){
	return complex(p1.re*p2.re-p1.im*p2.im,p1.re*p2.im+p1.im*p2.re);
}

complex operator/(complex p1,complex p2){
	double md=p2.re*p2.re+p2.im*p2.im;
	return complex((p1.re*p2.re+p1.im*p2.im)/md,(p1.im*p2.re-p1.re*p2.im)/md);
}

complex operator+(complex z,double a){
	return complex(z.re+a,z.im);
}

complex operator-(complex z,double a){
	return complex(z.re-a,z.im);
}

complex operator*(complex z,double a){
	return complex(z.re*a,z.im*a);
}

complex operator/(complex z,double a){
	return complex(z.re/a,z.im/a);
}

complex operator+(double a,complex z){
	return complex(z.re+a,z.im);
}

complex operator-(double a,complex z){
	return complex(a-z.re,-z.im);
}

complex operator*(double a,complex z){
	return complex(a*z.re,a*z.im);
}

complex operator/(double a,complex z){
	double md=z.re*z.re+z.im*z.im;
	return complex(a*z.re/md,-a*z.im/md);
}


complex operator+(complex z,int a){
	return complex(z.re+a,z.im);
}

complex operator-(complex z,int a){
	return complex(z.re-a,z.im);
}

complex operator*(complex z,int a){
	return complex(z.re*a,z.im*a);
}

complex operator/(complex z,int a){
	return complex(z.re/a,z.im/a);
}

complex operator+(int a,complex z){
	return complex(z.re+a,z.im);
}

complex operator-(int a,complex z){
	return complex(a-z.re,-z.im);
}

complex operator*(int a,complex z){
	return complex(a*z.re,a*z.im);
}

complex operator/(int a,complex z){
	double md=z.re*z.re+z.im*z.im;
	return complex(a*z.re/md,-a*z.im/md);
}

complex operator-(complex z){
	return complex(-z.re,-z.im);
}


void VBBinaryLensing::laguer(complex *a, int m, complex *x, int *its,double dzmax){
	static int iter,j;
	static double abx,abp,abm,err;
	static complex dx,x1,b,d,f,g,h,sq,gp,gm,g2;
	static double frac[MR+1]={0.0,0.5,0.25,0.75,0.13,0.38,0.62,0.88,1.0};

	for(iter=1;iter<=MAXIT;iter++){
		*its=iter;
		b=a[m];
		err=abs(b);
		d=f=complex (0.0,0.0);
		abx=abs(*x);
		for(j=m-1;j>=0;j--){
			f=(*x)*f+d;
			d=(*x)*d+b;
			b=(*x)*b+a[j];
			err=abs(b)+abx*err;
		}
		err*=dzmax;
		if(abs(b)<=err) return;
		g=d/b;
		g2=g*g;
		h=g2-2.0*f/b;
		sq=sqrt((m-1)*(m*h-g2));
		gp=g+sq;
		gm=g-sq;
		abp=abs(gp);
		abm=abs(gm);
		if(abp<abm){
			gp=gm;
			abp=abm;
		}
		dx=(abp >0.0 ? (m/gp) :  complex ((1+abx)*cos((double) iter),(1+abx)*sin((double) iter)));
		x1=*x - dx;
		if((x->re == x1.re) && (x->im == x1.im)){ 
			flagits++; 
			return;
		}
		if(iter % MT) *x=x1;
		else *x=*x- (frac[iter/MT])*dx;
	}
	//printf("too many iterations in laguer\n");
	//printf("%lg",dzmax);
//	printf("y_1=%lg y_2=%lg e=%lg a=%lg av=%lg\nphi=%lg phip=%lg phi0=%lg Om=%lf inc=%lg\n t0=%lg d3=%lg v3=%lg GM=%lg\n",y_1,y_2,e,a,av,phi,phip,phi0,Om,inc,t0,d3,v3,GM);
//	getchar();
	return;
	}

void VBBinaryLensing::zroots(complex *a,int m,complex *roots, int polish,double dzmax){
	static int its,j,jj;
	static complex x,b,c,ad[MAXM];

	flagits=-m;
	for(j=0;j<=m;j++) ad[j]=a[j];
	for(j=m;j>=1;j--){
		x=roots[j-1];
		laguer(ad,j,&x,&its,dzmax);
		roots[j-1]=x;
		b=ad[j];
		for(jj=j-1;jj>=0;jj--){
			c=ad[jj];
			ad[jj]=b;
			b=x*b+c;
		}
	}
	if(polish)
		for(j=0;j<m;j++){
			laguer(a,m,&roots[j],&its,dzmax);
		}
}


#define _ERRTOL 0.01
#define _TINY 1.5e-38
#define _BIG 3.0e37
#define _THIRD (1.0/3.0)
#define _C1 (1.0/24.0)
#define _C2 0.1
#define _C3 (3.0/44.0)
#define _C4 (1.0/14.0)

double FMIN(double x,double y){
	return (x<y)? x : y;
}

double FMAX(double x,double y){
	return (x>y)? x : y;
}

double SQR(double x){
	return x*x;
}

double rf(double x, double y, double z)
//Computes Carlson?s elliptic integral of the first kind, RF (x, y, z). x, y, and z must be nonnegative,
//and at most one can be zero. TINY must be at least 5 times the machine underflow limit,
//BIG at most one fifth the machine overflow limit.
{
	double alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,sqrtz,xt,yt,zt;
//	if (FMIN(FMIN(x,y),z) < 0.0 || FMIN(FMIN(x+y,x+z),y+z) < TINY || FMAX(FMAX(x,y),z) > BIG) nrerror("invalid arguments in rf");
	xt=x;
	yt=y;
	zt=z;
	do {
		sqrtx=sqrt(xt);
		sqrty=sqrt(yt);
		sqrtz=sqrt(zt);
		alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz;
		xt=0.25*(xt+alamb);
		yt=0.25*(yt+alamb);
		zt=0.25*(zt+alamb);
		ave=_THIRD*(xt+yt+zt);
		delx=(ave-xt)/ave;
		dely=(ave-yt)/ave;
		delz=(ave-zt)/ave;
	} while (FMAX(FMAX(fabs(delx),fabs(dely)),fabs(delz)) > _ERRTOL);
	e2=delx*dely-delz*delz;
	e3=delx*dely*delz;
	return (1.0+(_C1*e2-_C2-_C3*e3)*e2+_C4*e3)/sqrt(ave);
}

#define __TINY 1.0e-25
#define __BIG 4.5e21
#define __C1 (3.0/14.0)
#define __C2 (1.0/6.0)
#define __C3 (9.0/22.0)
#define __C4 (3.0/26.0)
#define __C5 (0.25*__C3)
#define __C6 (1.5*__C4)

double rd(double x, double y, double z)
//Computes Carlson?s elliptic integral of the second kind, RD(x, y, z). x and y must be nonnegative,
//and at most one can be zero. z must be positive. TINY must be at least twice the
//negative 2/3 power of the machine overflow limit. BIG must be at most 0.1 × ERRTOL times
//the negative 2/3 power of the machine underflow limit.
{
	double alamb,ave,delx,dely,delz,ea,eb,ec,ed,ee,fac,sqrtx,sqrty,
	sqrtz,sum,xt,yt,zt;
//	if (FMIN(x,y) < 0.0 || FMIN(x+y,z) < TINY || FMAX(FMAX(x,y),z) > BIG) nrerror("invalid arguments in rd");
	xt=x;
	yt=y;
	zt=z;
	sum=0.0;
	fac=1.0;
	do {
		sqrtx=sqrt(xt);
		sqrty=sqrt(yt);
		sqrtz=sqrt(zt);
		alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz;
		sum += fac/(sqrtz*(zt+alamb));
		fac=0.25*fac;
		xt=0.25*(xt+alamb);
		yt=0.25*(yt+alamb);
		zt=0.25*(zt+alamb);
		ave=0.2*(xt+yt+3.0*zt);
		delx=(ave-xt)/ave;
		dely=(ave-yt)/ave;
		delz=(ave-zt)/ave;
	} while (FMAX(FMAX(fabs(delx),fabs(dely)),fabs(delz)) > _ERRTOL);
	ea=delx*dely;
	eb=delz*delz;
	ec=ea-eb;
	ed=ea-6.0*eb;
	ee=ed+ec+ec;
	return 3.0*sum+fac*(1.0+ed*(-__C1+__C5*ed-__C6*delz*ee) +delz*(__C2*ee+delz*(-__C3*ec+delz*__C4*ea)))/(ave*sqrt(ave));
}

#define ___C1 (3.0/14.0)
#define ___C2 (1.0/3.0)
#define ___C3 (3.0/22.0)
#define ___C4 (3.0/26.0)
#define ___C5 (0.75*___C3)
#define ___C6 (1.5*___C4)
#define ___C7 (0.5*___C2)
#define ___C8 (___C3+___C3)

double rj(double x, double y, double z, double p)
//Computes Carlson?s elliptic integral of the third kind, RJ (x, y, z, p). x, y, and z must be
//nonnegative, and at most one can be zero. p must be nonzero. If p < 0, the Cauchy principal
//value is returned. TINY must be at least twice the cube root of the machine underflow limit,
//BIG at most one fifth the cube root of the machine overflow limit.
{
	double a,alamb,alpha,ans,ave,b,beta,delp,delx,dely,delz,ea,eb,ec,
	ed,ee,fac,pt,rcx,rho,sqrtx,sqrty,sqrtz,sum,tau,xt,yt,zt;
//	if (FMIN(FMIN(x,y),z) < 0.0 || FMIN(FMIN(FMIN(x+y,x+z),y+z),fabs(p)) < TINY || FMAX(FMAX(FMAX(x,y),z),fabs(p)) > BIG)	nrerror("invalid arguments in rj");
	sum=0.0;
	fac=1.0;
	if (p > 0.0) {
		xt=x;
		yt=y;
		zt=z;
		pt=p;
	} else {
		xt=FMIN(FMIN(x,y),z);
		zt=FMAX(FMAX(x,y),z);
		yt=x+y+z-xt-zt;
		a=1.0/(yt-p);
		b=a*(zt-yt)*(yt-xt);
		pt=yt+b;
		rho=xt*zt/yt;
		tau=p*pt/yt;
		rcx=rc(rho,tau);
	}
	do {
		sqrtx=sqrt(xt);
		sqrty=sqrt(yt);
		sqrtz=sqrt(zt);
		alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz;
		alpha=SQR(pt*(sqrtx+sqrty+sqrtz)+sqrtx*sqrty*sqrtz);
		beta=pt*SQR(pt+alamb);
		sum += fac*rc(alpha,beta);
		fac=0.25*fac;
		xt=0.25*(xt+alamb);
		yt=0.25*(yt+alamb);
		zt=0.25*(zt+alamb);
		pt=0.25*(pt+alamb);
		ave=0.2*(xt+yt+zt+pt+pt);
		delx=(ave-xt)/ave;
		dely=(ave-yt)/ave;
		delz=(ave-zt)/ave;
		delp=(ave-pt)/ave;
	} while (FMAX(FMAX(FMAX(fabs(delx),fabs(dely)), fabs(delz)),fabs(delp)) > _ERRTOL);
	ea=delx*(dely+delz)+dely*delz;
	eb=delx*dely*delz;
	ec=delp*delp;
	ed=ea-3.0*ec;
	ee=eb+2.0*delp*(ea-ec);
	ans=3.0*sum+fac*(1.0+ed*(-___C1+___C5*ed-___C6*ee)+eb*(___C7+delp*(-___C8+delp*___C4))+delp*ea*(___C2-delp*___C3)-___C2*delp*ec)/(ave*sqrt(ave));
	if (p <= 0.0) ans=a*(b*ans+3.0*(rcx-rf(xt,yt,zt)));
	return ans;
}

#define TINY_ 1.69e-38
#define SQRTNY_ 1.3e-19
#define BIG_ 3.e37
#define TNBG_ (TINY*BIG)
#define COMP1_ (2.236/SQRTNY)
#define COMP2_ (TNBG*TNBG/25.0)
#define THIRD_ (1.0/3.0)
#define C1_ 0.3
#define C2_ (1.0/7.0)
#define C3_ 0.375
#define C4_ (9.0/22.0)

double rc(double x, double y)
//Computes Carlson?s degenerate elliptic integral, RC(x, y). x must be nonnegative and y must
//be nonzero. If y < 0, the Cauchy principal value is returned. TINY must be at least 5 times
//the machine underflow limit, BIG at most one fifth the machine maximum overflow limit.
{
	double alamb,ave,s,w,xt,yt;
//	if (x < 0.0 || y == 0.0 || (x+fabs(y)) < TINY || (x+fabs(y)) > BIG || (y<-COMP1 && x > 0.0 && x < COMP2)) nrerror("invalid arguments in rc");
	if (y > 0.0) {
		xt=x;
		yt=y;
		w=1.0;
	} else {
		xt=x-y;
		yt = -y;
		w=sqrt(x)/sqrt(xt);
	}
	do {
		alamb=2.0*sqrt(xt)*sqrt(yt)+yt;
		xt=0.25*(xt+alamb);
		yt=0.25*(yt+alamb);
		ave=THIRD_*(xt+yt+yt);
		s=(yt-ave)/ave;
	} while (fabs(s) > _ERRTOL);
	return w*(1.0+s*s*(C1_+s*(C2_+s*(C3_+s*C4_))))/sqrt(ave);
}

double ellf(double phi, double ak)
//Legendre elliptic integral of the 1st kind F(?, k), evaluated using Carlson?s function RF. The
//argument ranges are 0 ? ? ? ?/2, 0 ? k sin ? ? 1.
{
	double rf(double x, double y, double z);
	double s;
	s=sin(phi);
	return s*rf(SQR(cos(phi)),(1.0-s*ak)*(1.0+s*ak),1.0);
}

double elle(double phi, double ak)
//Legendre elliptic integral of the 2nd kind E(?, k), evaluated using Carlson?s functions RD and
//RF . The argument ranges are 0 ? ? ? ?/2, 0 ? k sin ? ? 1.
{
	double cc,q,s;
	s=sin(phi);
	cc=SQR(cos(phi));
	q=(1.0-s*ak)*(1.0+s*ak);
	return s*(rf(cc,q,1.0)-(SQR(s*ak))*rd(cc,q,1.0)/3.0);
}

double ellpi(double phi, double en, double ak)
//Legendre elliptic integral of the 3rd kind ?(?, n, k), evaluated using Carlson?s functions RJ and
//RF . (Note that the sign convention on n is opposite that of Abramowitz and Stegun.) The
//ranges of ? and k are 0 ? ? ? ?/2, 0 ? k sin ? ? 1.
{
	double cc,enss,q,s;
	s=sin(phi);
	enss=-en*s*s;
	cc=SQR(cos(phi));
	q=(1.0-s*ak)*(1.0+s*ak);
	return s*(rf(cc,q,1.0)-enss*rj(cc,q,1.0,1.0+enss)/3.0);
}