// C++ translation of FORTRAN PLS algorithm coded in 1985 by Ildiko Frank

#include <fstream>
#include <iostream>
#include <cmath>

#define FORARR(p,dim1size,dim1,dim2) ( p[(dim1size)*(dim2)+(dim1)] )

void matrix( int *pnvar1, int *pnvar2, int *picomp, double *beta, 
	double *w1, double *w2, double *b, double *ro, double *z ) {
// ************************************************************

 int nvar1 = *pnvar1;
 int nvar2 = *pnvar2;
 int icomp = *picomp + 1;

 int temp= nvar1* nvar2;
  /* Initialize z and beta  */
 for (int i=0;i<temp;i++) beta[i]=0. ;

 for (int ii=0;ii< nvar1; ii++)
  { for (int i=0;i<nvar1;i++) z[i] = 0.;
    z[ii]= 1.0 ;
     /* start the summation over icomp */
    for (int ic=0;ic<icomp;ic++)
     {double ss = 0.0 ;
        /* Dummy latent variable */
       for (int i=0;i<nvar1;i++)  ss += z[i] * FORARR(w1,nvar1, i,ic) ;
       ss *= ro[ic];
       for (int k=0;k<nvar2;k++)
        FORARR(beta,nvar1, ii,k) += ss * FORARR(w2,nvar2, k,ic);
        /* Update z matrix */
       ss=0;
       for (int j=0;j<nvar1;j++) ss   += z[j] * FORARR(w1,nvar1, j,ic)  ;
       for (int j=0;j<nvar1;j++) z[j] -= ss * FORARR(b,nvar1, j,ic) * ro[ic] ;
     } /* ic loop */
   }   /* ii loop */
}

double r1unif( int *pcseed, int *ptseed, int ibyte[4], int *pfcn) {
// ******************************************************
 double r1unif;
 int i,k,icarry,i1,i2,j1,j2,it1,it2;
 static int first = true;
 static int cseed[6] ={0},
            tseed[32]={0},
            iscr[5],
            t[29] = {1,2,3,3,2,1,4,5,6,7,5,4,7,6,1,6,7,4,5,2,3,7,6,5,4,3,2,1,0},
            jcseed = 12345,
            jtseed = 1073   ;
 r1unif = 0.0 ;

 if ( first || *pfcn>0 )
  { if (*pfcn <= 0) { jcseed = abs(*pcseed); jtseed = abs(*ptseed); }
    first = false;
    cseed[0] = jcseed;
    for (i=0;i<5;i++)
     { cseed[i+1] = cseed[i]  /64;
       cseed[i]  -= cseed[i+1]*64;  }
    cseed[5] = cseed[5] % 4;
    if (jcseed && (cseed[0] % 2)==0 ) cseed[0] ++;
    tseed[0] = jtseed;
    for (i=0;i<11;i++)
     { tseed[i+1] = tseed[i]  /2;
       tseed[i]  -= tseed[i+1]*2; }
    for (i=11;i<32;i++) tseed[0] = 0;
    if (jtseed) tseed[0] = 1;
    if ( ! *pfcn ) return r1unif;
  }

 for (i= 0;i<17;i++) tseed[i] = abs ( tseed[i] - tseed[i+15] );
 for (i=17;i<32;i++) tseed[i] = abs ( tseed[i] - tseed[i-17] );
 cseed[5] = 13* cseed[5] + 55* cseed[4] + 16* cseed[3];
 cseed[4] = 13* cseed[4] + 55* cseed[3] + 16* cseed[2];
 cseed[3] = 13* cseed[3] + 55* cseed[2] + 16* cseed[1];
 cseed[2] = 13* cseed[2] + 55* cseed[1] + 16* cseed[0];
 cseed[1] = 13* cseed[1] + 55* cseed[0] ;
 cseed[0] = 13* cseed[0] ;

 k= -6;
 icarry =0;
 for (i=0;i<5;i++)
  { k += 6;
    cseed[i] += icarry;
    icarry    = cseed[i]/64;
    cseed[i] -= 64*icarry;
    i2 = cseed[i] / 8;
    i1 = cseed[i] - 8*i2;
    j1 = 4*tseed[k+2]+tseed[k+1]+tseed[k+1]+tseed[k  ];
    j2 = 4*tseed[k+5]+tseed[k+4]+tseed[k+4]+tseed[k+3];
    it1 = 28;
    if (i1>j1) it1 = (i1*i1-i1)/2+j1;
    if (i1<j1) it1 = (j1*j1-j1)/2+i1;
    it2 = 28;
    if (i2>j2) it2 = (i2*i2-i2)/2+j2;
    if (i2<j2) it2 = (j2*j2-j2)/2+i2;
    iscr[i] = 8*t[it2]+t[it1];
    r1unif = (r1unif + iscr[i])  / 64.0 ;
  } /* i loop */
 cseed[5] = (cseed[5]+icarry) % 4;
 j1 = tseed[30]+tseed[31]+tseed[31];
 it1 = abs(cseed[5]-j1);
 if (it1 == 1 && (cseed[5]+j1)==3 ) it1=3;
 r1unif = (r1unif + it1) / 4.0 ;
 if ( *pfcn ==1 ) return r1unif;
 ibyte[3] = iscr[0] + (iscr[1] % 4) * 64 ;
 ibyte[2] = iscr[1]/4+(iscr[2] % 16)* 16 ;
 ibyte[1] = iscr[2]/16+iscr[3]*4;
 ibyte[0] = iscr[4]+it1*64;
 return r1unif;
}

void ranums( double *x, int *pn) {
// *******************************************************
 int i;
 static int ibyte[4],
            icseed   = 0,
            itseed   = 0,
            ifcn     = 1 ;
 for (i=0;i< *pn; i++) x[i] = r1unif( &icseed, &itseed, ibyte, &ifcn);
}

void resid( int *ppat, int *pvar1, int *pvar2, double *x, double *y, 
	double *weyt, double *xbar, double *xscal, double *off, double *beta, 
	double *res, double *ypred, double *ss, int *perr ) {
// *****************************************************

 *perr=0;
  int nvar1 = *pvar1;
  int nvar2 = *pvar2;

 for (int j=0;j< nvar2;j++) {
    ss[j] = 0. ;
    int jj = j+nvar1;
    off[j] = xbar[jj];
    for (int k=0;k<nvar1;k++) {
        FORARR(beta,nvar1,k,j) *= (xscal[jj]/xscal[k]) ;
       off[j] -= FORARR(beta,nvar1,k,j) *xbar[k];
    }
  } /* j loop */
 for (int i=0;i< *ppat;i++) {
     if (weyt[i]>=0.0) continue;
     for (int j=0;j<nvar2;j++) {
            double s  = off[j];
            for (int k=0;k<nvar1;k++)
                s += FORARR(x,nvar1,k,i) * FORARR(beta,nvar1,k,j);
            FORARR(ypred,nvar2,j,i) = s;
            FORARR(res,nvar2+1,j,i) = FORARR(y    ,nvar2,j,i)
                               - FORARR(ypred,nvar2,j,i)   ;
            s = FORARR(res,nvar2+1,j,i);
            ss[j] += s*s * fabs(weyt[i]) ;
     } /* j loop */
   }   /* i loop */
}

void sort(double *v, int *a, int ii, int *jj) {
// ********************************************************
int i, ij, j, k, l, m, t, tt;
int il[20], iu[20];
double vt, vtt;

/* following two lines are only to convince lint that we aren't using
   uninitialized variables:                                             */
  tt = 123456789;
 vtt = 123456789.0;

/*     PUTS INTO A THE PERMUTATION VECTOR WHICH SORTS V INTO */
/*     INCREASING ORDER.  ONLY ELEMENTS FROM II TO JJ ARE CONSIDERED. */
/*     ARRAYS IU(K) AND IL(K) PERMIT SORTING UP TO 2**(K+1)-1 ELEMENTS */

/*     THIS IS A MODIFICATION OF CACM ALGORITHM #347 BY R. C. SINGLETON, */
/*     WHICH IS A MODIFIED HOARE QUICKSORT. */
// remains in F77=>C syntax
m = 1;
i = (ii);
j = (*jj);
L_10:
	if( i >= j )
		goto L_80;
L_20:
	k = i;
ij = (j+i)/2;
t = a[ij-1];
vt = v[ij-1];
if( v[i-1] <= vt )
	goto L_30;
a[ij-1] = a[i-1];
a[i-1] = t;
t = a[ij-1];
v[ij-1] = v[i-1];
v[i-1] = vt;
vt = v[ij-1];
L_30:
	l = j;
if( v[j-1] >= vt )
	goto L_50;
a[ij-1] = a[j-1];
a[j-1] = t;
t = a[ij-1];
v[ij-1] = v[j-1];
v[j-1] = vt;
vt = v[ij-1];
if( v[i-1] <= vt )
	goto L_50;
a[ij-1] = a[i-1];
a[i-1] = t;
t = a[ij-1];
v[ij-1] = v[i-1];
v[i-1] = vt;
vt = v[ij-1];
goto L_50;
L_40:
	a[l-1] = a[k-1];
a[k-1] = tt;
v[l-1] = v[k-1];
v[k-1] = vtt;
L_50:
	l = l - 1;
if( v[l-1] > vt )
	goto L_50;
tt = a[l-1];
vtt = v[l-1];
L_60:
	k = k + 1;
if( v[k-1] < vt )
	goto L_60;
if( k <= l )
	goto L_40;
if( l - i <= j - k )
	goto L_70;
il[m-1] = i;
iu[m-1] = l;
i = k;
m = m + 1;
goto L_90;
L_70:
	il[m-1] = k;
iu[m-1] = j;
j = l;
m = m + 1;
goto L_90;
L_80:
	m = m - 1;
if( m == 0 )
	return;
i = il[m-1];
j = iu[m-1];
L_90:
	if( j - i > 10 )
		goto L_20;
if( i == (ii) )
	goto L_10;
i = i - 1;
L_100:
	i = i + 1;
if( i == j )
	goto L_80;
t = a[i+1-1];
vt = v[i+1-1];
if( v[i-1] <= vt )
	goto L_100;
k = i;
L_110:
	a[k+1-1] = a[k-1];
v[k+1-1] = v[k-1];
k = k - 1;
if( vt < v[k-1] )
	goto L_110;
a[k+1-1] = t;
v[k+1-1] = vt;
goto L_100;
} /* end of func */

void qq( int *ppat, double *w, double *q1, int *perr) {
// ***************************************************

 *perr = 0;
 int npat = *ppat;
 double xnpat = 0. ;
 for (int i=0;i<npat;i++) xnpat += w[i];
 if (xnpat <= 0.) { *perr=1; return; }
 double u = -0.5 * w[0]/xnpat;
 for (int i=0;i<npat;i++) {
    u += w[i]/xnpat;
    double t = (u>0.5) ? sqrt( -2.0 * log(1.0-u)) : sqrt( -2.0 * log(u)) ;
    q1[i] = t - (2.515517 + 0.802853*t + 0.010328*t*t )
               /(1. + 1.432788*t + 0.189269*t*t + 0.001308*t*t*t) ;
    if (u<=0.5) q1[i] =  -q1[i];
  }
}

void see(double *ss, double *press, int *pnss, int *pnpress, int *pnpat, int *pncomp) {
// *******************************************************

    double RAT = *pnpat - *pncomp - 1;
    if (RAT <= 0.) RAT = 1.0;
    RAT = *pnpat / RAT;
    for (int i = 0; ss && i < *pnss; i++) ss[i] = sqrt(ss[i] * RAT);
    for (int i = 0; press && i < *pnpress; i++) press[i] = sqrt(press[i] * RAT);
 }

double powi( double fptv, int intv)
// *********************************************************
{
    double t= 1.0;
 double a = abs( intv);
 for (int i=0; i<a; i++) t *= fptv ;
 if (intv<0) t = 1.0 / t;		/* zero divide! */
 return t;
}

/* f77c.h - definitions to support FORTRAN translations to C
			made by RTC PLUS from COBALT BLUE (Feb, 1987) */

#define mod(a,b)	((a) % (b))
#define lmod(a,b)	((a) % (b))

	/* externals to use in F77 std lib macros TO AVOID side effects */
	/* static int		f77_i_, f77_j_; */
	/* static long		f77_l_, f77_m_; */
	/* static double   	f77_d_, f77_e_; */


int nobootp(double *x, double *y, int *npat) {
// *****************************************************************

for(int i=0; i <= *npat; i+=1 ) x[i] = y[i];
return 1;
} 

int crossp(double *x, double *y, int *ix, int *npat, int *iout, int *iex, 
	int *ic) {
// **********************************************************
int ii;
int i, i1, i2;	/* changed from long for prototyping 7/13/95 */

/*  SET WEIGHTS FOR CROSS-VALIDATION */

int lmin2 = *iex < (*ic)  ? *iex : (*ic);
i1 = ((*ic)-1)*(*iout) + lmin2;
lmin2 = *iex < *ic ? *iex : *ic;
i2 = (*ic)*(*iout) + lmin2;

for(i=0; i < *npat ; i+=1 ){
	ii = ix[i];
	x[ii] = y[ii];
	if( i >= i1 && i <= i2 ) x[ii] =  - x[ii];
}

return 1;
} /* end of func */

int pcpls( int *nvar1, int *nvar2, int *npat, int *icomp, 
	int *nitmax, double *eps, double *x, double *y, double *weyt,
	double *w1, double *w2, double *b, double *ro, double *u,
	double *v, double *vt, int *ierr ) {
// *********************************************************
#define X(I_,J_)        (x+(I_)*( *nvar1)+(J_))
#define Y(I_,J_)        (y+(I_)*( *nvar2)+(J_))
#define W1(I_,J_)       (w1+(I_)*( *nvar1)+(J_))
#define W2(I_,J_)       (w2+(I_)*( *nvar2)+(J_))
#define B(I_,J_)        (b+(I_)*( *nvar1)+(J_))
#define U(I_,J_)        (u+(I_)*( *npat)+(J_))
#define V(I_,J_)        (v+(I_)*( *npat)+(J_))

/*     CALCULATES THE ICOMP-TH PLS COMPONENT */

(*ierr) = 0;

/*     Initialize V */

for(int i=0; i < *npat; i+=1 )
                 *V((*icomp),i) =  *Y(i,0);
int iter = 0;

    L_100:
	;
iter = iter + 1;

/*     Weights of first block */

for(int j=0; j < *nvar1; j+=1 ){
    double ss = 0.;
    for(int i=0; i < *npat; i+=1 ) if( weyt[i] >  0. )
	    ss = ss +  *X(i,j) * *V( *icomp, i) * weyt[i];
    *W1( *icomp, j ) = ss;
}
double ss = 0.;
for(int k=0; k < *nvar1; k+=1 )
    ss = ss + powi((double) *W1( *icomp, k ),2);
if( ss <= 0. ) {
    (*ierr) =  - 1; return 0; }

ss = sqrt(ss);
for(int k=0; k < *nvar1; k+=1 ) *W1( *icomp, k ) /=  ss;

/*     Latent variable of first block */
for(int i=0; i < *npat; i+=1 ) if( weyt[i] > 0. ) {
    ss = 0.;
    for(int k=0; k < *nvar1; k+=1 )
    	ss += *X(i,k) * *W1(*icomp,k);
    *U( *icomp, i) = ss;
}

/*     Norm of U */
double sumu = 0.;
for(int i=0; i < *npat; i+=1 ) if( weyt[i] > 0. )
	sumu += powi((double) *U(*icomp,i),2);

/*     Weights of second block */
for(int j=0; j < *nvar2; j+=1 ){
    ss = 0.;
    for(int i=0; i < *npat; i+=1 ) if( weyt[i] > 0. )
		ss += *Y(i,j) * *U(*icomp,i) *weyt[i];
    *W2( *icomp, j ) = ss;
}

ss = 0.;
for(int k=0; k < *nvar2; k+=1 )
    ss += powi((double) *W2(*icomp,k),2);

if( ss <= 0. ){ (*ierr) =  - 2; return 0; }
ss = sqrt(ss);
for(int k=0; k < *nvar2; k+=1 ) *W2(*icomp, k) /= ss;

/*     Latent variable of second block */
for(int i=0; i < *npat; i+=1 ) if( weyt[i] > 0. ) {
	ss = 0.;
	for(int k=0; k < *nvar2; k+=1 )
        ss += *Y(i, k) * *W2(*icomp, k);
    vt[i] = ss;
}
/*     convergence? */
ss = 0.;
for(int i=0; i < *npat; i+=1 ) if( weyt[i] > 0. )
    ss += powi((double)(vt[i] - *V(*icomp,i)),2);
for(int i=0; i < *npat; i+=1 ) *V(*icomp,i) = vt[i];

if( ss > (*eps) && iter < (*nitmax) ) goto L_100;

/*     Inner relationship */
if( sumu <= 0. ){ (*ierr) =  - 3; return 0; }
ss = 0.;
for(int i=0; i < *npat; i+=1 ) if( weyt[i] > 0. )
    ss += *U( *icomp, i ) * *V( *icomp, i ) * weyt[i];

ro[(*icomp)] = ss/sumu;

/*     Norm of the model vector */
sumu = 0.;
for(int i=0; i < *npat; i+=1 ) if( weyt[i] > 0. )
    sumu += powi((double)( *U(*icomp,i) * ro[(*icomp)]),2) * weyt[i];

/*     Loadings */
if( sumu <= 0. ){ (*ierr) =  - 4; return 0; }
for(int j=0; j < *nvar1; j+=1 ){
    ss = 0.;
    for(int i=0; i < *npat; i+=1 ) if( weyt[i] > 0. )
	    ss += *X( i, j ) * *U(*icomp, i) * ro[(*icomp)] * weyt[i]/sumu;
    *B(*icomp,j) = ss;
}

/*     Residuals */
for(int i=0; i < *npat; i+=1 ) if( weyt[i] > 0. ){
    ss = ro[(*icomp)] * *U(*icomp,i);
    for(int j=0; j < *nvar1; j+=1 )
	    *X(i,j) -= *B(*icomp,j) * ss;
    for(int j=0; j < *nvar2; j+=1 )
	    *Y(i,j) -= *W2(*icomp,j) * ss;
}

return 1;
#undef V
#undef U
#undef B
#undef W2
#undef W1
#undef Y
#undef X

} /* end of func */

int randomp(double *x, double *y, int *ix, int *npat) {
// ********************************************************
 int i, ii;
/*  RANDOM NUMBER GENERATOR FOR DRAWING BOOTSTRAP SAMPLE */

ranums(x,npat);
for(i=0; i < *npat; i+=1 ) ix[i] = (*npat)*((int)x[i]) + 1;

for(i=0; i < *npat; i+=1 ) x[i] = 0.;

for(i=0; i < *npat; i+=1 ){
	ii = ix[i];
	x[ii] = x[ii] + y[ii];
}
return 1;
} /* end of func */

int bootpls( double *ss, double *r2, double *press, double *cr2, double *off, 
	double *beta, int *nopt, int *ncomp, int *iboot, int *nvar2,
	int *nvar1 )
// ********************************************************
{

/*     MEAN AND AVERAGE OF THE BOOTSTRAPED QUANTITIES */

/*  ZERO THE STORAGES */
for(int j=0; j < *nvar2; j+=1 ) {
    if (ss) { 
	FORARR( ss, *nvar2, *iboot, j ) = 0.0;
	FORARR( r2, *nvar2, *iboot, j ) = 0.0;
	FORARR( ss, *nvar2, (*iboot)+1, j ) = 0.0;
	FORARR( r2, *nvar2, (*iboot)+1, j ) = 0.0;
    }
    for(int k=0; press && k < *ncomp; k+=1 ){
	press[((*iboot))*(*ncomp)*(*nvar2+1) + k*(*nvar2+1) + j] = 0.;
	press[((*iboot)+1)*(*ncomp)*(*nvar2+1) + k*(*nvar2+1) + j] = 0.;
	cr2[((*iboot))*(*ncomp)*(*nvar2) + k*(*nvar2+1) + j] = 0.;
	cr2[((*iboot)+1)*(*ncomp)*(*nvar2) + k*(*nvar2+1) + j] = 0.;
    }
}

nopt[(*iboot)] = 0;
nopt[(*iboot)+1] = 0.;
for(int j=0; j < *nvar2; j+=1 ){
	FORARR( off, *nvar2, *iboot, j ) = 0.0;
	FORARR( off, *nvar2, (*iboot)+1, j ) = 0.0;
	for(int i=0; i < *nvar1; i+=1 ) {
	    beta[((*iboot))*(*nvar2)*(*nvar1) + j*(*nvar1) + i] = 0.;
	    beta[((*iboot)+1)*(*nvar2)*(*nvar1) + j*(*nvar1) + i] = 0.;
	}
}

/*  SUM THE ELEMENTS OF ALL VECTORS INTO THE IBOOT+1-TH ELEMENTS */
for(int ib=0; ib < *iboot; ib+=1 ) for(int k=0; k < *nvar2; k+=1 ){
    if (ss) { 
	FORARR( ss, *nvar2, *iboot, k ) = FORARR( ss, *nvar2, *iboot, k ) 
		+ FORARR( ss, *nvar2, ib, k );
	FORARR( r2, *nvar2, *iboot, k ) = FORARR( r2, *nvar2, *iboot, k ) 
		+ FORARR( r2, *nvar2, ib, k );
    }
    for(int j=0; press && j < *ncomp; j+=1 ){
	cr2[((*iboot))*(*ncomp)*(*nvar2) + j*(*nvar2)+k] =  
	    cr2[((*iboot))*(*ncomp)*(*nvar2)+j*(*nvar2)+k] + 
	    cr2[(ib)*(*ncomp)*(*nvar2)*(*nvar2)+k];
	press[((*iboot))*(*ncomp)*(*nvar2+1)+j*(*nvar2+1)+k] = 
	    press[((*iboot))*(*ncomp)*(*nvar2+1)+j*(*nvar2+1)+k] +  
	    press[ib*(*ncomp)*(*nvar2+1)+j*(*nvar2+1)+k];
    }	

    nopt[(*iboot)] = nopt[(*iboot)] + nopt[ib];
    for(int j=0; j < *nvar2; j+=1 ){
	FORARR( off, *nvar2, *iboot, j ) =
	    FORARR( off, *nvar2, *iboot, j ) + FORARR( off, *nvar2, ib, j);
	for(int i=0; i < *nvar1; i+=1 )
	    beta[((*iboot))*(*nvar2)*(*nvar1) + j*(*nvar1) + i] =  
		beta[((*iboot))*(*nvar2)*(*nvar1) + j*(*nvar1) + i] +  
		beta[ib*(*nvar2)*(*nvar1) + j*(*nvar1) +i];
    }
}

/*  CALCULATE AVERAGE */
double boo = (double) *iboot;
for(int k=0; k < *nvar2; k+=1 ){
    if (ss) {
	FORARR(ss, *nvar2, *iboot, k) /= boo;
	FORARR(r2, *nvar2, *iboot, k) /= boo;
    }
    for(int j=0; press && j < *ncomp; j+=1 ){
	press[((*iboot))*(*ncomp)*(*nvar2+1)+j*(*nvar2+1)+k] /= boo;
	cr2[((*iboot))*(*ncomp)*(*nvar2)+j*(*nvar2)+k] /= boo;
    }
}

nopt[(*iboot)] /= *iboot;
for(int j=0; j < *nvar2; j+=1 ){
    FORARR( off, *nvar2, *iboot, j ) /= boo;
    for(int i=0; i < *nvar1; i+=1 )
	beta[((*iboot))*(*nvar2)*(*nvar1) + j*(*nvar1) +i] /= boo;
}

/*  SUM THE SQUARED DEVIATION FROM THE AVERAGES IN THE IBOOT+2-TH ELEMENT */
for(int ib=0; ib < *iboot; ib+=1 ) {
  for(int k=0; k < *nvar2; k+=1 ){
    if (ss) { 
	FORARR( ss, *nvar2, *iboot+1, k ) =
 	    FORARR( ss, *nvar2, *iboot+1, k ) +
powi((double)( FORARR( ss, *nvar2, ib, k ) - FORARR( ss, *nvar2, *iboot, k)), 2);
	FORARR( ss, *nvar2, *iboot+1, k ) =
 	    FORARR( ss, *nvar2, *iboot+1, k ) +
powi((double)( FORARR( ss, *nvar2, ib, k ) - FORARR( ss, *nvar2, *iboot, k)), 2);
    }
    for(int j=0; press && j < *ncomp; j+=1 ){
	cr2[((*iboot)+1)*(*ncomp)*(*nvar2)+j*(*nvar2)+k] =  
	    cr2[((*iboot)+1)*(*ncomp)*(*nvar2)+j*(*nvar2)+k] + 
	    powi((double)( cr2[ib*(*ncomp)*(*nvar2)+j*(*nvar2)+k] -
	       cr2[((*iboot))*(*ncomp)*(*nvar2)+j*(*nvar2)+k]),2);
	press[((*iboot)+1)*(*ncomp)*(*nvar2+1)+j*(*nvar2+1)+k] =  
	    press[((*iboot)+1)*(*ncomp)*(*nvar2+1)+j*(*nvar2+1)+k] + 
	    powi((double)( press[ib*(*ncomp)*(*nvar2+1)+j*(*nvar2+1)+k] - 
		press[((*iboot))*(*ncomp)*(*nvar2+1)+j*(*nvar2+1)+k]),2);
    }
  }
  nopt[(*iboot)+1] = nopt[(*iboot)+1] + (int) powi((double)(nopt[ib] -nopt[(*iboot)]),2);
  for(int j=0; j < *nvar2; j+=1 ){
	FORARR(off, *nvar2, *iboot+1, j) =
	    FORARR( off, *nvar2, *iboot+1, j ) +
	    powi((double)( FORARR(off,*nvar2,ib,j) - FORARR(off,*nvar2,*iboot,j)),2);
	for(int i=0; i < *nvar1; i+=1 )
	    beta[((*iboot)+1)*(*nvar2)*(*nvar1) + j*(*nvar1) + i] = 
    		beta[((*iboot)+1)*(*nvar1)*(*nvar2) + j*(*nvar1) + i] + 
    		powi((double)( beta[ib*(*nvar2)*(*nvar1) + j*(*nvar1) + i] - 
		    beta[((*iboot))*(*nvar2)*(*nvar1) + j*(*nvar1) + i]),2);
  }
}

/*  CALCULATE VARIANCE */
for(int k=0; k < *nvar2; k+=1 ){
   if (ss) {
	FORARR( ss, *nvar2, *iboot+1, k )
	     = sqrt( FORARR( ss, *nvar2, *iboot+1, k )/boo);
	FORARR( r2, *nvar2, *iboot+1, k ) =
	    sqrt( FORARR( r2, *nvar2, *iboot+1, k )/boo );
   }
   for(int j=0; press && j < *ncomp; j+=1 ){
	press[((*iboot)+1)*(*ncomp)*(*nvar2+1)+j*(*nvar2+1)+k] = 
	    sqrt( press[((*iboot)+1)*(*ncomp)*(*nvar2+1)+j*(*nvar2+1)+k]/boo);
	cr2[((*iboot)+1)*(*ncomp)*(*nvar2)+j*(*nvar2)+k] = 
	    sqrt( cr2[((*iboot)+1)*(*ncomp)*(*nvar2)+j*(*nvar2)+k]/boo);
   }
}

nopt[(*iboot)+1] = (int)(sqrt(float(nopt[(*iboot)+1])/boo));
for(int j=0; j < *nvar2; j+=1 ){
   FORARR( off, *nvar2, *iboot+1, j ) = 
	sqrt( FORARR( off, *nvar2, *iboot+1, j )/boo);
   for(int i=0; i <= *nvar1; i+=1 )
	beta[((*iboot)+1)*(*nvar2)*(*nvar1) + j*(*nvar1)+i] = 
	    sqrt( beta[((*iboot)+1)*(*nvar2)*(*nvar1) + j*(*nvar1) + i]/boo);
}

return 1;
} /* end of func */


int plsjer( int *pnpat, int *pnvar1, int *pnvar2, int *pncomp, int *nitmax, 
	double *eps, int *iboot, int *icros, int *icent, double *x, double *y, 
	double *weyt, double *weytb, double *xbar, double *xscal, 
	double *off, double *beta, double *varnc, double *w1, double *w2, 
	double *b, double *u, double *v, double *ro, double *ss, double *r2, 
	double *press, double *cr2, int *nopt, double *ypred, double *res, 
	double *sss, double *ssy, int *ierr, double *scrtch, int *ix )
{
#define X(I_,J_)	(x+(I_)*(nvar1)+(J_))
#define Y(I_,J_)	(y+(I_)*(nvar2)+(J_))
#define OFF(I_,J_)	(off+(I_)*(nvar2)+(J_))
#define BETA(I_,J_,K_)	(beta+(I_)*(nvar2)*(nvar1)+(J_)*(nvar1)+(K_))
#define VARNC(I_,J_)	(varnc+(I_)*(nvar2)+(J_))
#define W1(I_,J_)	(w1+(I_)*(nvar1)+(J_))
#define W2(I_,J_)	(w2+(I_)*(nvar2)+(J_))
#define B(I_,J_)	(b+(I_)*(nvar1)+(J_))
#define U(I_,J_)	(u+(I_)*(npat)+(J_))
#define V(I_,J_)	(v+(I_)*(npat)+(J_))
#define SS(I_,J_)	(ss+(I_)*(nvar2)+(J_))
#define R2(I_,J_)	(r2+(I_)*(nvar2)+(J_))
#define PRESS(I_,J_,K_)	(press+(I_)*(ncomp)*(nvar2+1)+(J_)*(nvar2+1)+(K_))
#define CR2(I_,J_,K_)	(cr2+(I_)*(ncomp)*(nvar2)+(J_)*(nvar2)+(K_))
#define YPRED(I_,J_)	(ypred+(I_)*(nvar2)+(J_))
#define RES(I_,J_)	(res+(I_)*(nvar2+1)+(J_))
 /*
int npat,nvar1,nvar2,ncomp; 	// changed from long for prototyping 7/13/95
int i, ib, iboot1, ic, icomp, icros1, iex, iout, j, jj, kl, l,
	nn, novariance ;	// changed from long for prototyping 7/13/95
double atl, pmin,      s, sum, varmax, xnpat, xnpatt;
*/
/*     MAIN PLS DRIVER */

/* NPAT,NVAR1,NVAR2,NCOMP, */

int npat = *pnpat;
    int nvar1= *pnvar1;
    int nvar2 = *pnvar2;
    int ncomp = *pncomp;

(*ierr) = 0;
int icros1 = 0;
int iboot1 = 0;

int nn = npat*(nvar1 +nvar2 +1) + 1;

/*  calculate quantiles for the q-q plot */
qq(&npat,weyt,scrtch,ierr);
if( (*ierr) != 0 ) return 0;

//for(int i=0; i < npat; i+=1 ) *RES(i,nvar2) = scrtch[i];

/*  set flags for no cross-validation and no bootstrap */
if( (*icros) == 0 ){ (*icros) = 1; icros1 = 1; }
if( (*iboot) == 0 ){ (*iboot) = 1; iboot1 = 1; }

/*  bootstrap monster loop */
for(int ib=0; ib < *iboot; ib+=1 ) {
    for(int icomp=0; icomp < ncomp; icomp+=1 )
	for(int j=0; press && j < nvar2+1; j+=1 )
	    press[ib*ncomp*(nvar2+1)+(icomp)*(nvar2+1)+j] = 0.;

	/*  no bootstrap, copy the weights */
	if( iboot1 == 1 ) for (int i = 0; i < npat; i++) weytb[i] = weyt[i];
	/*  draw a bootstrap sample , new weights in weytb */
	    else randomp(weytb,weyt,ix,&npat);

	/*  set up an auxiliary vector for cross-validation pointing to non zero weight */
    int kl = 0;
	if( icros1 != 1 ){
	    for(int i=0; i < npat; i+=1 ) if( weytb[i] > 0. ) {
			kl = kl + 1;
			ix[kl-1] = i;
	    }
	    ranums(scrtch,&kl);
	    sort(scrtch,ix,1,&kl);
	    for(int i=0; i < npat; i+=1 ) scrtch[i] = 0.;

		/*  calculate how many samples to delete */
	}
    int iout = kl/(*icros);
    int iex = lmod(kl,(*icros));

	/*  calculate sum of squares */
    double atl = 0.0;
    double xnpatt = 0.;
	for(int j=0; j < nvar2; j+=1 ){
	    ssy[j] = 0.;
	    for(int i=0; i < npat; i+=1 ){
		    xnpatt = xnpatt + weytb[i];
		    atl = atl +  *Y( i, j ) *weytb[i];
        }
	    atl = atl/xnpatt;
	    for(int i=0; i < npat; i+=1 )
		    ssy[j] += powi((double)( *Y(i,j) - atl ),2) *weytb[i];
	}

	/*  cross-validation loop */
	for(int ic=0; ic < *icros; ic+=1 ){
			/*  no cross-validation, copy the weights */
	    if( icros1 == 1 ) for (int i = 0; i < npat; i++) scrtch[i] = weyt[i];
		/*  set weights for cross-validation in scrtch */
		else crossp(scrtch,weytb,ix,&kl,&iout,&iex,&ic);

		/*  SUM OF WEIGHTS */
	    double xnpat = 0.;
	    for(int i=0; i < npat; i+=1 ) if( scrtch[i] > 0. ) xnpat += scrtch[i];
	    if( xnpat <= 0. ){ (*ierr) = 2; return 0; }

		/*  MEAN VALUES AND SCALE */
		/*  1 - X;  2 - X-XBAR;   3 - (X-XBAR)/XSCAL;    4 - X/XSCAL; */
		/*  CALCULATE XBAR THE MEAN */
	    if( (*icent) > 1 ){
		    for(int j=0; j < nvar1; j+=1 ){
		        double s = 0.;
		        for(int i=0; i < npat; i+=1 ) if( scrtch[i] > 0. )
			        s += *X(i,j) * scrtch[i];
		        xbar[j] = s/xnpat;
		    }
		    for(int j=0; j < nvar2; ++j) {
		        double s = 0.;
		        for(int i=0; i < npat; i+=1 ) if( scrtch[i] > 0. )
			        s += *Y(i,j) * scrtch[i];
		        int jj = j + nvar1;
		        xbar[jj] = s/xnpat;
		    }
	    }
		/*   CALCULATE XSCAL THE VARIANCE */
	    if( (*icent) >= 3 ) {
		    int novariance = 0;
		    for(int j=0; j < nvar1; j+=1 ){
		        double s = 0.;
		        for(int i=0; i < npat; i+=1 ) if( scrtch[i] > 0. )
			        s += powi((double)(*X(i,j)- xbar[j]),2)*scrtch[i];
                xscal[j] = sqrt(s/(xnpat -1));
                if( s <= 0. ) novariance = novariance + 1;
            }
            if( novariance != 0 ) {
                double varmax = 1.e20;
                for(int j=0; j < nvar1; j+=1 )
                    if( xscal[j] > 0.0 && xscal[j] < varmax )
				            varmax = xscal[j];
                if( varmax == 1.e20 ){ (*ierr) = 2; return 0; }
                varmax = varmax/1000.;
                for(int j=0; j < nvar1; j+=1 ) if( xscal[j] <= 0.0 )
				        xscal[j] = varmax;
            }
            for(int j=0; j < nvar2; j+=1 ){
                int jj = j + nvar1;
                double s = 0.;
                for(int i=0; i < npat; i+=1 ) if( scrtch[i] > 0. )
			            s += powi((double)( *Y(i,j) - xbar[jj]), 2) * scrtch[i];
                xscal[jj] = sqrt(s/(xnpat -1));
                if( xscal[j] == 0 ){ (*ierr) = 2; return 0; }
            }
        }
		/*  ICENT = 1 NO CENTERING OR SCALING */
	    if( (*icent) == 1 ){
		    for(int i=0; i < npat; i+=1 ){
		        for(int j=0; j < nvar1; j+=1 ) {
                    int l = i*nvar1 + j +npat;
                    scrtch[l] = *X(i,j);
                }
		        for(int j=0; j < nvar2; j+=1 ){
			        int l = npat + npat*nvar1 + (i )*nvar2 + j;
			        scrtch[l] =  y[(i)*nvar2+j];
		        }
		    }

		    for(int j=0; j < nvar1+nvar2; j+=1 ){
		        xbar[j] = 0.;
		        xscal[j] = 1.;
		    }
	    }
		/*  ICENT = 2 */
	    if((*icent) == 2 ) for(int i=0; i < npat; i+=1 ){
		    for(int j=0; j < nvar1; j+=1 ){
			    int l = (i )*nvar1 + j + npat;
			    scrtch[l] =  *X(i,j) - xbar[j];
		    }
		    for(int xj=0; xj < nvar2; xj+=1 ){
			    int jj = xj + nvar1;
			    int xl = npat + npat*nvar1 + (i )*nvar2 + xj;
			    scrtch[xl] =  *Y(i,xj) - xbar[jj];
		    }
		    for(int j=0; j < nvar1+nvar2; j+=1 ) xscal[j] = 1.;
	    }
		/*  ICENT = 3 */
	    if( (*icent) == 3 ) for(int i=0; i < npat; i+=1 ){
		    for(int j=0; j < nvar1; j+=1 ){
			    int l = (i )*nvar1 + j + npat;
			    scrtch[l] = ( *X(i,j) - xbar[j])/xscal[j];
		    }
		    for(int j=0; j < nvar2; j+=1 ){
			    int jj = j + nvar1;
			    int l = npat + npat*nvar1 + (i )*nvar2 + j;
			    scrtch[l] = ( *Y(i,j) - xbar[jj])/xscal[jj];
		    }
	    }
		/*  ICENT = 4 */
	    if( (*icent) == 4 ){
		    for(int i=0; i < npat; i+=1 ){
		        for(int j=0; j < nvar1; j+=1 ){
			        int l = (i )*nvar1 + j + npat;
			        scrtch[l] =  *X(i,j)/ xscal[j];
		        }
		        for(int j=0; j < nvar2; j+=1 ){
		    	    int jj = j + nvar1;
		    	    int l = npat + npat*nvar1 + (i )*nvar2 + j;
		    	    scrtch[l] =  *Y(i,j) / xscal[jj];
		        }
		    }
		    for(int j=0; j < nvar1+nvar2; j+=1 ) xbar[j] = 0.;
	    }

		/*  CALCULATE MODEL OF NCOMP COMPONENTS */
	    for(int icomp=0; icomp < ncomp; icomp+=1 ){
            cout << "Extracting component " << icomp+1 << endl;
		    if (!pcpls(&nvar1,&nvar2,&npat,&icomp,nitmax,eps,&scrtch[npat],
			  &scrtch[npat*(nvar1+1)],scrtch,w1,w2,b,
			  ro,u,v,&scrtch[nn-1],ierr)) return 0;
		    if( (*ierr) != 0 ){
		        int xncomp = icomp - 1;
		        if( xncomp < 0 ) return 0;
		        (*ierr) = 0;
		        std::cout << "No more than " << xncomp << "components can be calculated\n";
		        goto L_222;
		    }
		    for(int j=0; j < nvar2; j+=1 ){
		        double sum = 0.;
		        int jj = j + nvar1;
		        for(int i=0; i < npat; i+=1 ) {
			        int l = npat + npat*nvar1 + (i )*nvar2 + j;
			        sum += powi((double)scrtch[l],2)*scrtch[i];
		        }
		        if( (*icent) >= 3 ) {
			        sum *= powi((double)xscal[jj],2);
			        *VARNC(icomp,j) = 1. - sum/ssy[j];
		        }
		    }
	    }
		L_222:
			;
		/*  load the data again for residual calculation */
	    for(int i=0; i < npat; i+=1 ){
		    for(int j=0; j < nvar1; j+=1 ){
		        int l = (i )*nvar1 + j + npat;
		        scrtch[l] = *X(i,j );
		    }
		    for(int j=0; j < nvar2; j+=1 ){
		        int l = (i )*nvar2 + j + npat*(nvar1 +1);
		        scrtch[l] =  *Y( i, j );
		    }
	    }
	/*  samples with negative weights are predicted, so set all of them to - */
	    if( icros1 == 1 ) for(int i=0; i < npat; i+=1 ) scrtch[i] =  - scrtch[i];

		/*  REGRESSION MATRIX AND PREDICTION */
	    for(int icomp=0; icomp < ncomp; icomp+=1 ){
		    if( iboot1 != 1 ){
			    matrix(&nvar1,&nvar2,&icomp,&beta[(ib)*nvar1*nvar2 ],
				  w1,w2,b,ro,&scrtch[nn-1]);
			    resid(&npat,&nvar1,&nvar2,&scrtch[npat],&scrtch[npat*(nvar1+1)+1-1],
                      scrtch,xbar,xscal,&off[nvar2*ib],
                      &beta[(ib)*nvar1*nvar2],res,ypred,sss,ierr);
		    }
		    else{
			    matrix(&nvar1,&nvar2,&icomp,&beta[(ic)*nvar2*nvar1],
				  w1,w2,b,ro,&scrtch[nn-1]);
			    resid(&npat,&nvar1,&nvar2,&scrtch[npat],
                 &scrtch[npat*(nvar1+1)],scrtch,xbar,xscal,&off[nvar2*ic],
                 &beta[(ic)*nvar1*nvar2],res,ypred,sss,ierr);
		    }
		    if( (*ierr) != 0 ) return 0;
		    if( icros1 == 0 )
				/*  sum of cross-validated squared residuals */
		        for(int j=0; press && j < nvar2; j+=1 )
			        *PRESS(ib,icomp,j) +=  sss[j];
	    }
	    if( icros1 == 1 ){
			/* sum of squared residuals and r squared */
		    for(int j=0;  j < nvar2; j+=1 ){
		        *SS(ib,j) = sss[j];
		        *R2(ib,j) = 1. - *SS(ib,j) /ssy[j];
		    }
		    for(int j=0; j < nvar2; j+=1 ) *SS( ib,j ) /= xnpatt;
	    }
	}

	/*  find the optimal number of components */
	if( icros1 != 1 ){
	    double pmin = 1.e10;
	    for(int icomp=0; icomp < ncomp; icomp+=1 ){
		double s = 0.;
		for(int j=0; j < nvar2; j+=1 )
		    s += *PRESS(ib,icomp,j);
            *PRESS(ib,icomp,nvar2+1) = s/( (double) nvar2);
            if (*PRESS(ib,icomp,nvar2+1) < pmin ){
				pmin = *PRESS(ib,icomp,nvar2+1);
				nopt[ib] = icomp;
		    }
	    }
	}
	/*  CROSS-VALIDATED R SQUARED */
	if( icros1 != 1 ){
	    for(int icomp=0; icomp < ncomp; icomp+=1 )
		    for(int j=0; j < nvar2; j+=1 )
		        *CR2(ib,icomp,j) = 1. - *PRESS(ib,icomp,j) /( ssy[ j]);
	        for(int j=0; j < nvar2+1; j+=1 ){
		        double undo;
		        if (npat-ncomp > 0) undo = (double) npat-ncomp;
			        else {undo = 1.0; ncomp = npat - 2; }
		        for(int icomp=0; icomp < ncomp; icomp+=1 )
                    *PRESS(ib,icomp,j) = *PRESS(ib,icomp,j) /xnpatt * ( undo / (double)(npat-icomp));
	        }
    }
}
if( iboot1 == 1 ) (*iboot) = 0;
if( icros1 == 1 ) (*icros) = 0;
/*  calculate the mean and the variance of the bootstrapped quantities */
if( iboot1 != 1 )
	bootpls(ss,r2,press,cr2,off,beta,nopt,&ncomp,iboot,&nvar2,&nvar1);

int tmpa_55 = ((*iboot) +2)*nvar2;
    int tmpa_56 = (nvar2 +1)*ncomp*((*iboot) +2);
see(ss,press,&tmpa_55,&tmpa_56,&npat,&ncomp);	/*SEE (sqrt) *///

return 1;
#undef RES
#undef YPRED
#undef CR2
#undef PRESS
#undef R2
#undef SS
#undef V
#undef U
#undef B
#undef W2
#undef W1
#undef VARNC
#undef BETA
#undef OFF
#undef Y
#undef X
} /* end of func */

/*
int main() {
// for testing
   double xs[9] = {3.0,2.0,3.0, 4.0,5.0,5.0, 6.0,7.0,6.0};
   double ys[3] = {3.0,4.0,5.0};
   int nOKMol = 3;
   int nOKfVal = 3;

   double ywts[nOKMol];
   int ierr;
   int nys = 1;
   int mxncomp = 2;
   int maxiter = 100;
   double eps =1.0e-4;
   int iboot = 0;
   int icros = 0;
   int icent = 3;
   double *weytb = new double[nOKMol];
   double *xbar = new double[nOKMol+1];
   double *xscal = new double[nOKMol+1];
   double *intcpt = new double[3];
   double *coeff = new double[2*nOKfVal];
   double *varnc = new double[mxncomp*nOKfVal];
   double *wtx = new double[mxncomp*nOKfVal];
   double *wty = new double[mxncomp];
   double *loadings = new double[mxncomp*nOKfVal];
   double *latentx = new double[mxncomp*nOKMol];
   double *latenty = new double[mxncomp*nOKMol];
   double *inner = new double[mxncomp];
   double *ssqRsdl = new double[mxncomp*6];
   double *r2 = new double[2];
   double *sdep = new double[nOKfVal*mxncomp];
   double *q2 = new double[2*mxncomp];
   int *optncomp = new int[2];
   double *ypred = new double[nOKMol];
   double *resdl = new double[2*nOKMol];
   double *sss = new double[2];
   double *ssy = new double[ nOKMol ];
   int *iscratch = new int[nOKMol];
   double *scratch = new double[nOKMol * (nOKfVal + 2) + (nOKfVal+1)];
   for (int iy = 0; iy < 3; ++iy) intcpt[iy] = 0.0;
   for (int iy = 0; iy < 2; ++iy) 
	{optncomp[iy] = 0; r2[iy]=sss[iy]=ssqRsdl[iy]=0.0;}
   for (int iy = 0; iy < nOKMol; ++iy) {iscratch[iy] = 0; 
	ssy[iy] = weytb[iy] = ypred[iy] = 0.0; ywts[iy] = 1.0;}
   for (int iy = 0; iy < 2*nOKMol; ++iy) resdl[iy] = 0.0;
   for (int iy = 0; iy < 2*nOKfVal; ++iy) coeff[iy] = 0.0; 
   for (int iy = 0; iy < mxncomp*nOKfVal; ++iy) 
	varnc[iy] = wtx[iy] = loadings[iy] = sdep[iy] = 0.0;
   for (int iy = 0; iy < mxncomp; ++iy) wty[iy] = inner[iy] = 0.0;
   for (int iy = 0; iy < mxncomp*nOKMol; ++iy) latentx[iy] = latenty[iy] = 0.0;
std::cout << nOKMol << ' ' << nOKfVal << ' ' << nOKMol * (nOKfVal + 2) + (nOKfVal+1) << std::endl;

   plsjer( &nOKMol, &nOKfVal, &nys, &mxncomp, &maxiter, &eps, &iboot,
        &icros, &icent, xs, ys, ywts, weytb, xbar, xscal, intcpt, coeff,
        varnc, wtx, wty, loadings, latentx, latenty, inner, ssqRsdl, r2,
           NULL, NULL, optncomp, ypred, resdl, sss, ssy, &ierr, scratch, iscratch);
std::cout << *r2 << ' ' << *ssqRsdl << "\n";
}

1201447-30072015
00001ZY40aiTiyS7T31eq1dj1m5ERE
Ml2EhOhukYoBZ47EDq5YKURQrGHC0F
NYhaxpIWRlGLX2UnxDWG8gVraYn7r0

1201447-16092015
00000ESP6dNrjLoE9poieD1NM"90o6
6hr3ChGO065!kEyofn6dVr1q41eY"x
ghyxDn"UTVF!6Z6a"vQ8OO6A4lfjMe
*/
