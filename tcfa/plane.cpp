#include <iostream>
#include <cmath>

using namespace std;

void UTL_GEOM_CROSS(
// ***************************************************************
double  *vector1,   // ptr to x,y,z components of first vector  
double  *vector2,   // ptr to x,y,z components of second vector 
double  *product)   // ptr to x,y,z components of cross product 
{

  *product     = *(vector1+1) * *(vector2+2) - *(vector1+2) * *(vector2+1);
  *(product+1) = *(vector1+2) * *vector2     - *vector1     * *(vector2+2);
  *(product+2) = *vector1     * *(vector2+1) - *(vector1+1) * *vector2;

}

double UTL_GEOM_VDIST (
// *********************************************************
const double *point1,        /* ptr to x,y,z of point 1 */
const double *point2)        /* ptr to x,y,z of point 2 */
{
  double diff1,        /* difference x1 - x2 */
         diff2,        /* difference y1 - y2 */
         diff3;        /* difference z1 - z2 */

  diff1 = *point1 - *point2;                 /* x1 - x2 */
  diff2 = *(point1 + 1) - *(point2 + 1);     /* y1 - y2 */
  diff3 = *(point1 + 2) - *(point2 + 2);     /* z1 - z2 */

  return (sqrt(diff1 * diff1 +
          diff2 * diff2 +             /* sum of squares */
          diff3 * diff3));
}

double UTL_GEOM_f2dVDIST (
// *********************************************************
const float *point1,        /* ptr to x,y,z of point 1 */
const float *point2)        /* ptr to x,y,z of point 2 */
{
  double diff1,        /* difference x1 - x2 */
         diff2,        /* difference y1 - y2 */
         diff3;        /* difference z1 - z2 */

  diff1 = *point1 - *point2;                 /* x1 - x2 */
  diff2 = *(point1 + 1) - *(point2 + 1);     /* y1 - y2 */
  diff3 = *(point1 + 2) - *(point2 + 2);     /* z1 - z2 */

  return ((double) sqrt(diff1 * diff1 +
          diff2 * diff2 +             /* sum of squares */
          diff3 * diff3));
}

double UTL_GEOM_VLEN (
// ****************************************************
double *vector)     /* ptr to x, y, and z components */
{
  double    vl;

  vl = (*vector * *vector) + (*(vector + 1) * *(vector + 1)) +
       (*(vector + 2) * *(vector + 2));
  if (vl > 1e-12)
     vl = sqrt (vl);
  else
     vl = 0.0;
  return (vl);
}

double UTL_GEOM_VNORM (
// *********************************************************
double *vector,   /* ptr to x,y,z components of vector */
double *nvector)  /* ptr to x,y,z components of normalized vector */
{
  double mag;              /* magnitude                  */
  int    i;                /* vector index               */

  mag = UTL_GEOM_VLEN (vector);
  if (mag > 0.000001)  {
    for (i=0; i<3; i++)
       *(nvector + i) = *(vector + i) / mag;
    }
  else  {
    for (i=0; i<3; i++)
      *(nvector + i) = 0.0;
    mag = 0.0;
    }

  return (mag);
}

static int isLinear(double *cord, int *planeAtoms, int testAtom, int *r_islinear )
// ********************************************************************************
/* returns true through int pionter isLinear if true, return value indicates if
2nd priority atom (testAtom) needs a reflection.  return value of True indicates reflection is needed. 
Is the plane along one of the axis, if so, it's simple to determine what is local RIGHT.
*/
{
	int i;
	double *cptr, *cptr2;
	double za, xa, zdiff, ztotal;
	double ytotal;
	double xdiff, xtotal;
	double zcords[3];
	double xcords[3];
	double ycords[3];
	double xdir, ydir;
	int linear;
    double dists[3], t;
    int max;

	*r_islinear = 0;

    /* Do the geometric test for linearity first */
    cptr = cord + (planeAtoms[0] * 3);
    cptr2 = cord + (planeAtoms[1] * 3);
    dists[0] = UTL_GEOM_VDIST(cptr, cptr2);
    max = 0;
    cptr2 = cord + (planeAtoms[2] * 3);
    dists[1] = UTL_GEOM_VDIST(cptr, cptr2);
    if (dists[1] > dists[0])
        max = 1;
    cptr = cord + (planeAtoms[1] * 3);
    dists[2] = UTL_GEOM_VDIST(cptr, cptr2);
    if (dists[2] > dists[max])
        max = 2;

    t = 0.0;
    for (i = 0; i < 3; i++) {
        if (i == max)
            t += dists[i];
        else
            t -= dists[i];
    }

    if (! (t > -0.001 && t < 0.001))
        return 0;

/* Now, the original algorithm, once true linearity has been confirmed */
	for ( xtotal = ytotal = ztotal = 0.0, i = 0; i < 3; i++ ) {
		cptr = cord + ( (planeAtoms[i] * 3));
		xcords[i] = *cptr;
		xtotal += *cptr;

		cptr++;
		ycords[i] = *cptr;
		ytotal += *cptr;

		cptr++;
		zcords[i] = *cptr;
		ztotal += *cptr;
	}
	ztotal /= 3.0;
	ytotal /= 3.0;
	xtotal /= 3.0;

	cptr = cord + (testAtom * 3);
	
/* first test Z, since this is the most common due to how we are aligning the compound 
If all the Z coords are close to one another, plane travelling along X axis.  */
	za = *(cptr+2);
	xdir = xcords[2] - xcords[0] + 0.2;
	for ( linear = 1, i = 0; i < 3 && linear; i++ ) {
		zdiff = ztotal - zcords[i];
		if ( zdiff > 0.2 || zdiff < -0.2 )
			linear = 0; 
	}
	if ( linear ) {

		if ( ( xdir > 0.0 && zdiff < -0.1 ) || ( xdir < 0.0 && zdiff > 0.1 ) )
			return 1;
		return 0;
	}

/* Not along the X or Y axis, test the Z axis */
	xa = *(cptr+1);
	ydir = ycords[2] - ycords[0];
	for ( linear = 1, i = 0; i < 3 && linear; i++ ) {
		xdiff = xtotal - xcords[i];
		if ( xdiff > 0.2 || xdiff < -0.2)
			linear = 0;
	}
	if ( linear ) {
		*r_islinear = 1;
		if ( (ydir > 0.0 && xa < xtotal ) || ( ydir < 0.0 && xa > xtotal ) )
			return 1;
		return 0;
	}
	return 0;
}

bool VECT_ISNOT_LOCAL_RIGHT(double *cord, int *planeAtoms, int testAtom, int testConnectedTo) {
// **********************************************************************************************
/* 
	THIS LOCAL GEOMETRY FUNCTION was developed to support topomer
   It determines what local direction the 4 vectors are going, and what dominates the orientation of the plane
	v1:  p2-p1  plane atom2 - atom1
        v2:  p3-p2  plane atom3 - atom2
	v3:  *-p2   testAtom - ( atom test atom is bonded to either atom2 or atom3. 
	v4:  p3-p1  General direction of the plane. 
	It also determines what controls the plane, the plane direction, and how to 
	determine if v3 is pointing to the local right. 
	returns: 
		0:  Don't need to reflect
		1:  Need to do a reflection.
Usage: 
		doReflection = VECT_ISNOT_LOCAL_RIGHT(double *cords, array of a1,a2,a3, atom id to test, test connected to 2 or 3

Usage 1: 
	after setting torsion for a0, a1, a2, and a3,  you want to test if the 2nd priority atom off of a2 is to the local right. 
	doReflection = VECT_ISNOT_LOCAL_RIGHT(double *cords, [a1,a2,a3], secondChoice[a2], 2 );
		where secondChoice is an array for each atoms 2nd choice. 
Usage 2: 
	After setting torsion for a0,a1,a2, and a3, you want to test if a3 ended up on the right, since a2-a3 is a bond in a ring. 
	doReflection = VECT_ISNOT_LOCAL_RIGHT(double *cords, [a0,a1,a2], a3, 3 );
Author: 
	Rob Jilek   4/2003
*/
	double absx, absy,absz; 
	int pmode;
	double v1[3],v2[3], v3[3];
	double nv1[3], nv2[3];
	double *c1, *c2, *c3, *c4;
	double ucv[3];
	double invert;
	int i, rc, islinear;

	rc = isLinear(cord, planeAtoms, testAtom, &islinear );
	if ( rc )
		return true;
	c1 = cord + (planeAtoms[0] * 3);
	c2 = cord + (planeAtoms[1] * 3);
	c3 = cord + (planeAtoms[2] * 3);
	c4 = cord + (testAtom * 3);

	for ( i = 0; i < 3; i++ ) {
		v1[i] = *(c2+i) - *(c1+i);
		v2[i] = *(c3+i) - *(c2+i);
		if ( testConnectedTo == 2 )
			v3[i] = *(c4+i) - *(c2+i);
		else
			v3[i] = *(c4+i) - *(c3+i);
		if (v1[i] < 0.1 && v1[i] > -0.1 ) v1[i] = 0.0;
		if (v2[i] < 0.1 && v2[i] > -0.1 ) v2[i] = 0.0;
		if (v3[i] < 0.1 && v3[i] > -0.1 ) v3[i] = 0.0;
	}

	UTL_GEOM_VNORM(v1, nv1);
	UTL_GEOM_VNORM(v2, nv2);
	UTL_GEOM_CROSS(nv1,nv2,ucv);

	absx = ucv[0];
	absy = ucv[1];
	absz = ucv[2];
	if ( absx < 0.0 ) absx *= -1.0;
	if ( absy < 0.0 ) absy *= -1.0;
	if ( absz < 0.0 ) absz *= -1.0;
	absz += 0.05; /* prefer Z */

	if ( absz >= absy && absz >= absx )
		pmode = 3;		/* we want Z to win for ties */
	else if ( absx > absy && absx > absz ) {
		if ( (absx - absz) < 0.05 ) pmode = 3;
		else
			pmode = 1;
	}
	else  /* absy must be the biggest */
	{
		if ( ( absy - absz ) < 0.1 )
			pmode = 3;
		else
			pmode = 2;
	}


	if ( pmode == 1 ) {
		invert = v1[2] * v3[0];
		if ( invert < 0.0 )
			return false;
		return 1;  /* Do reflection */
	}
	if ( pmode == 2 ) {
		invert = v1[0] * v3[1];
		if ( invert < 0.0 )
			return false;
		return true;  /* Do reflection. */
	}
	/* pmode must be 3 */
	invert = v1[0] * v3[2];
	if ( invert > 0.0 )
		return false;
	return true;  /* Do reflection */
}

int UTL_GEOM_SYMM_EIGENSYS (
double *a,
int    n,
double *d,
double *v,
int    *nrot)
{
   double   c, g, h, p, s, sm, t, tau, theta, thresh;
   int      i, ip, iq, j, k;
   bool converged = false;

   for (ip = 0;  ip < n;  ip++)  {
      for (iq = 0;  iq < n;  iq++)
         *(v + ip * n + iq) = (ip == iq) ? 1.0 : 0.0;
      }

   for (ip = 0;  ip < n;  ip++)  {
      *(d + ip) = *(a + ip * n + ip);
      }

   *nrot = 0;
   for (i = 0;  i < 50 && !converged;  i++)  {
      sm = 0.0;
      for (ip = 0;  ip < n - 1;  ip++)  {
         for (iq = ip + 1;  iq < n;  iq++)
            sm += abs(*(a + ip * n + iq));
      }
      if (sm == 0.0)
         converged = true;

      thresh = (i < 3) ? (0.2 * sm / ((double) (n * n))) : 0.0;

      for (ip = 0;  ip < n - 1;  ip++)  {
         for (iq = ip + 1;  iq < n;  iq++)  {
            g = 100.0 * abs(*(a + ip * n + iq));
            if ((i > 4) && (abs(*(d + ip)) + g == abs(*(d + ip)))
                        && (abs(*(d + iq)) + g == abs(*(d + iq))))
               *(a + ip * n + iq) = 0.0;
            else if (abs(*(a + ip * n + iq)) > thresh)  {
               h = *(d + iq) - *(d + ip);
               if (abs(h) + g == abs(h))
                  t = *(a + ip * n + iq) / h;  /* t = 1 / (2 * theta) */
               else  {
                  theta = 0.5 * h / *(a + ip * n + iq);  /* Eq 11.1.10 */
                  t = 1.0 / (abs(theta) + sqrt(1.0 + theta * theta));
                  if (theta < 0.0)  t = -t;
                  }
               c = 1.0 / sqrt(1.0 + t * t);
               s = t * c;
               tau = s / (1.0 + c);
               h = t * (*(a + ip * n + iq));
               *(d + ip) -= h;
               *(d + iq) += h;
               *(a + ip * n + iq) = 0.0;

               for (j = 0;  j < ip;  j++)  {
                  g = *(a + j * n + ip);
                  h = *(a + j * n + iq);
                  *(a + j * n + ip) = g - s * (h + g * tau);
                  *(a + j * n + iq) = h + s * (g - h * tau);
                  }

               for (j = ip + 1;  j < iq;  j++)  {
                  g = *(a + ip * n + j);
                  h = *(a + j * n + iq);
                  *(a + ip * n + j) = g - s * (h + g * tau);
                  *(a + j * n + iq) = h + s * (g - h * tau);
                  }

               for (j = iq + 1;  j < n;  j++)  {
                  g = *(a + ip * n + j);
                  h = *(a + iq * n + j);
                  *(a + ip * n + j) = g - s * (h + g * tau);
                  *(a + iq * n + j) = h + s * (g - h * tau);
                  }

               for (j = 0;  j < n;  j++)  {
                  g = *(v + j * n + ip);
                  h = *(v + j * n + iq);
                  *(v + j * n + ip) = g - s * (h + g * tau);
                  *(v + j * n + iq) = h + s * (g - h * tau);
                  }
               *nrot += 1;
               }
            }
         }
      }
   if (!converged) return false;

   for (ip = 0;  ip < n - 1;  ip++)  {
         for (iq = ip + 1;  iq < n;  iq++)
            *(a + ip * n + iq) = *(a + iq * n + ip);
         }

   for (i = 0;  i < n - 1;  i++)  {
         k = i;
         p = *(d + i);
         for (j = i + 1;  j < n;  j++)  {
            if (*(d + j) <= p)  {
               k = j;
               p = *(d + j);
               }
            }
         if (k != i)  {
            *(d + k) = *(d + i);
            *(d + i) = p;
            for (j = 0;  j < n;  j++)  {
               p = *(v + j * n + i);
               *(v + j * n + i) = *(v + j * n + k);
               *(v + j * n + k) = p;
               }
            }
         }
      return true;
}


bool geom_lsplane(double *coords, int *ids, int numIds, double *l, double *m, double *n, double *d )
// ***********************************************************************
{
        double  cent[3], eval[3], evec[3][3], mat[3][3], x, xsq, xy, xz,
                y, ysq, yz, z, zsq, sumx, sumy, sumz;
        int         nrot, i;
        double *cptr;
        double  u[3], v[3], ucv[3], mag;
        double *c1, *c2, *c3;

        *l = *m = *n = *d = 0.0;
        if ( !coords || !ids || numIds < 3 )
                return false;

        sumx = xsq = sumy = ysq = sumz = zsq = xy = xz = yz = 0.0;

        for ( i = 0; i < numIds; i++ )
        {
                cptr = coords + (ids[i] * 3);
                x = *cptr;
                y = *(cptr+1);
                z = *(cptr+2);

                sumx += x;
                xsq += x * x;

                sumy += y;
                ysq += y * y;

                sumz += z;
                zsq += z * z;

                xy += x * y;
                xz += x * z;
                yz += y * z;
        }
   cent[0] = sumx / (double) numIds;
   cent[1] = sumy / (double) numIds ;
   cent[2] = sumz / (double) numIds;

   mat[0][0] = xsq - sumx * cent[0];
   mat[0][1] = xy  - sumx * cent[1];
   mat[0][2] = xz  - sumx * cent[2];
   mat[1][0] = xy  - sumy * cent[0];
   mat[1][1] = ysq - sumy * cent[1];
   mat[1][2] = yz  - sumy * cent[2];
   mat[2][0] = xz  - sumz * cent[0];
   mat[2][1] = yz  - sumz * cent[1];
   mat[2][2] = zsq - sumz * cent[2];

   if (!UTL_GEOM_SYMM_EIGENSYS ((double *) mat, 3, eval, (double *) evec, &nrot))  return false;

   *l = evec[0][0];
   *m = evec[1][0];
   *n = evec[2][0];
   *d = (*l * cent[0] + *m * cent[1] + *n * cent[2]);

        c1 = coords + (ids[0] * 3);
        c2 = coords + (ids[1] * 3);
        c3 = coords + (ids[2] * 3);

        for ( i = 0; i < 3; i++ )
        {
                u[i] = *(c2+i) - *(c1+i);       /* 2 - 1 */
                v[i] = *(c3+i) - *(c2+i);       /* 3 - 2 */
        }
        UTL_GEOM_CROSS (u, v, ucv);

        mag = *d;
        mag += *l * ucv[0];
        mag += *m * ucv[1];
        mag += *n * ucv[2];
        if (mag < 0.0 )
        {
                *l *= -1.0;
                *m *= -1.0;
                *n *= -1.0;
                *d *= -1.0;
        }
   return true;
}

void UTL_GEOM_REFLECT (
// *********************************************************8
double l,
double m,
double n,
double d,
double *xyz)
{
   double h;

   h = l * xyz[0]  +  m * xyz[1]  +  n * xyz[2]  - d;
   xyz[0] -= 2.0 * l * h;
   xyz[1] -= 2.0 * m * h;
   xyz[2] -= 2.0 * n * h;

   return;
}

/*
int main() {
   double mat[9] = {2.23943,-.614402,0.0,-.61440,.551018,0.0,0.0,0.0,0.0};
   double eval[3];
   double evec[9];
   int nrot;
   int ans = UTL_GEOM_SYMM_EIGENSYS(mat,3,eval,evec,&nrot);
   cout << ans << endl;
   for (int i = 0; i < 3; ++i) cout << eval[i] << " ";
   cout << endl;
   for (int i = 0; i < 9; ++i) cout << evec[i] << " ";
   
}
*/
   
