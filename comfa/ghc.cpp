#include "openeye.h"

#include <fstream>
#include <iostream>
#include <cmath>
#include "oeplatform.h"
#include "oesystem.h"
#include "oechem.h"

#define Index2(x, y)     (x > y) ? (((x*(x-1))/2)+y) : (((y*(y-1))/2)+x)
#define sq_indx(x, y, z)  ((x*y)+z)

#define HUCK_THRESH             0.01
#define GAST_SIGMA_MAX_ITER     6
#define GAST_SIGMA_DAMP         0.5
#define SIGMA_ORB               1
#define HYDROGEN_TYPE           13
#define GAST_HYD_ELECT          20.02

using namespace std;
using namespace OEPlatform;
using namespace OEChem;
using namespace OESystem;

typedef struct _AtomSys {
    /* Used by charge calc routines    */
    int natoms;
    int *atoms;
} AtomSys;

typedef struct _GastOrb {
    int orb_type;
    /* SIGMA or PI orbital          */
    int type;
    /* atom type                    */
    double anion;
    /* electroneg - anion           */
    double neutral;
    /* electroneg - neutral         */
    double cation;
    /* electroneg - cation          */

    double coefA;
    /* parabolic parameters         */
    double coefB;
    double coefC;
} GastOrb;

typedef struct _HuckAtom {
    int at1;
    /* first atom type              */
    int bt1;
    /* first bond type              */
    int type;
    /* atom type                    */
    int bt2;
    /* second bond type             */
    int at2;
    /* second atom type             */
    double nel;
    /* number of donated electrons  */
    double h;              /* Coulomb integral             */
} HuckAtom;

typedef struct _HuckBond {
    int at1;
    /* first atom of bond           */
    int btype;
    /* bond type                    */
    int at2;
    /* second atom of bond          */
    double res_int;        /* resonance integral           */
} HuckBond;

typedef struct _HuckNBR {
    int num_neigh;
    /* #bound atoms on atom         */
    int type[10];
    /* type(s) of bound atom(s)     */
    int bond[10];       /* type(s) of bond(s)           */
} HuckNBR;


typedef struct _ChControl {
    int NumGastOrb;
    /* #orbital/type Gast. params   */
    GastOrb *Gast_orb;
    /* the actual information       */

    int NumHuckAtom;
    /* #atoms with Huckel params    */
    HuckAtom *Huck_atoms;
    /* the actual information       */
    int NumHuckBonds;
    /* #bonds with Huckel params    */
    HuckBond *Huck_bonds;   /* the actual information       */
} ChControl;

HuckAtom HA[35] =
        {{13, 1, 1,  1, 0,  2.0, 2.0},
         {0,  0, 2,  0, 0,  1.0, 0.0},
         {0,  0, 2,  2, 7,  1.0, 0.1},
         {0,  0, 2,  1, 14, 1.0, 0.2},
         {0,  0, 2,  1, 15, 1.0, 0.2},
         {0,  0, 2,  1, 16, 1.0, 0.1},
         {0,  0, 2,  1, 17, 1.0, 0.1},
         {0,  0, 3,  0, 0,  1.0, 0.0},
         {0,  0, 3,  1, 14, 1.0, 0.2},
         {0,  0, 3,  1, 15, 1.0, 0.2},
         {0,  0, 3,  1, 16, 1.0, 0.1},
         {0,  0, 3,  1, 17, 1.0, 0.1},
         {0,  0, 3,  5, 11, 1.0, 0.1},
         {0,  0, 4,  0, 0,  1.0, 0.0},
         {1,  1, 5,  1, 1,  2.0, 1.6},
         {13, 1, 5,  1, 1,  2.0, 1.7},
         {13, 1, 5,  1, 13, 2.0, 1.8},
         {0,  0, 31, 0, 0,  2.0, 0.3},
         {0,  0, 5,  0, 0,  2.0, 2.0},
         {0,  0, 6,  0, 0,  1.0, 0.4},
         {0,  0, 7,  0, 0,  1.0, 1.0},
         {0,  0, 8,  0, 0,  2.0, 2.0},
         {0,  0, 9,  0, 0,  1.0, 0.7},
         {0,  0, 10, 0, 0,  2.0, 2.0},
         {0,  0, 11, 0, 0,  1.0, 0.4},
         {0,  0, 12, 0, 0,  1.0, -1.0},
         {0,  0, 14, 0, 0,  2.0, 2.2},
         {0,  0, 15, 0, 0,  2.0, 2.4},
         {0,  0, 16, 0, 0,  2.0, 2.6},
         {0,  0, 17, 0, 0,  2.0, 2.0},
         {0,  0, 18, 0, 0,  1.0, 0.350},
         {0,  0, 19, 0, 0,  2.0, 1.6},
         {0,  0, 28, 0, 0,  2.0, 1.1},
         {0,  0, 32, 5, 0,  1.0, 0.7},
         {0,  5, 33, 5, 0,  0.0, -1.0}};

HuckBond HB[132] =
        {{1,  1, 2,  0.700},
         {1,  1, 3,  0.700},
         {2,  1, 2,  1.000},
         {2,  2, 2,  1.000},
         {2,  5, 2,  1.000},
         {2,  1, 3,  1.000},
         {2,  1, 4,  0.970},
         {2,  2, 4,  1.180},
         {2,  1, 31, 0.700},
         {2,  1, 5,  0.900},
         {2,  1, 6,  0.900},
         {2,  2, 6,  1.100},
         {2,  1, 7,  1.000},
         {2,  1, 8,  0.900},
         {2,  2, 9,  2.000},
         {2,  5, 32, 2.000},
         {2,  1, 10, 0.700},
         {2,  1, 12, 0.900},
         {2,  1, 14, 0.500},
         {2,  1, 15, 0.700},
         {2,  1, 16, 0.800},
         {2,  1, 17, 0.200},
         {2,  2, 18, 1.200},
         {2,  2, 19, 0.900},
         {2,  1, 19, 0.900},
         {2,  5, 19, 0.900},
         {33, 5, 19, 0.900},
         {2,  1, 28, 0.900},
         {2,  4, 28, 0.900},
         {3,  1, 3,  0.900},
         {3,  5, 3,  1.000},
         {3,  1, 4,  0.900},
         {3,  1, 31, 0.700},
         {3,  1, 5,  0.900},
         {3,  1, 6,  0.900},
         {3,  1, 7,  0.800},
         {3,  1, 8,  0.900},
         {3,  1, 10, 0.700},
         {3,  5, 10, 0.800},
         {3,  5, 11, 1.000},
         {3,  1, 12, 1.000},
         {3,  1, 14, 0.500},
         {3,  1, 15, 0.700},
         {3,  1, 16, 0.800},
         {3,  1, 17, 0.200},
         {3,  1, 19, 0.900},
         {3,  5, 19, 0.900},
         {3,  1, 28, 0.900},
         {3,  4, 28, 0.900},
         {4,  1, 4,  1.040},
         {4,  2, 4,  1.240},
         {4,  3, 4,  1.400},
         {4,  1, 5,  1.000},
         {4,  1, 6,  1.000},
         {4,  2, 6,  1.200},
         {4,  3, 7,  1.500},
         {4,  1, 8,  1.000},
         {4,  2, 9,  1.300},
         {4,  1, 10, 0.800},
         {4,  1, 12, 1.000},
         {4,  1, 14, 0.600},
         {4,  1, 15, 0.700},
         {4,  1, 16, 0.800},
         {4,  1, 17, 0.300},
         {4,  1, 19, 0.900},
         {4,  1, 28, 0.900},
         {4,  4, 28, 0.900},
         {31, 1, 6,  0.700},
         {5,  1, 6,  0.800},
         {5,  1, 7,  0.800},
         {6,  1, 6,  0.800},
         {6,  2, 6,  1.100},
         {6,  1, 7,  1.100},
         {6,  2, 7,  1.100},
         {6,  3, 7,  1.100},
         {6,  1, 8,  1.000},
         {6,  2, 9,  1.200},
         {6,  1, 10, 0.300},
         {6,  1, 11, 1.000},
         {6,  2, 11, 1.000},
         {6,  5, 11, 1.000},
         {6,  1, 12, 1.000},
         {6,  2, 18, 0.600},
         {6,  1, 19, 0.900},
         {6,  1, 28, 0.900},
         {6,  4, 28, 0.900},
         {7,  1, 7,  1.200},
         {7,  2, 7,  1.200},
         {7,  3, 7,  1.200},
         {7,  1, 8,  0.800},
         {7,  2, 9,  1.200},
         {7,  1, 10, 0.500},
         {7,  2, 18, 0.800},
         {8,  1, 8,  0.500},
         {8,  1, 10, 0.600},
         {8,  1, 11, 0.800},
         {8,  5, 11, 0.800},
         {8,  1, 12, 1.000},
         {8,  1, 14, 0.400},
         {8,  1, 15, 0.500},
         {8,  1, 16, 0.600},
         {8,  1, 17, 0.300},
         {8,  1, 19, 0.500},
         {8,  1, 28, 0.500},
         {9,  2, 10, 1.300},
         {9,  2, 12, 1.400},
         {32, 5, 12, 1.400},
         {9,  2, 14, 0.900},
         {9,  2, 15, 0.800},
         {9,  2, 16, 0.700},
         {9,  2, 17, 1.000},
         {9,  2, 18, 0.800},
         {10, 1, 10, 0.300},
         {10, 1, 11, 0.500},
         {10, 2, 18, 0.800},
         {10, 1, 19, 0.600},
         {10, 1, 28, 0.600},
         {11, 1, 19, 1.000},
         {11, 5, 19, 1.000},
         {11, 1, 28, 1.000},
         {12, 2, 18, 0.800},
         {12, 1, 19, 1.000},
         {12, 1, 28, 1.000},
         {3,  3, 6,  1.100},
         {3,  3, 8,  0.900},
         {3,  3, 10, 0.700},
         {6,  3, 8,  1.000},
         {6,  3, 10, 0.300},
         {8,  3, 10, 0.600},
         {5,  3, 6,  0.800},
         {6,  3, 6,  0.800},
         {5,  3, 7,  0.800}};

GastOrb GO[49] =
        {{1, 13, 0.37,  7.17,  12.85, 0, 0, 0},
         {1, 1,  0.68,  7.98,  19.04, 0, 0, 0},
         {1, 2,  0.98,  8.79,  19.62, 0, 0, 0},
         {1, 3,  0.98,  8.79,  19.62, 0, 0, 0},
         {1, 4,  1.67,  10.39, 20.57, 0, 0, 0},
         {1, 5,  2.08,  11.54, 23.72, 0, 0, 0},
         {1, 6,  2.57,  12.87, 24.87, 0, 0, 0},
         {1, 11, 2.57,  12.87, 24.87, 0, 0, 0},
         {1, 7,  3.71,  15.68, 27.11, 0, 0, 0},
         {1, 8,  2.65,  14.18, 28.49, 0, 0, 0},
         {1, 9,  3.75,  17.07, 31.33, 0, 0, 0},
         {1, 16, 3.12,  14.66, 30.82, 0, 0, 0},
         {1, 15, 2.66,  11.00, 22.04, 0, 0, 0},
         {1, 14, 2.77,  10.08, 19.71, 0, 0, 0},
         {1, 17, 2.90,  9.90,  18.82, 0, 0, 0},
         {1, 10, 2.39,  10.14, 20.65, 0, 0, 0},
         {2, 29, 0.00,  6.60,  20.64, 0, 0, 0},
         {2, 30, 0.00,  6.60,  20.64, 0, 0, 0},
         {1, 29, 2.39,  10.14, 20.65, 0, 0, 0},
         {1, 30, 2.39,  12.00, 24.00, 0, 0, 0},
         {2, 9,  1.23,  10.09, 24.69, 0, 0, 0},
         {2, 8,  0.00,  7.91,  29.52, 0, 0, 0},
         {2, 18, 1.38,  7.73,  17.70, 0, 0, 0},
         {2, 10, 0.00,  6.60,  20.64, 0, 0, 0},
         {2, 6,  0.89,  7.95,  20.35, 0, 0, 0},
         {2, 11, 0.89,  7.95,  20.35, 0, 0, 0},
         {2, 5,  0.00,  4.54,  23.72, 0, 0, 0},
         {2, 2,  -0.39, 5.60,  17.47, 0, 0, 0},
         {2, 3,  -0.39, 5.60,  17.47, 0, 0, 0},
         {2, 16, 3.12,  7.34,  30.84, 0, 0, 0},
         {2, 15, 2.30,  6.50,  21.68, 0, 0, 0},
         {2, 14, 0.00,  5.20,  19.36, 0, 0, 0},
         {2, 17, 0.00,  8.81,  17.62, 0, 0, 0},
         {2, 4,  0.05,  5.64,  16.85, 0, 0, 0},
         {2, 7,  0.83,  7.92,  20.38, 0, 0, 0},
         {1, 12, 1.62,  8.90,  18.10, 0, 0, 0},
         {1, 18, 2.72,  10.88, 21.69, 0, 0, 0},
         {1, 19, 2.46,  12.32, 24.86, 0, 0, 0},
         {2, 19, 0.00,  3.57,  20.34, 0, 0, 0},
         {1, 25, 1.06,  5.47,  11.65, 0, 0, 0},
         {1, 28, 2.46,  12.32, 24.86, 0, 0, 0},
         {2, 28, 0.00,  3.57,  20.35, 0, 0, 0},
         {1, 31, 0.00,  0.00,  23.72, 0, 0, 0},
         {2, 31, 0.00,  0.00,  23.72, 0, 0, 0},
         {1, 32, 3.75,  17.07, 31.33, 0, 0, 0},
         {2, 32, 1.23,  10.09, 24.69, 0, 0, 0},
         {1, 33, 0.98,  8.79,  19.62, 0, 0, 0},
         {2, 33, -0.39, 5.60,  17.47, 0, 0, 0}};

ChControl *ChCntl = NULL;

#define DET(a11, a12, a13, a21, a22, a23, a31, a32, a33) \
                         (a11*a22*a33 + a12*a23*a31 + a13*a21*a32 -\
                          a12*a21*a33 - a11*a23*a32 - a13*a22*a31)

void get_parabola(double x1, double y1, double x2, double y2,
                  double x3, double y3, double *a, double *b,
                  double *c) {
// *****************************************************
    double x1s, x2s, x3s, w;

    x1s = x1 * x1;
    x2s = x2 * x2;
    x3s = x3 * x3;

    w = DET(1.0, x1, x1s, 1.0, x2, x2s, 1.0, x3, x3s);

    if (fabs(w) < 0.0001) *a = *b = *c = 0.0;
    else {
        *a = DET(y1, x1, x1s, y2, x2, x2s, y3, x3, x3s) / w;
        *b = DET(1.0, y1, x1s, 1.0, y2, x2s, 1.0, y3, x3s) / w;
        *c = DET(1.0, x1, y1, 1.0, x2, y2, 1.0, x3, y3) / w;
    }
}

void ChgInit() {
// ****************************************************

    ChCntl = new _ChControl[1];
    ChCntl->NumGastOrb = 49;
    ChCntl->Gast_orb = GO;
    ChCntl->NumHuckAtom = 35;
    ChCntl->Huck_atoms = HA;
    for (int i = 0; i < 35; i++) ChCntl->Huck_atoms[i].h = -ChCntl->Huck_atoms[i].h;  //weird!!
    ChCntl->NumHuckBonds = 132;
    ChCntl->Huck_bonds = HB;
    for (int i = 0; i < 123; i++) ChCntl->Huck_bonds[i].res_int = -ChCntl->Huck_bonds[i].res_int;  //weird!!
    for (int indx = 0; indx < ChCntl->NumGastOrb; indx++)
        get_parabola(-1.0, ChCntl->Gast_orb[indx].anion,
                     0.0, ChCntl->Gast_orb[indx].neutral,
                     1.0, ChCntl->Gast_orb[indx].cation,
                     &ChCntl->Gast_orb[indx].coefA,
                     &ChCntl->Gast_orb[indx].coefB,
                     &ChCntl->Gast_orb[indx].coefC);
}

bool get_ATypes(const OEGraphMol mol, int *SybAtypes) {
// *******************************************************
    for (OEIter<OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        unsigned int anum = atom->GetAtomicNum();
        int t = 0;
        int ocount = 0;
        unsigned int hyb = atom->GetDegree();

        switch (anum) {
            case 1:
                t = 13;
                break;
            case 6:
                if (atom->IsAromatic()) t = 3;
                else {
                    if (hyb == 4) t = 1;
                    if (hyb == 3) t = 2;
                    if (hyb == 2) t = 4;
                }
                break;
            case 7:
                if (atom->IsAromatic() && (OEAtomIsInRingSize(*atom, 6) || atom->GetExplicitDegree() == 2)) t = 11;
                else {
                    if (hyb == 4) t = 31;
                    if (hyb == 3) {
                        if (atom->GetHvyDegree() == 4) t = 31;
                        else {
                            for (OEIter<OEAtomBase> nbor = atom->GetAtoms(); nbor; ++nbor)
                                if (nbor->GetDegree() != 4 && nbor->GetAtomicNum() == 6) {
                                    t = 19;
                                    for (OEIter<OEAtomBase> nbor2 = nbor->GetAtoms(); nbor2; ++nbor2)
                                        if (nbor2->GetAtomicNum() == 8 || nbor2->GetAtomicNum() == 16)
                                            t = 28;
                                }
                            if (t == 0) t = 5;
                        }
                    }
                    if (hyb == 2) t = 6;
                    if (hyb == 1) t = 7;
                }
                break;
            case 8:
                if (hyb == 2) t = 8;
                else t = 9;
                break;
            case 16:
                for (OEIter<OEAtomBase> nbor = atom->GetAtoms(); nbor; ++nbor)
                    if (nbor->GetAtomicNum() == 8) ocount++;
                if (atom->GetDegree() == 3 && ocount == 1) t = 29;
                if (atom->GetDegree() == 4 && ocount == 2) t = 30;
                if (atom->GetDegree() == 2) t = 10;
                if (atom->GetDegree() == 1) t = 18;
                if (t == 0) t = 18;
                break;
            case 15:
                t = 12;
                break;
            case 35:
                t = 14;
                break;
            case 17:
                t = 15;
                break;
            case 9:
                t = 16;
                break;
            case 53:
                t = 17;
                break;
            case 11:
                t = 21;
                break;
            case 19:
                t = 22;
                break;
            case 20:
                t = 23;
                break;
            case 3:
                t = 24;
                break;
            case 13:
                t = 25;
                break;
            case 14:
                t = 27;
                break;
        }
        if (t == 0) {
            cout << "WARNING: No Gast-Huck parameters for atomic number " << anum << '\n';
            return false;
        }
        SybAtypes[atom->GetIdx()] = t;
    }
    return true;
}

double getRI(int type1, int btype, int type2) {
// **************************************
    double RI = 0.0;
    for (int indx = 0; indx < ChCntl->NumHuckBonds; indx++)
        if (((type1 == ChCntl->Huck_bonds[indx].at1) &&
             (btype == ChCntl->Huck_bonds[indx].btype) &&
             (type2 == ChCntl->Huck_bonds[indx].at2)) ||
            ((type2 == ChCntl->Huck_bonds[indx].at1) &&
             (btype == ChCntl->Huck_bonds[indx].btype) &&
             (type1 == ChCntl->Huck_bonds[indx].at2)))
            RI = ChCntl->Huck_bonds[indx].res_int;
    return RI;
}

void getResonIntgl(const OEGraphMol mol, const int *ATypes, const int *BTypes,
                   double *ResonIntgl) {
// ******************************************************

    for (int i = 0; i < ((int) mol.NumAtoms() * ((int) mol.NumAtoms() - 1)) / 2; ++i)
        ResonIntgl[i] = 0.0;
    for (OEIter<OEBondBase> bond = mol.GetBonds(); bond; ++bond) {
        int btype = BTypes[bond->GetIdx()];
        int at1 = (bond->GetBgn())->GetIdx();
        int type1 = ATypes[at1];
        int at2 = (bond->GetEnd())->GetIdx();
        int type2 = ATypes[at2];
        ResonIntgl[Index2(at1, at2)] = getRI(type1, btype, type2);
    }
}

bool diag_hamil(int n, double *a, double *v, double *d) {
// *************************************************************
    double c, g, h, p, s, sm, t, tau, theta, thresh;
    int i, ip, iq, j, k;

    for (ip = 0; ip < n; ip++)
        for (iq = 0; iq < n; iq++)
            *(v + ip * n + iq) = (ip == iq) ? 1.0 : 0.0;
    for (ip = 0; ip < n; ip++)
        *(d + ip) = *(a + ip * n + ip);

    for (i = 0; i < 50; i++) {
        sm = 0.0;
        for (ip = 0; ip < n - 1; ip++)
            for (iq = ip + 1; iq < n; iq++)
                sm += fabs(*(a + ip * n + iq));

        if (sm == 0.0) goto done;

        thresh = (i < 3) ? (0.2 * sm / ((double) (n * n))) : 0.0;
        for (ip = 0; ip < n - 1; ip++) {
            for (iq = ip + 1; iq < n; iq++) {
                g = 100.0 * fabs(*(a + ip * n + iq));

                if ((i > 4) && (fabs(*(d + ip)) + g == fabs(*(d + ip)))
                    && (fabs(*(d + iq)) + g == fabs(*(d + iq))))
                    *(a + ip * n + iq) = 0.0;
                else if (fabs(*(a + ip * n + iq)) > thresh) {
                    h = *(d + iq) - *(d + ip);
                    if (fabs(h) + g == fabs(h))
                        t = *(a + ip * n + iq) / h;  /* t = 1 / (2 * theta) */
                    else {
                        theta = 0.5 * h / *(a + ip * n + iq);  /* Eq 11.1.10 */
                        t = 1.0 / (fabs(theta) + sqrt(1.0 + theta * theta));
                        if (theta < 0.0) t = -t;
                    }
                    c = 1.0 / sqrt(1.0 + t * t);
                    s = t * c;
                    tau = s / (1.0 + c);
                    h = t * (*(a + ip * n + iq));
                    *(d + ip) -= h;
                    *(d + iq) += h;
                    *(a + ip * n + iq) = 0.0;

                    for (j = 0; j < ip; j++) {
                        g = *(a + j * n + ip);
                        h = *(a + j * n + iq);
                        *(a + j * n + ip) = g - s * (h + g * tau);
                        *(a + j * n + iq) = h + s * (g - h * tau);
                    }
                    for (j = ip + 1; j < iq; j++) {
                        g = *(a + ip * n + j);
                        h = *(a + j * n + iq);
                        *(a + ip * n + j) = g - s * (h + g * tau);
                        *(a + j * n + iq) = h + s * (g - h * tau);
                    }
                    for (j = iq + 1; j < n; j++) {
                        g = *(a + ip * n + j);
                        h = *(a + iq * n + j);
                        *(a + ip * n + j) = g - s * (h + g * tau);
                        *(a + iq * n + j) = h + s * (g - h * tau);
                    }
                    for (j = 0; j < n; j++) {
                        g = *(v + j * n + ip);
                        h = *(v + j * n + iq);
                        *(v + j * n + ip) = g - s * (h + g * tau);
                        *(v + j * n + iq) = h + s * (g - h * tau);
                    }
                }
            }
        }
    }
    return false;

    done:
    for (ip = 0; ip < n - 1; ip++) {
        for (iq = ip + 1; iq < n; iq++)
            *(a + ip * n + iq) = *(a + iq * n + ip);
    }
    for (i = 0; i < n - 1; i++) {
        k = i;
        p = *(d + i);
        for (j = i + 1; j < n; j++) {
            if (*(d + j) <= p) {
                k = j;
                p = *(d + j);
            }
        }
        if (k != i) {
            *(d + k) = *(d + i);
            *(d + i) = p;
            for (j = 0; j < n; j++) {
                p = *(v + j * n + i);
                *(v + j * n + i) = *(v + j * n + k);
                *(v + j * n + k) = p;
            }
        }
    }
    return true;
}

void recordHuck(const OEGraphMol mol, int *HAggs, bool *seen, int *npi, const OEAtomBase *nxt,
                bool inHuck, const double *ResInt, const int *ATypes, const int *BTypes, const bool *Hcands) {
// ***********************************************************
    std::vector<OEAtomBase *> nbAts;
    std::vector<OEBondBase *> nbBds;
    std::vector<bool> Huckbd;
    for (OEIter <OEBondBase> bond = nxt->GetBonds(); bond; ++bond) {
        if (seen[bond->GetIdx()]) continue;
        OEAtomBase *nbor = bond->GetNbr(nxt);
        nbAts.push_back(nbor);
        nbBds.push_back(bond);
        bool isHuck = ResInt[Index2(nbor->GetIdx(), nxt->GetIdx())] != 0.0;
        Huckbd.push_back(isHuck);
    }
    int nNbrs = (int) nbAts.size();
    if (nNbrs == 0) return;

    int *NbrOrd = new int[nNbrs];
    int ino = 0;
    for (int i = 0; i < nNbrs; i++)
        if (Huckbd[i]) {
            NbrOrd[ino] = i;
            ino++;
        }
    for (int i = 0; i < nNbrs; i++)
        if (!Huckbd[i]) {
            NbrOrd[ino] = i;
            ino++;
        }

    bool nowHuck = inHuck;
    for (int nxAt = 0; nxAt < nNbrs; ++nxAt) {
        OEChem::OEAtomBase *nbAt = nbAts[NbrOrd[nxAt]];
//        int ii = nxt->GetIdx();
//        int io = nbAt->GetIdx();
        OEChem::OEBondBase *nbBd = nbBds[NbrOrd[nxAt]];
//        cout << nbBd->GetIdx() << endl;
        if (Huckbd[NbrOrd[nxAt]]) {
            if (!inHuck) {
                if (!nowHuck)
                    (*npi)++; // if the 1st encounter with this HAgg
                HAggs[nbAt->GetIdx()] = *npi;
            }
            else HAggs[nbAt->GetIdx()] = HAggs[nxt->GetIdx()];
            nowHuck = true;
        }
        else
            // if (!inHuck || (ResInt[ Index2(nxt->GetIdx(),nbAt->GetIdx()) ] == 0.0))
            //               else if (!inHuck || !nxt->IsInRing() || ! nbAt->IsInRing() || nbAt->GetExplicitDegree() < 3)
//        else if (!inHuck || getRI(ATypes[ nbAt->GetIdx()], BTypes[ nbBd->GetIdx() ], ATypes[nxt->GetIdx()]) == 0.0)
            nowHuck = false;

        seen[nbBd->GetIdx()] = true;
        if (nowHuck) recordHuck(mol, HAggs, seen, npi, nbAt, nowHuck, ResInt, ATypes, BTypes, Hcands);
    }
}

bool getCoulArray(const OEGraphMol mol, const int *ATypes, const int *BTypes,
                  double *Coulomb, double *Nel, double *ResInt, double *charges) {
// ******************************************************

// types of all neighbor atoms & bonds to each atom
    int *atomIdx = new int[(int) mol.NumAtoms()];
    for (int i = 0; i < (int) mol.NumAtoms(); i++) {
        atomIdx[i] = 0;
        Coulomb[i] = Nel[i] = 0.0;
    }

    double *occ = new double[((int) mol.NumAtoms()) * 8];
    HuckNBR *Huck_nbr = new _HuckNBR[((int) mol.NumAtoms())];

    for (OEIter<OEBondBase> bond = mol.GetBonds(); bond; ++bond) {
        int at1 = bond->GetBgnIdx();
        int at2 = bond->GetEndIdx();
        Huck_nbr[at1].type[atomIdx[at1]] = ATypes[at2];
        Huck_nbr[at1].bond[atomIdx[at1]] = BTypes[bond->GetIdx()];
        atomIdx[at1] += 1;

        Huck_nbr[at2].type[atomIdx[at2]] = ATypes[at1];
        Huck_nbr[at2].bond[atomIdx[at2]] = BTypes[bond->GetIdx()];
        atomIdx[at2] += 1;
    }
    for (int indx = 0; indx < (int) mol.NumAtoms(); indx++)
        Huck_nbr[indx].num_neigh = atomIdx[indx];

    for (int iat = 0; iat < (int) mol.NumAtoms(); iat++) {
        int best_match = 0;
        for (int indx = 0; indx < ChCntl->NumHuckAtom; indx++) {
            int right_match = -1;
            int left_match = -1;
            for (int indx2 = 0; indx2 < Huck_nbr[iat].num_neigh && right_match != 2; indx2++)
                if (ChCntl->Huck_atoms[indx].type == ATypes[iat]) {
                    if ((ChCntl->Huck_atoms[indx].at1 == Huck_nbr[iat].type[indx2]) &&
                        (ChCntl->Huck_atoms[indx].bt1 == Huck_nbr[iat].bond[indx2]))
                        right_match = 2;
                    else if (((ChCntl->Huck_atoms[indx].at1 == Huck_nbr[iat].type[indx2]) &&
                              (!ChCntl->Huck_atoms[indx].bt1)) ||
                             ((ChCntl->Huck_atoms[indx].bt1 == Huck_nbr[iat].bond[indx2]) &&
                              (!ChCntl->Huck_atoms[indx].at1))) {
                        if (right_match < 1) right_match = 1;
                    }
                    else if (!ChCntl->Huck_atoms[indx].at1 &&
                             (!ChCntl->Huck_atoms[indx].bt1)) {
                        if (right_match < 0) right_match = 0;
                    }
                }
            if (right_match == -1) continue;
            for (int indx2 = 0; indx2 < Huck_nbr[iat].num_neigh && left_match != 2; indx2++) {
                if ((ChCntl->Huck_atoms[indx].at2 == Huck_nbr[iat].type[indx2]) &&
                    (ChCntl->Huck_atoms[indx].bt2 == Huck_nbr[iat].bond[indx2]))
                    left_match = 2;

                else if (((ChCntl->Huck_atoms[indx].at2 == Huck_nbr[iat].type[indx2]) &&
                          (!ChCntl->Huck_atoms[indx].bt2)) ||
                         ((ChCntl->Huck_atoms[indx].bt2 == Huck_nbr[iat].bond[indx2]) &&
                          (!ChCntl->Huck_atoms[indx].at2))) {
                    if (left_match < 1) left_match = 1;
                }
                else if (!ChCntl->Huck_atoms[indx].at2 &&
                         (!ChCntl->Huck_atoms[indx].bt2)) {
                    if (left_match < 0) left_match = 0;
                }
            }
            if (left_match == -1) continue;
            if ((1 + right_match + left_match) > best_match) {
                Coulomb[iat] = ChCntl->Huck_atoms[indx].h;
                Nel[iat] = ChCntl->Huck_atoms[indx].nel;
                best_match = 1 + right_match + left_match;
            }
        }
    }
    int *HAggs = new int[mol.NumAtoms()];
    for (unsigned int i = 0; i < mol.NumAtoms(); i++) HAggs[i] = 0;
    bool *seen = new bool[mol.NumBonds()];
    for (unsigned int i = 0; i < mol.NumBonds(); i++) seen[i] = false;
    bool *Hcands = new bool[mol.NumAtoms()];
    for (unsigned int i = 0; i < mol.NumAtoms(); i++) Hcands[i] = false;
    for (OEIter<OEBondBase> bi = mol.GetBonds(); bi; ++bi) {
        int at1 = bi->GetBgn()->GetIdx();
        int at2 = bi->GetEnd()->GetIdx();
        if (ResInt[Index2(at1, at2)] != 0.0)
            Hcands[at1] = Hcands[at2] = true;
    }
    int nPiSys = 0;
    HAggs[0] = 1;

    for (OEIter<OEBondBase> bi = mol.GetBonds(); bi; ++bi)
        if (!seen[bi->GetIdx()]) {
            int at1 = bi->GetBgn()->GetIdx();
            int at2 = bi->GetEnd()->GetIdx();
            if (ResInt[Index2(at1, at2)] == 0.0) {
                seen[bi->GetIdx()] = true;
                continue;
            }
            recordHuck(mol, HAggs, seen, &nPiSys, bi->GetBgn(), false, ResInt, ATypes, BTypes, Hcands);

        }

    AtomSys *pi_sys = new AtomSys[nPiSys];
    for (int indx = 0; indx < nPiSys; ++indx) {
        int natoms = 0;
        for (unsigned int i = 0; i < mol.NumAtoms(); i++) if (HAggs[i] == indx + 1) natoms++;
        pi_sys[indx].natoms = natoms;
        pi_sys[indx].atoms = new int[natoms];
        int j = 0;
        for (unsigned int i = 0; i < mol.NumAtoms(); i++)
            if (HAggs[i] == indx + 1) {
                pi_sys[indx].atoms[j] = i;
                j++;
            }
    }

    for (int indx = 0; indx < nPiSys; ++indx) {
        int natoms = pi_sys[indx].natoms;
        double *hmat = new double[pi_sys[indx].natoms * pi_sys[indx].natoms];
        for (int i = 0; i < pi_sys[indx].natoms * pi_sys[indx].natoms; i++)
            hmat[i] = 0.0;
        double *coef = new double[pi_sys[indx].natoms * pi_sys[indx].natoms];
        for (int i = 0; i < pi_sys[indx].natoms * pi_sys[indx].natoms; i++)
            coef[i] = 0.0;
        double *eig = new double[pi_sys[indx].natoms * 8];
        for (int i = 0; i < pi_sys[indx].natoms * 8; i++)
            eig[i] = 0.0;
        for (int indx2 = 0; indx2 < pi_sys[indx].natoms; indx2++)
            hmat[sq_indx(pi_sys[indx].natoms, indx2, indx2)] =
                    Coulomb[pi_sys[indx].atoms[indx2]];
        for (int indx2 = 0; indx2 < pi_sys[indx].natoms; indx2++) {
            int iat = pi_sys[indx].atoms[indx2];
            int jat = 0;
            for (int indx3 = 0; indx3 < pi_sys[indx].natoms; indx3++)
                if (indx3 != indx2) {
                    jat = pi_sys[indx].atoms[indx3];

                    if (ResInt[Index2(iat, jat)] != 0.0)
                        hmat[sq_indx(pi_sys[indx].natoms, indx2, indx3)] =
                        hmat[sq_indx(pi_sys[indx].natoms, indx3, indx2)] =
                                ResInt[Index2(iat, jat)];
                }
        }
        if (!diag_hamil(pi_sys[indx].natoms, hmat, coef, eig)) return false;

        double nelect = 0.0;
        for (int indx2 = 0; indx2 < pi_sys[indx].natoms; indx2++)
            nelect += Nel[pi_sys[indx].atoms[indx2]] -
                      charges[pi_sys[indx].atoms[indx2]];
        if (nelect > (2.0 * (double) pi_sys[indx].natoms)) {
            cout << "Too many e- in pi calculation\n";
            return false;
        }

        for (int i = 0; i < (int) mol.NumAtoms() * 8; i++) occ[i] = 0.0;
        int norb = 0;
        int num_left = (int) nelect;
        int ndegen = 0;
        do {
            /* Check for degeneracy */

            for (ndegen = 1;
                /* (fabs(eig[norb+ndegen]-eig[norb]) < HUCK_THRESH); */
                 (fabs(eig[norb + ndegen] - eig[norb + ndegen - 1]) < HUCK_THRESH) && (norb + ndegen) < natoms;
                 ndegen++);

            if (num_left >= ((double) ndegen * 2.0))
                for (int indx2 = 0; indx2 < ndegen; indx2++) {
                    occ[norb++] = 2.0;
//                    num_left -= 2.0;
                    num_left -= 2;
                }
            else {    /* have to partition e- over orbitals */
                double q = (double) num_left / (double) ndegen;
                for (int indx2 = 0; indx2 < ndegen; indx2++) {
                    occ[norb++] = q;
                    num_left -= int(q);
                }
            }
        } while (num_left > 0.0 && norb < (int) mol.NumAtoms());
        for (int indx2 = 0; indx2 < pi_sys[indx].natoms; indx2++) {
            double q = 0.0;

            for (int indx3 = 0; indx3 < pi_sys[indx].natoms; indx3++) {
                double c = coef[sq_indx(pi_sys[indx].natoms, indx2, indx3)];
                q += occ[indx3] * c * c;
            }
            charges[pi_sys[indx].atoms[indx2]] = Nel[pi_sys[indx].atoms[indx2]] - q;
        }
    }
    return true;
}

bool getGastSigma(const OEGraphMol mol, int *ATypes, double *AtChg) {
// *******************************************************
    int maxType = 0;
    for (int indx = 0; indx < ChCntl->NumGastOrb; indx++)
        if ((ChCntl->Gast_orb[indx].orb_type == SIGMA_ORB) &&
            (ChCntl->Gast_orb[indx].type > maxType))
            maxType = ChCntl->Gast_orb[indx].type;
    int gastSigmaMaxType = maxType;
    bool *gastSigmaExists = new bool[maxType + 1];
    double *gastSigmaCoefA = new double[maxType + 1];
    double *gastSigmaCoefB = new double[maxType + 1];
    double *gastSigmaCoefC = new double[maxType + 1];

    for (int indx = 0; indx < ChCntl->NumGastOrb; indx++)
        if (ChCntl->Gast_orb[indx].orb_type == SIGMA_ORB) {
            int idx2 = ChCntl->Gast_orb[indx].type;
            gastSigmaExists[idx2] = true;
            gastSigmaCoefA[idx2] = ChCntl->Gast_orb[indx].coefA;
            gastSigmaCoefB[idx2] = ChCntl->Gast_orb[indx].coefB;
            gastSigmaCoefC[idx2] = ChCntl->Gast_orb[indx].coefC;
        }
    double *Chi = new double[((int) mol.NumAtoms())];
    double *savChg = new double[((int) mol.NumAtoms())];
    for (int i = 0; i < ((int) mol.NumAtoms()); i++) savChg[i] = AtChg[i];

    double current_damp = 1.0;
    for (int iter = 0; iter < GAST_SIGMA_MAX_ITER; iter++) {
        current_damp *= GAST_SIGMA_DAMP;
        for (int i = 0; i < ((int) mol.NumAtoms()); i++) Chi[i] = 0.0;

        for (int indx = 0; indx < ((int) mol.NumAtoms()); indx++) {
            int aType1 = ATypes[indx];
            if (aType1 >= 0 && aType1 <= gastSigmaMaxType && gastSigmaExists[aType1]) {
                double coefa = gastSigmaCoefA[aType1];
                double coefb = gastSigmaCoefB[aType1];
                double coefc = gastSigmaCoefC[aType1];
                double chrg = AtChg[indx];
                Chi[indx] = coefa + (coefb * chrg) + (coefc * chrg * chrg);
            }
        }
        for (OEIter<OEBondBase> bond = mol.GetBonds(); bond; ++bond) {
            int aType1 = ATypes[bond->GetBgnIdx()];
            int aType2 = ATypes[bond->GetEndIdx()];
            if (aType1 < 0 || aType1 > gastSigmaMaxType || !gastSigmaExists[aType1] ||
                aType2 < 0 || aType2 > gastSigmaMaxType || !gastSigmaExists[aType2])
                continue;
            double coefa = gastSigmaCoefA[aType1];
            double coefb = gastSigmaCoefB[aType1];
            double coefc = gastSigmaCoefC[aType1];
            double coefa2 = gastSigmaCoefA[aType2];
            double coefb2 = gastSigmaCoefB[aType2];
            double coefc2 = gastSigmaCoefC[aType2];

            double delta_electroneg = Chi[bond->GetBgnIdx()] - Chi[bond->GetEndIdx()];
            double cat_electroneg = 0.0;
            int atom_type = 0;
            if (delta_electroneg < 0.0) {
                cat_electroneg = coefa + coefb + coefc;
                atom_type = aType1;
            }
            else {
                cat_electroneg = coefa2 + coefb2 + coefc2;
                atom_type = aType2;
            }
            if (atom_type == HYDROGEN_TYPE) cat_electroneg = GAST_HYD_ELECT;
            double damped_charge = current_damp * delta_electroneg / cat_electroneg;
            AtChg[bond->GetBgnIdx()] -= damped_charge;
            AtChg[bond->GetEndIdx()] += damped_charge;
        }
    }
    return true;
}

void getFormalCharges(const OEGraphMol mol, const int *ATypes, double *AtChg) {
// **********************************************
    static int type_valences[] = {-1, 4, 3, 3, 2, 4, 3, 2, 4, 3, 4, 3, 4, 1, 4, 4, 4, 4, 3, 3, 1, 0, 0, 0, 0, 3, 0, 4,
                                  3, 4, 4, 4, 3, 3};
    static int type_lonep[] = {-1, 0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 1, 0, 0, 3, 3, 3, 3, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                               1, 0, 0, 2, 0};
    static double def_charge[] = {0.0, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
                                  0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.00, 1.00, 2.00, 1.00, 0.00, 0.00, 0.00,
                                  0.00, 0.00, 0.00, 1.00, -0.50, 1.00};
    for (OEIter<OEAtomBase> ai = mol.GetAtoms(); ai; ++ai) {
        int idx = ATypes[ai->GetIdx()];
        double cCharge = def_charge[idx] - (double) (type_valences[idx] - type_lonep[idx]);
        cCharge += (double) ai->GetExplicitDegree();
        AtChg[ai->GetIdx()] = cCharge;
    }
}

bool GHCharges(const OEGraphMol mol, double *AtChg) {
// *******************************************************

    if (ChCntl == NULL) ChgInit();

    int *ATypes = new int[(int) mol.NumAtoms()];
    if (!get_ATypes(mol, ATypes)) return false;
    int *BTypes = new int[(int) mol.NumBonds()];
    for (OEIter<OEBondBase> bond = mol.GetBonds(); bond; ++bond) {
        if (bond->GetIntType() == 1 &&
            ((ATypes[bond->GetBgnIdx()] == 28 && ATypes[bond->GetEndIdx()] == 1) ||
             (ATypes[bond->GetEndIdx()] == 28 && ATypes[bond->GetBgnIdx()] == 1)))
            BTypes[bond->GetIdx()] = 4; // amide
        else {
            if (bond->IsAromatic()) BTypes[bond->GetIdx()] = 5;
            else BTypes[bond->GetIdx()] = bond->GetIntType();
        }
    }

    getFormalCharges(mol, ATypes, AtChg);

    double *ri = new double[((int) mol.NumAtoms() * ((int) mol.NumAtoms() - 1)) / 2];
    getResonIntgl(mol, ATypes, BTypes, ri);

    double *Coulomb = new double[((int) mol.NumAtoms())];
    double *Nel = new double[((int) mol.NumAtoms())];
    if (!getCoulArray(mol, ATypes, BTypes, Coulomb, Nel, ri, AtChg)) return false;
    if (!getGastSigma(mol, ATypes, AtChg)) return false;

    return true;
}
/*
int main() {
  oemolistream ifs("/Users/richardcramer/redo/cdk2/gh.sdf");
  OEGraphMol mol;
  while (OEReadMolecule(ifs,mol)) {
      OEAddExplicitHydrogens(mol);
	double *chg = new double[ (int) mol.NumAtoms() ];
	if (!GHCharges(mol,chg)) return 1;
	delete chg;
  }
    return 0;
}
*/
