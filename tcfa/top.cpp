/*
#include "openeye.h"

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include "string.h"
#include "oeplatform.h"
#include "oesystem.h"
#include "oechem.h"
#include "oedepict.h"

using namespace std;
using namespace OEPlatform;
using namespace OEChem;
using namespace OESystem;
using namespace OEDepict;
using namespace OEMath;

#include "../plane/plane.cpp"
unsigned int Crings[200];
*/

typedef struct top_priority_def {
        int atomIdx;
        int covID;
        int branchCount;
        double branchWeight;
        double branchDepthWeight;
} TopPriority;


void getBondDist( const OEGraphMol mol, unsigned int ABase, int *BDist,
	bool *AXclude ) {
// ***********************************************************

    for (unsigned int i = 0; i < mol.NumAtoms(); ++i) BDist[i] = 0;
    BDist[ ABase ] = 1;
    if (AXclude) for (unsigned int i = 0; i < mol.NumAtoms(); ++i)
	if (AXclude[i]) BDist[i] = -1; 
    
    bool added = true;
    for (int level = 1; added && level <= (int) mol.NumAtoms(); ++level) {
	added = false;
        for (OEIter<OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) 
	    if (BDist[atom->GetIdx()] == level) 
		for (OEIter<OEAtomBase> nbor = atom->GetAtoms(); nbor; ++nbor) {
		    if (abs(BDist[nbor->GetIdx()]) > 0) continue;
		    BDist[nbor->GetIdx()] = level + 1;
		    added = true;
		}
    }
}

void RingPuckerLR( OEGraphMol mol, double *coo, int a1, int a2, int a3,
	bool *detecting ) {
// ***********************************************************

// cout << a1 << " " << a2 << " " << a3 <<endl;
    if (detecting) *detecting = false;
    	else if (Crings[a2] + Crings[a3] == 0 || Crings[a2] != Crings[a3] ) return; // shouldm't happen !
    vector<int> rids;
    for (unsigned int i = 0; i < mol.NumAtoms(); ++i) if (Crings[i] == Crings[a1])
	rids.push_back(i);
    if ((int) rids.size() == 3) return; // 3 membered rings can't pucker
    OEAtomBase *a2pr = mol.GetAtom(OEHasAtomIdx(a2));    
    OEAtomBase *a3pr = mol.GetAtom(OEHasAtomIdx(a3));    
    int planeAt = -1; 
//    int a4 = -1;
    if (detecting) planeAt = a1;
    else {
     int nfd = 0;
     for (int nra = 0; nra < (int) rids.size(); ++nra) {
	if (a2pr->GetBond( mol.GetAtom(OEHasAtomIdx(rids[nra])))) {
	    ++nfd;
	    planeAt = nra;
	}
	if (a3pr->GetBond( mol.GetAtom(OEHasAtomIdx(rids[nra])))) {
	    ++nfd;
//	    a4 = nra;
	}
	if (nfd == 2) break;
     }
    }
    int ids[3];
    ids[0] = planeAt;
    ids[1] = a2;
    ids[2] = a3;
    if (planeAt < a3) {
	ids[0] = a3;
	ids[2] = planeAt;
    }
    double l, m, n, d, center[3], x, y, z, h;
    if (!geom_lsplane(coo, ids, 3, &l, &m, &n, &d)) return;
    for (int i = 0; i < 3; i++) center[i] = 0.0;
    double testh = 0.0; 
    double toth = 0.0; 
//    double patomh = 0.0; 
    double *cptr;
    for (int nra = 0; nra < (int) rids.size(); ++nra) {
    	cptr = coo + rids[nra] * 3;
	x = *cptr; y = *(cptr+1); z = *(cptr+2);
	h = l*x + m*y + n*z - d;
	if (rids[nra] == a3 || rids[nra] == a2 || rids[nra] == planeAt) {
	    testh += h;
	}
//	if (rids[nra] == a4) patomh = h;
	toth += h;
	center[0] += x * h;
	center[1] += y * h;
	center[2] += z * h;
    }
// cout << planeAt << " " << a2 << " " << a3 << " " << toth << " " << (int) rids.size() << endl;
    if (fabs(toth) < 0.4) return; // planar ring system, so nothing to do
    if (detecting) {
	*detecting = true;
	return;
    }
    for (int i = 0; i < 3; i++) center[i] /= toth;
    double *c1, *c2, *c3;
    c1 = coo+a1*3; c2 = coo+a2*3; c3 = coo+a3*3;
    double dihed = OEGeom3DTorsion(c1,c2,c3,center);
    while (dihed < 0.0) dihed += 6.2831853072;
    while (dihed > 6.2831853072) dihed -= 6.2831853072;
    if (dihed <= 3.1415926536) return; // ring bulk already is local right

    ids[0] = a1; ids[1] = a2; ids[2] = a3;
    if (!geom_lsplane(coo, ids, 3, &l, &m, &n, &d)) return;
    OEAtomBase *stem = mol.GetAtom(OEHasAtomIdx(a2));
    OEAtomBase *root = mol.GetAtom(OEHasAtomIdx(a3));
    for (OEIter<OEAtomBase> atom = OEGetSubtree(stem,root); atom; ++atom)
	UTL_GEOM_REFLECT( l, m, n, d, &coo[ 3*atom->GetIdx() ]);
}

void reflectAtoms( OEGraphMol mol, double *coo, int npt, int *aplane, vector<int> FrameAts, bool doAll,
	bool simple)
/* =================================================================================== */
/* reflects atms through the plane defined by the atoms whose indexes (base 0 )are in aplane, by modifying values in coo */
{

double cent[3], eval[3], evec[3][3], mat[3][3], x, xsq, xy, xz,
               y, ysq, yz, z, zsq, *cx, *cy, *cz, l, m, n, d, *xyz, h;
int na, nrot, elem;

//cout << "Ref";
/* Now perform the sums to determine the parameters of the plane      */
/* equation.                                                          */
   x = xsq = y = ysq = z = zsq = xy = xz = yz = 0.0;
   for (na = 0; na < npt; na++ ) {
      cx = coo + 3 * ( aplane[ na ] );
      x   += *cx;
      xsq += (*cx) * (*cx);
      cy = cx + 1;
      y   += *cy;
      ysq += (*cy) * (*cy);
      cz = cy + 1;
      z   += *cz;
      zsq += (*cz) * (*cz);
      xy  += (*cx) * (*cy);
      xz  += (*cx) * (*cz);
      yz  += (*cy) * (*cz);
   }
   cent[0] = x / (double) npt;
   cent[1] = y / (double) npt;
   cent[2] = z / (double) npt;

   mat[0][0] = xsq - x * cent[0];
   mat[0][1] = xy  - x * cent[1];
   mat[0][2] = xz  - x * cent[2];
   mat[1][0] = xy  - y * cent[0];
   mat[1][1] = ysq - y * cent[1];
   mat[1][2] = yz  - y * cent[2];
   mat[2][0] = xz  - z * cent[0];
   mat[2][1] = yz  - z * cent[1];
   mat[2][2] = zsq - z * cent[2];

/* calculate the plane */
   if (!UTL_GEOM_SYMM_EIGENSYS ((double *)mat, 3, eval, (double *) evec, &nrot))  	return ;

   l = evec[0][0];
   m = evec[1][0];
   n = evec[2][0];
   d = (l * cent[0] + m * cent[1] + n * cent[2]);

	if (simple) {}
/* perform reflection for the requested atom IDs within the input coordinate sets */
   bool *AXclude = new bool[mol.NumAtoms() ];
   for (int i = 0; i < (int) mol.NumAtoms(); ++i) AXclude[i] = false;
   for (int i = 0; i < (int) FrameAts.size(); ++i) AXclude[FrameAts[i]] = true;
	if (simple) {
		for ( elem = 0; elem < (int)mol.NumAtoms(); elem++ )
			if (AXclude[elem]) {
				xyz = coo + (elem * 3);
				h = l * xyz[0]  +  m * xyz[1]  +  n * xyz[2]  - d;
				xyz[0] -= 2.0 * l * h;
				xyz[1] -= 2.0 * m * h;
				xyz[2] -= 2.0 * n * h;
			}
		return;
	}
//cout << "Go";
   int ats2do[mol.NumAtoms()];
   if (!doAll) {
   	getBondDist( mol, aplane[0], ats2do, AXclude);
   }
   for ( elem = 0; elem < (int)mol.NumAtoms(); elem++ ) {
        if (!doAll && ( ats2do[elem] <= 0) ) continue;
//cout << elem << " ";
        xyz = coo + (elem * 3);
        h = l * xyz[0]  +  m * xyz[1]  +  n * xyz[2]  - d;
        xyz[0] -= 2.0 * l * h;
        xyz[1] -= 2.0 * m * h;
        xyz[2] -= 2.0 * n * h;
   }
}

int getFromAtom(double *cord, int *atomdist, int baseAtom, int toAtom, int *fromAtom, int *fromTies, int atom, int natoms)
// ******************************************************************
{
    int rightPlane[3];

        if ( atomdist[atom] == 1 ||
           baseAtom < 0 || baseAtom >= natoms ||
           toAtom < 0 || toAtom >= natoms ||
           atom < 0 || atom >= natoms
        ) return -1;

    if (fromTies[atom] > -1) {
        rightPlane[0] = baseAtom;
        rightPlane[1] = atom;
        rightPlane[2] = toAtom;
        if (VECT_ISNOT_LOCAL_RIGHT(cord, rightPlane, fromAtom[atom], 2))
            return fromTies[atom];
    }
    return fromAtom[atom];
}

void fixFromAtoms( OEGraphMol mol, int now, int *toAt, int *fromAt ) {
// *****************************************************
    int t2, oldfrom;
    while (toAt[now] >= 0) {
	if (!mol.GetAtom(OEHasAtomIdx(now))->IsInRing()) return;
	t2 = toAt[now];
	if (!mol.GetAtom(OEHasAtomIdx(t2))->IsInRing()) return;
	oldfrom = fromAt[t2];
	if (oldfrom != now) fromAt[t2] = now;
	now = t2;
    }
}
void getTieTorsions( double *Tcoo, int torDa, int torDb, int torA, int torB, int torC, double *r_t1, double *r_t2 )
// *****************************************************
{
   double t1 = OEGeom3DTorsion(
        &Tcoo[torA],&Tcoo[torB],&Tcoo[torC],&Tcoo[torDa]);
   while (t1 < 0.0) t1 += 6.2831853072;
   while (t1 > 6.2831853072) t1 -= 6.2831853072;
   double t2 = OEGeom3DTorsion(
        &Tcoo[torA],&Tcoo[torB],&Tcoo[torC],&Tcoo[torDb]);
   while (t2 < 0.0) t2 += 6.2831853072;
   while (t2 > 6.2831853072) t2 -= 6.2831853072;
   *r_t1 = t1;
   *r_t2 = t2;
}

int priorityCompare(const void *vnrec, const void *vtrec ) {
// **********************************************************
        TopPriority *nrec;
        TopPriority *trec;
        int rc;
        double drc;

        nrec = (TopPriority *) vnrec;
        trec = (TopPriority *) vtrec;
        rc = trec->branchCount - nrec->branchCount;   /* larger is higher priority */
        if ( rc ) return rc;

        drc = trec->branchWeight - nrec->branchWeight;
        if ( drc > 0.001 || drc < -0.001 ) {
                if ( drc > 0.0 ) return 1;
                return -1;
        }
        drc = trec->branchDepthWeight - nrec->branchDepthWeight;
        if ( drc > 0.001 || drc < -0.001 ) {
                if ( drc > 0.0 ) return 1;
                return -1;
        }
        return 0;
}

int ScoreBrnchs( const OEGraphMol mol, int *BDist, int aroot,
	int *BrchAts, int nBranch, double *AtWts, int *nxtBest,
	int *ties, int *symties ) {
// *********************************************************

    TopPriority TPr[nBranch];
    int covered[mol.NumAtoms()], atDepth[mol.NumAtoms()], 
	whBrch[mol.NumAtoms()], coverCt[mol.NumAtoms()];

    for (int nbr = 0 ; nbr < nBranch; nbr++) {
	TPr[nbr].atomIdx = TPr[nbr].covID = TPr[nbr].branchCount = 0;
	TPr[nbr].branchWeight = TPr[nbr].branchDepthWeight = 0.0;
//		TPr[nbr].cordDist = 0.0;
    }
    for (int iat = 0; iat < (int) mol.NumAtoms(); ++iat) 
	covered[iat] = atDepth[iat] = whBrch[iat] = coverCt[iat] = 0;
    int depth = 1;
    int covMask = 1;
    int minDis = BDist[aroot];
    for (int nbr = 0; nbr < nBranch; ++nbr) {
	covered[ BrchAts[ nbr ] ] = TPr[nbr].covID = covMask;
	atDepth[ BrchAts[ nbr ] ] = depth;
	whBrch[ BrchAts[ nbr ] ] = nbr + 1;
	TPr[nbr].atomIdx = BrchAts[ nbr ];
	covMask *= 2;
    }
// enumerate each branch
    int iBrch = nBranch;
    while (iBrch > 0) {
	iBrch = 0;
	int nxtDep = depth + 1;
	for (int iat = 0; iat < (int) mol.NumAtoms(); ++iat) {
	    if (atDepth[iat] == depth && BDist[iat] > minDis) {
                OEAtomBase *apr = mol.GetAtom(OEHasAtomIdx(iat));
                for (OEIter<OEAtomBase> nbor = apr->GetAtoms(); nbor; ++nbor)
			if (BDist[nbor->GetIdx()] > minDis) {
		    if (atDepth[nbor->GetIdx()] == 0) {
		    	covered[nbor->GetIdx()] = covered[iat] | covered[nbor->GetIdx()];
			coverCt[nbor->GetIdx()] += 1;
			atDepth[nbor->GetIdx()] = depth+1;
			whBrch[nbor->GetIdx()] = whBrch[iat];
			iBrch++;
		    }
		    else if (atDepth[nbor->GetIdx()] == nxtDep) {
			if (whBrch[nbor->GetIdx()] != whBrch[iat]) {
			    covered[nbor->GetIdx()] = 0;
			    atDepth[nbor->GetIdx()] = -1;
			}
			else {
			    covered[nbor->GetIdx()] = covered[iat] | covered[nbor->GetIdx()];
                            coverCt[nbor->GetIdx()] += 1;
                            atDepth[nbor->GetIdx()] = depth+1;
                            whBrch[nbor->GetIdx()] = whBrch[iat];
                            iBrch++;
			}
		    }
		}
	    }
	}
	depth++;
    }
// score each branch
    for (int iat = 0; iat < (int) mol.NumAtoms(); ++iat) if (covered[iat]) {
	covMask = 1;
    	for (int nbr = 0; nbr < nBranch; ++nbr) {
	    if (covMask & covered[iat]) {
		TPr[nbr].branchCount += 1;
		TPr[nbr].branchWeight += AtWts[iat];
		if (atDepth[iat] > 0) TPr[nbr].branchDepthWeight +=
			AtWts[iat] * (double) atDepth[iat]; 
	    }
	    covMask *= 2;
	}
    }
// order the brnaches and record all results
    qsort( TPr, nBranch, sizeof(TPr)/nBranch, priorityCompare );
    int best1 = TPr[0].atomIdx;
    double diff1 = -1.0;
    double diff2 = -1.0;
    if (nBranch > 1) {
	nxtBest[aroot] = TPr[1].atomIdx;
	if (TPr[0].branchCount == TPr[1].branchCount)
	    diff1 = TPr[0].branchWeight - TPr[1].branchWeight;
	    diff2 = TPr[0].branchDepthWeight - TPr[1].branchDepthWeight;
	    if (diff1 < 0.0) diff1 *= -1.0;
	    if (diff2 < 0.0) diff2 *= -1.0;
	    if (diff1 < 0.001 && diff2 < 0.001) {
		ties[aroot] = TPr[1].atomIdx;
		symties[aroot] = TPr[1].atomIdx;
	    }
		
    }
    if (dbg2) {
	TClog << aroot << endl;
	for (int nbr = 0; nbr < nBranch; ++nbr) 
	TClog << " " << nbr << " " << TPr[nbr].atomIdx << " " << TPr[nbr].covID << " " << TPr[nbr].branchCount << " " << TPr[nbr].branchWeight << " " << TPr[nbr].branchDepthWeight << endl;
    }
    return best1;
}

void getAtWts( const OEGraphMol mol, double *AtWts ) {
// ***********************************************************
    for (OEIter<OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) 
	AtWts[atom->GetIdx()] = OEGetAverageWeight(atom->GetAtomicNum()) + 
		1.008 * atom->GetImplicitHCount();
}

void getAB( const OEGraphMol mol, int *BDist, bool *isTors, bool *isTorRoot, 
	bool *isProChiral ) {
// *****************************************************************
    
   for (unsigned int i = 0; i < mol.NumBonds(); ++i) isTors[i] = false;
   for (unsigned int i = 0; i < mol.NumAtoms(); ++i) isTorRoot[i] = false;
   for (OEIter<OEBondBase> bond = mol.GetBonds(); bond; ++bond)
		if (bond->IsRotor()) {
	isTors[bond->GetIdx()] = true;
 // cout << bond->GetEndIdx() << " " << bond->GetBgnIdx() << "  ";
  	int bRoot = BDist[ bond->GetBgnIdx() ];
	if (BDist[ bond->GetEndIdx() ] < bRoot) 
		bRoot = BDist[ bond->GetEndIdx() ];
	isTorRoot[ bRoot ] = true;
   }
   for (unsigned int i = 0; i < mol.NumAtoms(); ++i) isProChiral[i] = false;
   for (OEIter< OEAtomBase > nbor = mol.GetAtoms(); nbor; ++nbor ) {
	if (dbg2) {
	    TClog << nbor->GetIdx() << " " << nbor->GetAtomicNum() << " " << nbor->GetDegree() << "  (";
   	    for (OEIter< OEAtomBase > obor = nbor->GetAtoms(); obor; ++obor ) TClog << obor->GetIdx() << " " ;
   	    TClog << ")\n";
   	}
// include prochiral atosms
	if ((nbor->GetDegree() == 4) 
            && (nbor->GetDegree() - nbor->GetExplicitDegree() < 2)
// don't mess with assigned stereocenters!
		&& !nbor->HasStereoSpecified())
	isProChiral[nbor->GetIdx()] = true;
   }
}

bool rankBranches( OEGraphMol mol, int *BDist, double *AtWts,
	int *toAtom, int *nxtBest, int *ties, int *symties, int *fromAts,
	int *whBrch, int *fromties) {
// ******************************************************

// score all potential path choices -- toAtoms and fromAtoms
   int BrchAts[10];
   int mainID = -1;
   for (int dist = 1; dist < (int) mol.NumAtoms(); ++dist ) {
	for (int iat = 0; iat < (int) mol.NumAtoms(); ++iat ) {
	    if (BDist[iat] != dist) continue;
	    int nxtDist = dist+1;
	    int nBranch = 0;
	    OEAtomBase *apr = mol.GetAtom(OEHasAtomIdx(iat));
	    for (OEIter<OEAtomBase> nbor = apr->GetAtoms(); nbor; ++nbor) 
			if (BDist[nbor->GetIdx()] == nxtDist) {
		    mainID = nbor->GetIdx();	
		    BrchAts[nBranch] = mainID;
		    ++nBranch;
		}
	    if (nBranch == 1) toAtom[iat] = mainID;
	    else if (nBranch > 0)
		toAtom[iat] = ScoreBrnchs( mol, BDist, iat, BrchAts, nBranch,
			AtWts, nxtBest, ties, symties );
	}
   }
   for (int iat = 0; iat < (int) mol.NumAtoms(); ++iat ) {
	int dist = BDist[iat];
	int nxtDist = dist - 1;
	for (int i = 0; i < (int) mol.NumAtoms(); ++i ) whBrch[i] = 0;
	OEAtomBase *apr = mol.GetAtom(OEHasAtomIdx(iat));
	int nBranch = 0;
	for (OEIter<OEAtomBase> nbor = apr->GetAtoms(); nbor; ++nbor)
        		if (BDist[nbor->GetIdx()] == nxtDist) {
	    BrchAts[nBranch] = nbor->GetIdx();
	    whBrch[nbor->GetIdx()] = whBrch[nbor->GetIdx()] | (nBranch+1)*2;
	    nBranch++;
	}
	if (nBranch == 0) {fromAts[iat] = -1; continue;}
	if (nBranch == 1) {fromAts[iat] = BrchAts[0]; continue;} 
// multiple fromAts/paths (ring). Find shortest path (ties broken by local right)
	dist --;
	nxtDist = dist - 1;
	while (dist > 0) {
	  for (int kat = 0; kat < (int) mol.NumAtoms(); ++kat ) {
	    if (BDist[kat] != dist) continue;
	    if (whBrch[kat] == 0) continue;
	    OEAtomBase *apr2 = mol.GetAtom(OEHasAtomIdx(kat));
	    for (OEIter<OEAtomBase> nbor = apr2->GetAtoms(); nbor; ++nbor)
               	if (BDist[nbor->GetIdx()] == nxtDist) 
		    whBrch[nbor->GetIdx()] = whBrch[nbor->GetIdx()] |
			whBrch[kat]; 
	  }
	  dist --;
	  nxtDist = dist - 1;
	}
	int aid1 = -1;
	int aid2 = -1;
	for (int jat = 0; jat < (int) mol.NumAtoms(); ++jat ) {
	    if (BDist[jat] != 1) continue;
	    for (int nbr = 0; nbr < nBranch; ++nbr)
		if (whBrch[jat] & (nbr*2)) {
		     if (aid1 == -1) aid1 = BrchAts[nbr-1];
			else aid2 = BrchAts[nbr-1];
		}
	}
	fromAts[iat] = aid1;
	fromties[iat] = aid2;
    }

   return true;
}

/*
int main() {
    OEGraphMol mol;
    OESmilesToMol(mol,"CC(=O)OC1C(Sc2ccccc2N(C1=O)CCN(C)C)c3ccc(cc3)OC");
  oemolostream ofs;
  if (!ofs.open("last1.mol2"))
    OEThrow.Fatal("Unable to create 'last1.mol2'");
  OEWriteMolecule(ofs, mol);
  ofs.close();
  OEDetermineRingSystems(mol,Crings);
*/

void topomerize( OEGraphMol mol, double *coord, int *AMap, bool *Posed,
	const int cAnch[2] )
// *********************************************************************
{
    int *BDist = new int[mol.NumAtoms()];
    getBondDist(mol,cAnch[0],BDist,NULL);

    double *AtWts = new double[mol.NumAtoms()];
    getAtWts(mol,AtWts);
//    for (unsigned int i = 0; i < mol.NumAtoms(); ++i) cout << BDist[i] << " " << AtWts[i] << " ";
//    cout << endl;

    bool *isTors = new bool[mol.NumBonds()];
    bool *isTorRoot = new bool[mol.NumAtoms()];
    bool *isProChir = new bool[mol.NumAtoms()];
    getAB( mol, BDist, isTors, isTorRoot, isProChir );
//    for (unsigned int i = 0; i < mol.NumAtoms(); ++i) cout << isTorRoot[i] << " " << isProChir[i] << "  ";
    int toAtom[mol.NumAtoms()],nxtBest[mol.NumAtoms()],ties[mol.NumAtoms()],
	symties[mol.NumAtoms()], fromAts[mol.NumAtoms()],
	whBrch[mol.NumAtoms()], fromties[mol.NumAtoms()];
   
    for (int iat = 0; iat < (int) mol.NumAtoms(); ++iat ) 
	toAtom[iat] = nxtBest[iat] = ties[iat] = symties[iat] = 
		fromAts[iat] = fromties[iat] = -1;
// score and record all potential decisions
    rankBranches( mol, BDist, AtWts, toAtom, nxtBest, ties, symties, 
	fromAts,  whBrch, fromties);

    if (dbg2) 
    	for (int iat = 0; iat < (int) mol.NumAtoms(); ++iat )
TClog << iat << " " << toAtom[iat] << " " << fromAts[iat] << " " << nxtBest[iat] << " " << ties[iat] << " " << symties[iat] << " " << whBrch[iat] << " " << fromties[iat] << endl;

// adjust conformations (in distance order of atoms)
    int atoms[mol.NumAtoms()];
    for (int iat = 0; iat < (int) mol.NumAtoms(); ++iat ) atoms[iat] = 0;
    for (OEIter<OEBondBase> bond=mol.GetBonds(OEIsRotor()); bond; ++bond) {
	atoms[bond->GetBgnIdx()] = 1;
	atoms[bond->GetEndIdx()] = 1;
    }
    int maxDist = 0;
    for (int iat = 0; iat < (int) mol.NumAtoms(); ++iat )
	if (atoms[iat] > 0 && BDist[iat] > maxDist) maxDist = BDist[iat];
    for (int dist = 1; dist <= maxDist; ++dist)
	for (int iat = 0; iat < (int) mol.NumAtoms(); ++iat ) {
	    if (atoms[iat] == 0) continue;
	    if (Posed[iat]) continue;
// cout << iat << endl;
	    if (AMap[iat] != -1) continue;
	    if (BDist[iat] != dist) continue;
	    if (dbg2) TClog << "Topomer " << iat << ": ";
//	    bool tRing2 = false;
	    int a1 = iat;
            OEAtomBase *apr = mol.GetAtom(OEHasAtomIdx(iat));
            for (OEIter<OEAtomBase> nbor = apr->GetAtoms(); nbor; ++nbor) {
		if (BDist[nbor->GetIdx()] != dist+1) continue;
		OEBondBase *bpt = apr->GetBond(nbor);
		if (!bpt->IsRotor()) continue;
// have a relevant bond ..
		int a2 = nbor->GetIdx();
		int a0 = getFromAtom(coord, BDist, cAnch[0], toAtom[iat],
			fromAts, fromties, iat, (int) mol.NumAtoms() ); 
		int a3 = toAtom[a2];
		if ( a0 == -1 || a1 == -1 || a2 == -1 || a3 == -1 ) continue;
		if (a0 == a1 || a0 == a2 || a0 ==a3 || a1 == a2 || a2 == a3 ) continue;
		bool chkPuck = false;
		int nRingAt = 0;
		if (Crings[a1] != Crings[a2]) {
		    if (Crings[a1]) ++nRingAt;
		    if (Crings[a2]) ++nRingAt;
		    if (Crings[a2] && Crings[a3] == Crings[a2]) {
			chkPuck = true;
//		   	if (a2 == nxtBest[a1]) tRing2 = true; 
		    }
		}
		double torsion = (nRingAt ? 
		   ( nRingAt == 1 ? 1.5707963268 : 1.0471975512) 
		   : 3.1415926536 );
		if (bpt->GetOrder() > 1) torsion = 3.1415926536;
		if (iat != cAnch[0]) {
		    if (chkPuck && ties[a2] >= 0 && nxtBest[a2] >= 0) {
			double tor1, tor2;
			getTieTorsions(coord,toAtom[a2],ties[a2],a0,a1,a2,
				&tor1,&tor2);
			bool doPuck = 
			   (tor1 > 3.1415926536 && tor2 < 3.1415926536);
			if (doPuck) {
			    nxtBest[a2] = a3;
			    toAtom[a2] = a3 = ties[a2];
			    fixFromAtoms( mol, a2, toAtom, fromAts );
			}
		    }
		    OESetPackedCoords(mol,coord);
		    OESetTorsion(mol, mol.GetAtom(OEHasAtomIdx(a0)),
                	mol.GetAtom(OEHasAtomIdx(a1)),
                	mol.GetAtom(OEHasAtomIdx(a2)),
                	mol.GetAtom(OEHasAtomIdx(a3)), torsion);
		    if (dbg2) TClog << a0 << "-" << a1 << "-" << a2 << "-" << a3 << "-> << " << 57.2*torsion << "; " ;
		    OEGetPackedCoords(mol,coord);
		}
		if (chkPuck) {
		    if (dbg2) TClog << "Ring pucker check: ";
		    int plane[3];
		    plane[0] = a0;
		    plane[1] = a1;
		    plane[2] = a2;
		    std::vector<int> frameAts;
		    frameAts.push_back(a1);
		    frameAts.push_back(a2);
		    if (VECT_ISNOT_LOCAL_RIGHT(coord,plane,a3,3)) {
			int aplane[3];
			aplane[0] = a1; aplane[1] = a2; aplane[2] = a0;
			reflectAtoms(mol, coord, 3, aplane, 
				frameAts, false, false);
			if (dbg2) TClog << " to LR; " ;
		    }
		    RingPuckerLR( mol, coord, a1, a2, a3, NULL );
		}
		OEAtomBase *stem = mol.GetAtom(OEHasAtomIdx(a1));
		OEAtomBase *root = mol.GetAtom(OEHasAtomIdx(a0));
		bool doProChiral = (a1 >= 0 && a0 >= 0 && nxtBest[a1] >= 0) &&
			stem->GetDegree() >= 4 &&
			!(stem->HasStereoSpecified()) &&
			(stem->GetDegree() - stem->GetExplicitDegree() <= 1);
		if (doProChiral) {
		    if (dbg2) TClog << "  Prochiral: ";
                    int plane[3];
                    plane[0] = a0;
                    plane[1] = a1;
                    plane[2] = toAtom[a1];
		    if (VECT_ISNOT_LOCAL_RIGHT(coord,plane,nxtBest[a1],2)) {
			std::vector<int> ringAts;
				if (dbg2) TClog << "  Inverting ";
			int Nbrs[10];
			int nNbrs = 0;
			if (Crings[a1]) {
			    for (int iat2 = 0; iat2 < (int) mol.NumAtoms(); 
					++iat2 ) 
				if (Crings[iat2] == Crings[a1]) 
				    ringAts.push_back(iat2);
			    for (OEIter<OEAtomBase> atom = stem->GetAtoms(); 
					atom; ++atom)
				if (((int)atom->GetIdx()) != a0 && 
					Crings[atom->GetIdx()] == Crings[a1]) {
				    Nbrs[nNbrs] = atom->GetIdx();
				    ++nNbrs;
				}
			    if (nNbrs >= 2) {
			    	plane[0] = a1;
				if (Nbrs[0] == fromAts[a1]) {
				    plane[1] = Nbrs[0];
				    plane[2] = Nbrs[1];
				}
				else {
				    plane[2] = Nbrs[0];
				    plane[1] = Nbrs[1];
				}
				if (VECT_ISNOT_LOCAL_RIGHT(coord,plane,
						nxtBest[a2],2)) {
				    reflectAtoms(mol, coord, 3, plane,
                                	ringAts, false, false);
				    if (dbg2) TClog << " Inverted (ring); " ;
				}
			    }
			}
			else {
			    vector<int> ats2refl;
			    for (OEIter<OEAtomBase> atom = stem->GetAtoms();
                                        atom; ++atom) 
				if (((int)atom->GetIdx()) != a0)
				  for (OEIter<OEAtomBase> atom2 = 
					OEGetSubtree(stem,root); atom2; ++atom2)
				    ats2refl.push_back(atom2->GetIdx());
			    reflectAtoms(mol, coord, 3, plane, ats2refl, false, false);
			    if (dbg2) TClog << " Inverted (chain); " ;
			}
		    }
		}
		if (dbg2) TClog << endl;
	    }
	    if (dbg2) TClog << endl;
    }
/*
    delete BDist;
    delete AtWts;
    delete isTors;
    delete isTorRoot;
    delete isProChir;
*/
}
