// Performs "Template Alignment" (see "Template CoMFA" publications involving Richard D. Cramer, this code's author)

#include "openeye.h"

#include <fstream>
#include <iostream>
#include <vector>
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

#include "plane.cpp"

const char *InterfaceData =
        "!PARAMETER -tmpl\n"
                "  !TYPE string\n"
                "  !DEFAULT tmpl.sdf\n"
                "  !BRIEF 3D template file\n"
                "!END\n"
                "!PARAMETER -in\n"
                "  !TYPE string\n"
                "  !DEFAULT in.sdf\n"
                "  !BRIEF 2D input file\n"
                "!END\n"
                "!PARAMETER -out\n"
                "  !TYPE string\n"
                "  !DEFAULT out.sdf\n"
                "  !BRIEF 3D output file\n"
                "!END\n"
                "!PARAMETER -dis\n"
                "  !TYPE double\n"
                "  !DEFAULT 0.05\n"
                "  !BRIEF max distance for copying atom coords\n"
                "!END\n"
                "!PARAMETER -log\n"
                "  !TYPE int\n"
                "  !DEFAULT 0\n"
                "  !BRIEF detail level in TC.log (0-2)\n"
                "!END\n"
                "!PARAMETER -pad\n"
                "  !TYPE bool\n"
                "  !DEFAULT false\n"
                "  !BRIEF dummy mol in output if tcfa fails\n"
                "!END\n"
                "!PARAMETER -featwt\n"
                "  !TYPE double\n"
                "  !DEFAULT 0.0\n"
                "  !BRIEF weight for special SMARTS features\n"
                "!END\n"
                "!PARAMETER -n4syb\n"
                "  !TYPE int\n"
                "  !DEFAULT 0\n"
                "  !BRIEF maximum # to model\n"
                "!END\n"
                "!PARAMETER -MCSS1\n"
                "  !TYPE int\n"
                "  !DEFAULT 0\n"
                "  !BRIEF during exact MCS, atype equiv definition (0-1)\n"
                "!END\n";

bool dbg = false;
bool dbg2 = false;
double TCatDisThreshold = (double) 0.05;
ofstream TClog;

int MCSeqvLvl = 0;
unsigned int nCrings, *Crings = NULL, *nTrings, **Trings; // ring system data
double **Tbdist;    // bond separation arrays
bool *RotBds = NULL;
double *T2CatDist = NULL;  // distance between atoms (after ring positioning)
double featwt = 0.0;
bool padMolOutput = false;

#include "top.cpp"

#define nKTYPES 4
const char *kTypes[nKTYPES];
double kScore[nKTYPES];
bool ***tmplK;


void prt(string foo) {
    std::ofstream o;

    o.open("dbg.txt",std::ios_base::app);
    o << foo;
}

void hout(const vector<int> History) {
    cout << " >History: ";
    for (int i = 0; i < (int) History.size(); ++i) cout << History[i] << " ";
    cout << endl;
}

void wmol(OEGraphMol mol) {
    std::string smi;
    OECreateCanSmiString(smi, mol);
    cout << smi << endl;

    oemolostream ofs;
    if (!ofs.open("last1.sdf"))
        OEThrow.Fatal("Unable to create 'last1.sdf'");
    OEWriteMolecule(ofs, mol);
    ofs.close();
}

void cchk(OEGraphMol molC) {
// ***********************************************

    double *coo = new double[3 * (int) molC.NumAtoms()];
    OEGetPackedCoords(molC, coo);
}

double TCScoreASim(const OEAtomBase *target, const OEAtomBase *pattern) {
// *************************************************************
// returns the degree of the dissimilarity (total or individual) 
// between MCS-matched cand and template atom(s)

    double nbad = 0.0;
    if (target->GetAtomicNum() != pattern->GetAtomicNum())
        nbad += 1.0;
    if (target->GetDegree() != pattern->GetDegree())
        nbad += 0.6;
    if (target->GetFormalCharge() != pattern->GetFormalCharge())
        nbad += 0.8;
    if (target->IsInRing() != pattern->IsInRing())
        nbad += 0.25;
    if (OEAtomGetSmallestRingSize(target) !=
        OEAtomGetSmallestRingSize(pattern))
        nbad += 0.5;
    if (target->IsAromatic() != pattern->IsAromatic())
        nbad += 0.25;
    if (target->HasStereoSpecified() &&
        (target->IsAromatic() != pattern->IsAromatic()))
        nbad += 1.0;
    return nbad;
}

void bestRingPucker(OEGraphMol molC, double *Ccoo, const double *Tcoo, const int *AMap,
                    const unsigned int trunkID, const unsigned int stemID, const unsigned int sproutID) {
// ***********************************************************
    // is this a ring entry?
    if (Crings[sproutID] != Crings[stemID]) return;
    if (Crings[stemID] == Crings[trunkID]) return;

    // is T/C dist OK already?
    double C2Tdist = 0.0;
    int nT2C = 0;
    for (int i = 0; i < (int) molC.NumAtoms(); i++)
        if (Crings[i] == Crings[stemID] && AMap[i] != -1) {
            C2Tdist += UTL_GEOM_VDIST(&Ccoo[i * 3], &Tcoo[AMap[i] * 3]);
            nT2C++;
        }
    if (sqrt(C2Tdist / ((double) nT2C)) < TCatDisThreshold) return;

    // identify atom of origin ring closure
    int rgCloseID = -1;
    for (OEIter<OEAtomBase> nbor = molC.GetAtom(OEHasAtomIdx(stemID))->GetAtoms(); nbor; ++nbor)
        if (nbor->GetIdx() != sproutID && Crings[nbor->GetIdx()] == Crings[stemID]) {
            rgCloseID = nbor->GetIdx();
            break;
        }
    if (rgCloseID == -1) return;

    // is ring origin nonplanar?
    bool NplrO = false;
    int ats[4];
    ats[0] = trunkID;
    ats[1] = stemID;
    ats[2] = sproutID;
    ats[3] = rgCloseID;
    double l, m, n, d;
    if (!geom_lsplane(Ccoo, ats, 4, &l, &m, &n, &d)) return;
    double oopD = 0.0;
    for (int i = 0; i < 4; i++) {
        oopD += Ccoo[ats[i] * 3] * l;
        oopD += Ccoo[ats[i] * 3 + 1] * m;
        oopD += Ccoo[ats[i] * 3 + 2] * n;
        oopD -= d;
        if (fabs(oopD) > TCatDisThreshold) {
            NplrO = true;
            break;
        }
    }

    // is ring puckered?
    bool nonPlanarRing;
    RingPuckerLR(molC, Ccoo, rgCloseID, stemID, sproutID, &nonPlanarRing);

    if (!NplrO && !nonPlanarRing) return;

    // initialize ring pucker checking
    double bestC2T = C2Tdist;
    double *bestCoo = new double[3 * molC.NumAtoms()];
    for (unsigned int i = 0; i < 3 * molC.NumAtoms(); i++) bestCoo[i] = Ccoo[i];
    double *wkgCoo = new double[3 * molC.NumAtoms()];
    for (unsigned int i = 0; i < 3 * molC.NumAtoms(); i++) wkgCoo[i] = Ccoo[i];

    vector<int> ringAts;
    for (int i = 0; i < (int) molC.NumAtoms(); i++) if (Crings[i] == Crings[stemID]) ringAts.push_back(i);

    vector<int> subtreeAts;
    for (OEIter<OEAtomBase> atom =
            OEGetSubtree(molC.GetAtom(OEHasAtomIdx(trunkID)), molC.GetAtom(OEHasAtomIdx(stemID)));
         atom; ++atom)
        subtreeAts.push_back(atom->GetIdx());
    // trying non-planar ring origin
    if (dbg2) TClog << "  Improved C2T fit via new ring reflections: ";
    if (NplrO) {
        reflectAtoms(molC, wkgCoo, 4, ats, subtreeAts, false, true);
        C2Tdist = 0.0;
        for (int i = 0; i < (int) molC.NumAtoms(); i++)
            if (Crings[i] == Crings[stemID] && AMap[i] != -1)
                C2Tdist += UTL_GEOM_VDIST(&wkgCoo[i * 3], &Tcoo[AMap[i] * 3]);
        if (C2Tdist < bestC2T) {
            if (dbg2) TClog << "at 1st rg atom: ";
            bestC2T = C2Tdist;
            for (unsigned int i = 0; i < 3 * molC.NumAtoms(); i++) bestCoo[i] = wkgCoo[i];
        }
    }
// trying non-planar ring
    if (nonPlanarRing) {
        for (unsigned int i = 0; i < 3 * molC.NumAtoms(); i++) wkgCoo[i] = Ccoo[i];
        ats[0] = rgCloseID;
        reflectAtoms(molC, wkgCoo, 3, ats, subtreeAts, false, true);
        C2Tdist = 0.0;
        for (int i = 0; i < (int) molC.NumAtoms(); i++)
            if (Crings[i] == Crings[stemID] && AMap[i] != -1)
                C2Tdist += UTL_GEOM_VDIST(&wkgCoo[i * 3], &Tcoo[AMap[i] * 3]);
        if (C2Tdist < bestC2T) {
            if (dbg2) TClog << "at 1st rg bond: ";
            bestC2T = C2Tdist;
            for (unsigned int i = 0; i < 3 * molC.NumAtoms(); i++) bestCoo[i] = wkgCoo[i];
        }
// and both
        if (NplrO) {
            ats[0] = trunkID;
            reflectAtoms(molC, wkgCoo, 3, ats, subtreeAts, false, true);
            C2Tdist = 0.0;
            for (int i = 0; i < (int) molC.NumAtoms(); i++)
                if (Crings[i] == Crings[stemID] && AMap[i] != -1)
                    C2Tdist += UTL_GEOM_VDIST(&wkgCoo[i * 3], &Tcoo[AMap[i] * 3]);
            if (C2Tdist < bestC2T) {
                if (dbg2) TClog << "at both 1st atom AND 1st bond: ";
                bestC2T = C2Tdist;
                for (unsigned int i = 0; i < 3 * molC.NumAtoms(); i++) bestCoo[i] = wkgCoo[i];
            }
        }
        if (dbg2) TClog << '\n';
    }

// whatever reflection is best, return it
    for (unsigned int i = 0; i < (int) 3 * molC.NumAtoms(); i++) Ccoo[i] = bestCoo[i];

}

void TCsetTorsion(OEGraphMol molC, double *Ccoo, const int *AMap,
                  std::vector<int> History, const unsigned int stemID,
                  const unsigned int sproutID, bool usingT, const OEGraphMol molT,
                  const double *Tcoo) {
// ***************************************************************
// applies a torsion from template to candidate

    if (((unsigned int) History.size()) < 3) return;
    unsigned int trunkID = (unsigned int) History[(unsigned int) History.size() - 2];
    if (trunkID == sproutID || trunkID == stemID) return;
    if (AMap[trunkID] == -1) return;

    OEAtomBase *Astem = molC.GetAtom(OEHasAtomIdx(stemID));
    OEAtomBase *Atrunk = molC.GetAtom(OEHasAtomIdx(trunkID));
    OEBondBase *Broto = Astem->GetBond(Atrunk);
    if (!Broto->IsRotor()) return;

    unsigned int rootID = (unsigned int) History[(unsigned int) History.size() - 3];
    if (AMap[rootID] == -1) return;
    unsigned int torA = (unsigned int) AMap[rootID];
    unsigned int torB = (unsigned int) AMap[trunkID];
    unsigned int torC = (unsigned int) AMap[stemID];
    unsigned int torD = (unsigned int) AMap[sproutID];
    double dihed = OEGeom3DTorsion(&Tcoo[3 * torA], &Tcoo[3 * torB], &Tcoo[3 * torC],
                                   &Tcoo[3 * torD]);
    if (dbg2)
        TClog << "Tmpl Dihed 1 <=" << torA << "-" << torB << "-" << torC << "-" << torD << " : " <<
              57.2957795 * dihed << endl;
    if (usingT) {
        int nDihed = 1;
// include any other attachments to torC in dihed as average
        if (molT.GetAtom(OEHasAtomIdx(torC))->GetExplicitDegree() > 2)
            for (OEIter<OEAtomBase> atom =
                    molT.GetAtom(OEHasAtomIdx(torC))->GetAtoms(); atom; ++atom) {
                unsigned int nxtD = atom->GetIdx();
                if (nxtD == torD || nxtD == torB) continue;
                double nxtDihed = OEGeom3DTorsion(
                        &Tcoo[3 * torA], &Tcoo[3 * torB], &Tcoo[3 * torC], &Tcoo[3 * nxtD]);
                while (nxtDihed < 0.0) nxtDihed += 3.1415926536;
                if (atom->GetDegree() == 3) {
                    dihed += nxtDihed - 3.1415926536;  // ~180 degrees
                    nDihed++;
                } else if (atom->GetDegree() == 4) {    // 15 degree tolerance
                    if (fabs(nxtDihed - 2.0943951024) < 0.2617993878) {
                        dihed += nxtDihed - 2.0943951024;  // ~120 dgrees
                        nDihed++;
                    }
                    if (fabs(nxtDihed - 4.1887902048) < 0.2617993878) {
                        dihed += nxtDihed - 4.1887902048;  // ~240 degrees
                        nDihed++;
                    }
                    if (dbg2)
                        TClog << "  Dihed  <=" << torA << "-" << torB << "-" << torC << "-" << torD << " : " <<
                              57.2957795 * dihed << " (" << nDihed + 1 << ")" << endl;
                }
            }
        dihed /= (double) nDihed;
    }
    if (dbg2)
        TClog << "Cand Dihed " << rootID << "-" << trunkID << "-" << stemID << "-" << sproutID << " -> " <<
              57.2957795 * dihed << endl;

// coords back&forth so as to also rotate all distal-to-sprout atoms ..
    OESetPackedCoords(molC, Ccoo);
    OESetTorsion(molC, molC.GetAtom(OEHasAtomIdx(rootID)),
                 molC.GetAtom(OEHasAtomIdx(trunkID)),
                 molC.GetAtom(OEHasAtomIdx(stemID)),
                 molC.GetAtom(OEHasAtomIdx(sproutID)), dihed);
    OEGetPackedCoords(molC, Ccoo);

    bestRingPucker(molC, Ccoo, Tcoo, AMap, trunkID, stemID, sproutID);
}

void doCandValGeom(const unsigned int sprout, const unsigned int stem,
                   OEGraphMol molC, double *coo,
                   const double *Tcoo, const int *AMap, std::vector<int> History,
                   bool *Posed, OEGraphMol molT, double *origCoo) {
// *******************************************************
#define thCOO 3.0
#define thBLG 0.2
#define thANG .26

    if (AMap[sprout] == -1) return;
    if (dbg) TClog << "Using candidate geometry for " << stem << '-' << sprout << " bond\n";

    OESetPackedCoords(molC, coo);
    OEAtomBase *sproutA = molC.GetAtom(OEHasAtomIdx(sprout));
    OEAtomBase *stemA = molC.GetAtom(OEHasAtomIdx(stem));
    Posed[sprout] = true;
    if (Crings[sprout] > 0)
        for (unsigned int i = 0; i < molC.NumAtoms(); i++)
            if (Crings[i] == Crings[sprout])
                Posed[i] = true;

    int trunk = (int) History[History.size() - 2];
    if (AMap[trunk] == -1) return;

    double Cbdlen = UTL_GEOM_VDIST(&origCoo[sprout * 3], &origCoo[stem * 3]);
    double Tbdlen = UTL_GEOM_VDIST(&Tcoo[AMap[sprout] * 3], &Tcoo[AMap[stem] * 3]);
// TClog << trunk << "-" << stem << "-" << sprout << " " << AMap[trunk] << "-" << AMap[stem] << "-" << AMap[sprout] << endl;
    double Cang = OEGeom3DAngle(&origCoo[trunk * 3],
                                &origCoo[stem * 3], &origCoo[sprout * 3]);
    double Tang = OEGeom3DAngle(&Tcoo[AMap[trunk] * 3],
                                &Tcoo[AMap[stem] * 3], &Tcoo[AMap[sprout] * 3]);
    if (UTL_GEOM_VDIST(&coo[sprout * 3], &Tcoo[AMap[sprout] * 3]) > TCatDisThreshold ||
        fabs(Cbdlen - Tbdlen) > thBLG ||
        fabs(Cang - Tang) > thANG) {
        OESetPackedCoords(molC, coo);
        OESetDistance(molC, stemA, sproutA, Cbdlen);
        OEAtomBase *trunkA = molC.GetAtom(OEHasAtomIdx((unsigned int) trunk));
        OESetAngle(molC, trunkA, sproutA, stemA, Cang);
        OEGetPackedCoords(molC, coo);
        if (dbg)
            TClog << "  " << Tbdlen << " bleng=>" << Cbdlen << "; " << 57.3 * Tang << " valAng=>" <<
                  57.2957795 * Cang << endl;
    }
    TCsetTorsion(molC, coo, AMap, History, stem, sprout, false, molT, Tcoo);
}

void markRingAtoms(OEGraphMol molC, const unsigned int sproutID,
//	const int *AMap, bool *Visited, bool *Posed ) {
                   const int *AMap, bool *Posed) {
// **************************************************************
// all ring atoms will be posed, whether ok2copy or not
    if (dbg) TClog << "Ring Atoms: ";
    for (unsigned int ai = 0; ai < molC.NumAtoms(); ai++)
        if (Crings[ai] == Crings[sproutID]) {
            Posed[ai] = true;
// but need to unpose all unmapped side chains of this ring, for topomerization
            for (OEIter<OEAtomBase> atom =
                    molC.GetAtom(OEHasAtomIdx(ai))->GetAtoms(); atom; ++atom)
                if (AMap[atom->GetIdx()] == -1 &&
                    (Crings[atom->GetIdx()] != Crings[sproutID])) {
                    Posed[ai] = false;
                    break;
                }
            if (dbg) TClog << " " << ai;
        }
    if (dbg) TClog << endl;
}

bool try2copy(OEGraphMol molC, const OEGraphMol molT,
              const std::vector<int> History, const unsigned int stemID,
              const unsigned int sproutID, double *Ccoo, const double *Tcoo,
//	const int *AMap, bool *Posed , bool *Visited, double *origCoo ) {
              const int *AMap, bool *Posed, double *origCoo) {
// ****************************************************
    if (dbg)
        TClog << "Copying " << sproutID << " >ring " << Crings[sproutID] << " <= " << stemID << " >ring " <<
              Crings[stemID] << endl;

    if (AMap[sproutID] == -1) {
        if (dbg) TClog << "  Unmapped atom\n";
        return false;  // topomerize any unmapped atom
    }

//    Visited[sproutID] = Posed[sproutID] = true;
    Posed[sproutID] = true;

    if (UTL_GEOM_VDIST(&Ccoo[stemID * 3], &Tcoo[AMap[stemID] * 3]) > TCatDisThreshold) {
// can't copy sprout coords if stem coords wrong (from previous ring)
        if (dbg)
            TClog << "  Stem diff: " << UTL_GEOM_VDIST(&Ccoo[stemID * 3], &Tcoo[AMap[stemID] * 3]) << " " <<
                  UTL_GEOM_VDIST(&Ccoo[sproutID * 3], &Tcoo[AMap[sproutID] * 3]) << endl;

        doCandValGeom(sproutID, stemID, molC, Ccoo, Tcoo, AMap, History,
                      Posed, molT, origCoo);
        return false;
    }
    bool ok2copy = (UTL_GEOM_VDIST(&Ccoo[sproutID * 3], &Tcoo[AMap[sproutID] * 3])
                    <= TCatDisThreshold);
    if (!ok2copy) {
// try changing trunk-stem torsional rotation only 
        TCsetTorsion(molC, Ccoo, AMap, History, stemID, sproutID,
                     false, molT, Tcoo);
        ok2copy = (UTL_GEOM_VDIST(&Ccoo[sproutID * 3], &Tcoo[AMap[sproutID] * 3])
                   <= TCatDisThreshold);
        if (!ok2copy) {
            doCandValGeom(sproutID, stemID, molC, Ccoo, Tcoo, AMap, History,
                          Posed, molT, origCoo);
            return false;
        }
    }

    if (Crings[sproutID] > 0 && Crings[sproutID] == Crings[stemID]) {

// new ring entry (2nd ring atom), entire ring system is completely processed
// ring to be T=>C copied entirely if
//	smallest ring systems of same size AND
//	hybridization of this same sproutID entry point the same AND
//	(all atoms have same hybridization in T & C --
//	  assumed, conservatively, unmapped T ring atoms in the same order)
// if not entirely copied, only the torsional setting of the stem->sprout
//	bond is TC copied

//	markRingAtoms( molC, sproutID, AMap, Visited, Posed );
        markRingAtoms(molC, sproutID, AMap, Posed);
        OEAtomBase *Catm = molC.GetAtom(OEHasAtomIdx(sproutID));
        OEAtomBase *Tatm = molT.GetAtom(OEHasAtomIdx((unsigned int) AMap[sproutID]));

        ok2copy = ((OEAtomGetSmallestRingSize(Catm) ==
                    OEAtomGetSmallestRingSize(Tatm)) &&
                   (Catm->GetDegree() == Tatm->GetDegree()));
        if (ok2copy) {
// checking hybridization, ring atom by ring atom. 
// if an unmapped ring atom, ok2copy is faldse;
            for (unsigned int ai = 0; ai < molC.NumAtoms(); ai++)
                if (Crings[ai] == Crings[sproutID]) {
// need to keep checking ?
                    if (ok2copy) {
                        int tmplAt = AMap[ai];
                        if (tmplAt == -1) ok2copy = false;
                        else if (molC.GetAtom(OEHasAtomIdx(ai))->GetDegree() !=
                                 molT.GetAtom(OEHasAtomIdx((unsigned int) tmplAt))->GetDegree())
                            ok2copy = false;
                    }
                }
        }
        if (!ok2copy)
            doCandValGeom(sproutID, stemID, molC, Ccoo, Tcoo,
                          AMap, History, Posed, molT, origCoo);
        else
            for (unsigned int ai = 0; ai < molC.NumAtoms(); ai++)
                if (Crings[ai] == Crings[sproutID]) {
                    for (int xi = 0; xi < 3; ++xi)
                        Ccoo[3 * ai + xi] = Tcoo[3 * AMap[ai] + xi];
//		   TClog << ai << " <- " << AMap[ai] << endl;
                }
    }

// otherwise, still ok2copy, so just copy the atom's coords
    else {
        for (int xi = 0; xi < 3; ++xi)
            Ccoo[3 * sproutID + xi] = Tcoo[3 * AMap[sproutID] + xi];
//	TClog << sproutID << " <- " << AMap[sproutID] << " " << Ccoo[60] << " " << Tcoo[57] << endl;
    }

    return ok2copy;
}

void matchTorPS(const unsigned int sprout, const unsigned int stem,
                OEGraphMol molC, double *coo,
                const double *Tcoo, const int *AMap, std::vector<int> History,
//	bool *Posed, bool *Visited, OEGraphMol molT  ) {
                bool *Posed, OEGraphMol molT) {
// *******************************************************
// simply inherit exixting candidate valence geometry, 
// but if newly defined torsion is a rotor, match it to the template

    if (Crings[sprout] > 0 && Crings[sprout] == Crings[stem])
//	markRingAtoms( molC, sprout, AMap, Visited, Posed );
        markRingAtoms(molC, sprout, AMap, Posed);

    TCsetTorsion(molC, coo, AMap, History, stem, sprout, false, molT, Tcoo);
}

int AlignNextAtoms(OEGraphMol molC, const OEGraphMol molT, const OEAtomBase *stem,
                   const OEAtomBase *trunk, bool *Visited, bool *Posed, const int *AMap,
                   double *Ccoo, const double *Tcoo, std::vector<int> History,
                   const bool OK2copy, double *origCoo) {
// **********************************************************
//    hout(History);
    Visited[stem->GetIdx()] = true;
    bool copyOK = OK2copy;

// determining the processing order of the neighbors
// ... the criteria
    std::vector<OEAtomBase *> nbAts;
    std::vector<bool> BdInRing;
    std::vector<bool> AtMapped;
    for (OEIter<OEBondBase> bond = stem->GetBonds(); bond; ++bond) {
        OEAtomBase *nbor = bond->GetNbr(stem);
        if (nbor->GetIdx() == trunk->GetIdx())
/*        for (std::vector<int>::iterator hID = History.begin(); hID != History.end(); hID++)
            if (nbor->GetIdx() == *hID)
*/
            continue;
        nbAts.push_back(nbor);
        BdInRing.push_back(bond->IsInRing());
        AtMapped.push_back(AMap[nbor->GetIdx()] >= 0);
    }
// ... scoring each neighbor
    int nNbrs = (int) nbAts.size();
    int priority[nNbrs];
    for (int i = 0; i < nNbrs; ++i) {
        priority[i] = 0;
        if (AtMapped[i]) priority[i] += 2;
        if (BdInRing[i]) priority[i] += 1;
    }
// ... processing order is by descending score
    int NbrOrd[nNbrs];
    int ino = 0;
    for (int score = 3; score > -1; --score)
        for (int i = 0; i < nNbrs; ++i)
            if (priority[i] == score) {
                NbrOrd[ino] = i;
                ino++;
            }
    if (dbg2) {
        TClog << "From " << stem->GetIdx() << ":  ";
        for (int i = 0; i < nNbrs; ++i) {
            int nxAt = NbrOrd[i];
            TClog << "Next AtID " << nbAts[nxAt]->GetIdx() << " (Priority=" << priority[nxAt] << ",";
            if (AtMapped[nxAt]) TClog << "Mapped,";
            if (BdInRing[nxAt]) TClog << "Cyclic,";
            if (OK2copy) TClog << "GeomOK,";
        }
        TClog << ")" << endl;
/*
        cout << "From " << stem->GetIdx() << ":  ";
        for (int i = 0; i < nNbrs; ++i) {
            int nxAt = NbrOrd[i];
            cout << nbAts[nxAt]->GetIdx() << " " << priority[nxAt] <<
            "<=(" << AtMapped[nxAt] << "," << BdInRing[nxAt] << "," << OK2copy << ")  ";
        }
        cout << endl;
*/
    }

    for (int nxAt = 0; nxAt < nNbrs; ++nxAt) {
        OEChem::OEAtomBase *nbor = nbAts[NbrOrd[nxAt]];
// no back-tracking
        if (nbor->GetIdx() == trunk->GetIdx()) {
            if (dbg) TClog << nbor->GetIdx() << " matches Trunk " << trunk->GetIdx() << '\n';
            continue;
        }
/*
        bool seen = false;
        for (std::vector<int>::iterator hID = History.begin(); hID != History.end(); hID++)
            if (nbor->GetIdx() == (unsigned int) *hID) {
                seen = true;
                break;
            }
        if (seen) continue;
*/
        // ring closure ?
        if (Visited[nbor->GetIdx()]) {
            if (dbg) TClog << nbor->GetIdx() << " already visited\n";
            continue;
        }
        if (!Posed[nbor->GetIdx()]) {
            if (copyOK)
                copyOK = try2copy(molC, molT, History,
                                  stem->GetIdx(), nbor->GetIdx(),
//			Ccoo, Tcoo, AMap, Posed, Visited, origCoo);
                                  Ccoo, Tcoo, AMap, Posed, origCoo);
            if (!copyOK)
                if (AMap[nbor->GetIdx()] > -1) {
//                OEBondBase *Broto = stem->GetBond(trunk);
//                if (Broto->IsRotor())
                    matchTorPS(nbor->GetIdx(), stem->GetIdx(), molC, Ccoo, Tcoo,
                               AMap, History, Posed, molT);
//                else doCandValGeom(nbor->GetIdx(), stem->GetIdx(), molC, Ccoo,
//                                        Tcoo, AMap, History, Posed,  molT, origCoo);
                }
        } else if (dbg2) TClog << " already posed\n";
// depth first recursion ..
        History.push_back(nbor->GetIdx());
        if (dbg2) {
            TClog << " History: ";
            for (int i = 0; i < (int) History.size(); ++i) TClog << History[i] << " ";
            TClog << endl;
        }
        int nxtatm = AlignNextAtoms(molC, molT, nbor, stem, Visited, Posed,
                                    AMap, Ccoo, Tcoo, History, copyOK, origCoo);
        History.pop_back();
        if (nxtatm <= 0) return nxtatm;
    }
    return 1;
}

bool SuperPositionMol(OEMolBase &molC, double *misfitDist,
                      const OEAtomBase *Anchor1, const OEAtomBase *Anchor2,
                      const int *AMap, double *Ccoords, double *Tcoords, bool *Posed) {
// *****************************************************************
// Identifies a third atom attached to one of the anchoring atoms,
// necessary to establish the third dimension, and if necessary the
// chirality of the "root atom". 
// Performs the rigid-body superposition of the candidate onto the template, 
// including inversion of the root atm, if needed.

    double refc[9], fitc[9], rot[9], trans[3];
    int Tids[3], Cids[3];

    Cids[0] = Anchor1->GetIdx();
    Cids[1] = Anchor2->GetIdx();

// need a 3rd atom, mapped, preferably attached to Anchor1, else Anchor2
    bool fd3 = false;
//    bool bFlip = false;
    for (OEIter<OEAtomBase> nbor = Anchor1->GetAtoms(); nbor; ++nbor)
        if ((int) nbor->GetIdx() != Cids[1] &&
            AMap[nbor->GetIdx()] != -1) {
            Cids[2] = nbor->GetIdx();
            fd3 = true;
            break;
        }
    if (!fd3)
        for (OEIter<OEAtomBase> nbor = Anchor2->GetAtoms(); nbor; ++nbor)
            if (nbor->GetIdx() != Anchor1->GetIdx() &&
                AMap[nbor->GetIdx()] != -1) {
                Cids[2] = nbor->GetIdx();
                fd3 = true;
//                bFlip = true;
                break;
            }
    if (!fd3) return false;
    Posed[Cids[2]] = true;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            fitc[i * 3 + j] = Ccoords[Cids[i] * 3 + j];

    for (int i = 0; i < 3; i++) Tids[i] = AMap[Cids[i]];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            refc[i * 3 + j] = Tcoords[Tids[i] * 3 + j];

    if (dbg2) {
        TClog << "Superposing:\n";
        for (int i = 0; i < 3; i++) TClog << "  " << Cids[i] << " onto " << Tids[i] << endl;
    }
    double RMSD = OERMSD(refc, fitc, 3, true, rot, trans);
    if (dbg) TClog << "RMSD=" << RMSD;
    if (RMSD > TCatDisThreshold) {
        if (dbg2) TClog << " Superposition failed" << endl;
        return false;
    }
    if (dbg2) TClog << endl;

    OERotate(molC, rot);
    OETranslate(molC, trans);
    OEGetPackedCoords(molC, Ccoords);

    // identify the bonded sequence within Cids
    int stemID = -1;
    int trunkID, sproutID;
    for (int i = 1; i < 3; i++)
        if (!molC.GetAtom(OEHasAtomIdx((unsigned int) Cids[i]))->GetBond(
                molC.GetAtom(OEHasAtomIdx((unsigned int) Cids[0])))) {
            stemID = i == 1 ? Cids[2] : Cids[1];
            trunkID = Cids[0];
            sproutID = i == 1 ? Cids[1] : Cids[2];
            break;
        }
    if (stemID == -1) {
        stemID = Cids[0];
        sproutID = Cids[1];
        trunkID = Cids[2];
    }

// establish best "4th" coord (C=>T of "stemiD"-attached atms) by reflection through superposition plane
    vector<int> atAttach;
    for (OEIter<OEAtomBase> nbor = molC.GetAtom(OEHasAtomIdx((unsigned int) stemID))->GetAtoms(); nbor; ++nbor)
        if (nbor->GetIdx() != (unsigned int) sproutID && nbor->GetIdx() != (unsigned int) trunkID && AMap[nbor->GetIdx()] != -1)
            atAttach.push_back(nbor->GetIdx());
    if (atAttach.size() > 0) {
        int aplane[3];
        aplane[0] = trunkID;
        aplane[1] = stemID;
        aplane[2] = sproutID;

        double dist1 = 0.0;
        for (int i = 0; i < (int) atAttach.size(); i++) {
            int j = atAttach[i];
            dist1 += UTL_GEOM_VDIST(&Ccoords[j * 3], &Tcoords[AMap[j] * 3]);
        }
        vector<int> dummy;
        reflectAtoms(molC, Ccoords, 3, aplane, dummy, true, false);
        double dist2 = 0.0;
        for (int i = 0; i < (int) atAttach.size(); i++) {
            int j = atAttach[i];
            dist2 += UTL_GEOM_VDIST(&Ccoords[j * 3], &Tcoords[AMap[j] * 3]);
        }
        if (dist2 < dist1) {
            OESetPackedCoords(molC, Ccoords);
            if (dbg2) TClog << " (Reflection thru superposition plane)";
        } else OEGetPackedCoords(molC, Ccoords);
    }


// if either bond is in a ring system, all the other atoms in that ring system 
// now are also positioned
    int inRing = 0;
    unsigned int ringIDs[3] = {0, 0, 0};
    for (unsigned int cid = 0; cid < 3; ++cid)
        if (Crings[Cids[cid]] > 0) {
            ++inRing;
            ringIDs[cid] = Crings[Cids[cid]];
        }

    bool anchCyclic = false;
    if (inRing >= 2) {
// note that the Cids sequence may be either 3-1-2 or 1-2-3
        unsigned int ringID = 0;
        anchCyclic = true;
        if (ringIDs[0] == ringIDs[1]) {
            ringID = ringIDs[0];
            anchCyclic = true;
        } else if (ringIDs[0] == ringIDs[2]) ringID = ringIDs[0];
        else if (ringIDs[1] == ringIDs[2]) ringID = ringIDs[1];
        if (ringID != 0)
            for (unsigned int i = 0; i < molC.NumAtoms(); ++i)
                if (Crings[i] == ringID) {
                    Posed[i] = true;
                    for (OEIter<OEAtomBase> atom =
                            molC.GetAtom(OEHasAtomIdx(i))->GetAtoms(); atom; ++atom)
                        if (AMap[atom->GetIdx()] == -1 &&
                            (Crings[atom->GetIdx()] != ringIDs[0])) {
                            Posed[i] = false;
                            break;
                        }
                }
    }

//    if (RotBds) {delete RotBds; RotBds = NULL;}
    RotBds = new bool[molC.NumBonds()];
    for (unsigned int i = 0; i < molC.NumBonds(); ++i) RotBds[i] = false;
    for (OEIter<OEBondBase> ai = molC.GetBonds(); ai; ++ai)
        if (ai->IsRotor()) RotBds[ai->GetIdx()] = true;

    if (anchCyclic) {

        bestRingPucker(molC, Ccoords, Tcoords, AMap, (const unsigned int) trunkID, (const unsigned int) stemID,
                       (const unsigned int) sproutID);

        int tmp = trunkID;
        trunkID = sproutID;
        sproutID = tmp;
        bestRingPucker(molC, Ccoords, Tcoords, AMap, (const unsigned int) trunkID, (const unsigned int) stemID,
                       (const unsigned int) sproutID);
    }


    // computing distances between permanently positioned atoms
//    if (T2CatDist) {delete T2CatDist; T2CatDist = NULL;}
    T2CatDist = new double[molC.NumAtoms()];
    unsigned int nDist = 0;
    double totDist = 0.0;
    for (unsigned int i = 0; i < 3; ++i) {
        T2CatDist[Cids[i]] = UTL_GEOM_VDIST(&Ccoords[Cids[i] * 3],
                                            &Tcoords[Tids[i] * 3]);
        totDist += T2CatDist[Cids[i]];
        nDist += 1;
    }
    *misfitDist = totDist / (double) nDist;

    return true;
}

bool AlignMol(OEMolBase &molC, const OEGraphMol molT,
              const int cAnch[2], int *AMap) {
// ******************************************************************
// directs the complete alignment of a candidate, for a particular anchoring

    bool Visited[molC.NumAtoms()];
    bool Posed[molC.NumAtoms()];
    for (unsigned int i = 0; i < molC.NumAtoms(); i++)
        Visited[i] = Posed[i] = false;
    for (unsigned int i = 0; i < 2; i++)
        Visited[cAnch[i]] = Posed[cAnch[i]] = true;

    if (dbg) TClog << "Aligning with Anchor Atoms: " << cAnch[0] << " " << cAnch[1] << endl;
    OEAtomBase *Anchor1 = molC.GetAtom(OEHasAtomIdx((unsigned int) cAnch[0]));
    OEAtomBase *Anchor2 = molC.GetAtom(OEHasAtomIdx((unsigned int) cAnch[1]));

    double Tcoords[3 * ((int) molT.NumAtoms())];
    OEGetPackedCoords(molT, Tcoords);

    double Ccoords[3 * ((int) molC.NumAtoms())];
    OEGetPackedCoords(molC, Ccoords);

    double misfitDist = 0.0;
    if (!SuperPositionMol(molC, &misfitDist, Anchor1, Anchor2, AMap,
                          Ccoords, Tcoords, Posed)) {
        return false;
    }

    std::vector<int> History;

// provide a root so that an initial torsion can be set

    for (OEIter<OEAtomBase>
                 nbor = molC.GetAtom(OEHasAtomIdx((unsigned int) cAnch[0]))->GetAtoms(); nbor; ++nbor)
        if ((int) nbor->GetIdx() != cAnch[1] && AMap[nbor->GetIdx()] >= 0)
            History.push_back(nbor->GetIdx());

    History.push_back(cAnch[0]);
    History.push_back(cAnch[1]);

    /*
    for (int i = 0; i < 3; i++) History.push_back(embedIDs[i]);
    Anchor1 = molC.GetAtom(OEHasAtomIdx(embedIDs[1]));
    Anchor2 = molC.GetAtom(OEHasAtomIdx(embedIDs[2]));
*/
    bool OK2copy = (misfitDist < TCatDisThreshold);
    if (dbg) TClog << misfitDist << "\nDirection 1\n";

    double origCoo[3 * ((int) molC.NumAtoms())];
    OEGetPackedCoords(molC, origCoo);
    if (AlignNextAtoms(molC, molT, Anchor2, Anchor1, Visited, Posed, AMap, Ccoords, Tcoords,
                       History, OK2copy, origCoo) < 0)
        return false;

    History.clear();
    for (OEIter<OEAtomBase>
                 nbor = molC.GetAtom(OEHasAtomIdx((unsigned int) cAnch[1]))->GetAtoms(); nbor; ++nbor)
        if ((int) nbor->GetIdx() != cAnch[0] && AMap[nbor->GetIdx()] >= 0) {
            History.push_back(nbor->GetIdx());
//            break;
        }

    History.push_back(cAnch[1]);
    History.push_back(cAnch[0]);
    if (dbg) TClog << "\nDirection 2\n";
    if (AlignNextAtoms(molC, molT, Anchor1, Anchor2, Visited, Posed, AMap, Ccoords, Tcoords,
                       History, OK2copy, origCoo) < 0)
        return false;

    topomerize(molC, Ccoords, AMap, Posed, cAnch);

    OESetPackedCoords(molC, Ccoords);

    return true;
}

class LabelMCSid : public OEDisplayAtomPropBase {
public:
    LabelMCSid() {}

    string operator()(const OEChem::OEAtomBase &atom) const {
        return string(to_string((int) (1000.0 *
                                       (atom.GetPartialCharge() + 0.0001))));
    }

    OEDisplayAtomPropBase *CreateCopy() const {
        return new LabelMCSid(*this);
    }
};

void TCShowMolIdx(OEGraphMol mol, string title, bool translate) {
// *************************************************
// generates 2D image of molecule with atom IDs for MCS output

    OEPrepareDepiction(mol);

    OE2DMolDisplayOptions opts(500, 300, OEScale::AutoScale);

    if (translate) {
        LabelMCSid atomlabel;
        opts.SetAtomPropertyFunctor(atomlabel);
    } else opts.SetAtomPropertyFunctor(OEDisplayAtomIdx());

    OE2DMolDisplay disp(mol, opts);

    OERenderMolecule(title + ".png", disp);
}

void MatchOut(const OEMatchBase *miG, int cid, int tid) {
// ******************************************************
    OEIter<OEMatchPair<OEAtomBase> > apr;
    TClog << "\n CANDIDATE: " << cid << " to TEMPLATE: " << tid + 1 << "\n  tmpl atoms: ";
    for (apr = miG->GetAtoms(); apr; ++apr)
        TClog << apr->pattern->GetIdx() << ' ';
    TClog << endl;
    TClog << "  cand atoms: ";
    for (apr = miG->GetAtoms(); apr; ++apr)
        TClog << apr->target->GetIdx() << ' ';
    TClog << endl;
    //create match subgraph
    OEGraphMol m;
    OESubsetMol(m, miG, true);
    TClog << " smiles of MCSS = " << OEMolToSmiles(m) << endl;
}

double TCScoreMCS(const OEMatchBase *matchG, const int atid) {
// *************************************************************
// returns the total or individual cand <-> template atom MCS mismatch 
    OEIter<OEMatchPair<OEAtomBase> > apr;
    double nbad = 0.0;
    for (apr = matchG->GetAtoms(); apr; ++apr)
        if (atid == -1 || atid == (int) apr->target->GetIdx())
            nbad += TCScoreASim(apr->target, apr->pattern);
    return nbad;
}

bool atomsBonded(const OEGraphMol mol, unsigned int at1, unsigned int at2) {
// ***********************************************************
    OEAtomBase *Anchor1 = mol.GetAtom(OEHasAtomIdx(at1));
    for (OEIter<OEAtomBase> nbor = Anchor1->GetAtoms(); nbor; ++nbor)
        if (nbor->GetIdx() == at2) return true;
    return false;
}

bool getAnchor(int cAnch[2], const vector<double> TCScore,
               const vector<int> CandAts, const vector<int> TmplAts,
               bool *badRoots, int *nxtRoot, const OEGraphMol Templ,
               int tmplID) {
// *************************************************************
// seeks two bonded atoms within the candidate, one the "root" which must be 
// from the MCS and is attached to as many MCS atoms as possible, and
// which then maps highest among such in the characteristicness and centrality
// of its template counterpart

    struct tmAtScores {
        bool aromatic;
        int nattach;
        int valence;
        int nNbr;
        int nMustNbr;
        double score;
        double blScore;
    };
    int nMatch = (int) CandAts.size();

    int bestMusts[2];
    int isMust[Templ.NumAtoms()];

    if (dbg2) TClog << " Anchors (ID>base+Xtra): ";
    cAnch[1] = -2;
    for (unsigned int i = 0; i < Templ.NumAtoms(); i++) isMust[i] = -1;
    for (int m = 0; m < nMatch; m++) isMust[TmplAts[m]] = m;
    tmAtScores tas[nMatch];
    for (int m = 0; m < nMatch; m++) {
        OEAtomBase *apr = Templ.GetAtom(OEHasAtomIdx((unsigned int) TmplAts[m]));
        tas[m].valence = apr->GetValence();
        tas[m].nattach = apr->GetHvyDegree();
        tas[m].aromatic = apr->IsAromatic();
        tas[m].nMustNbr = 0;
        tas[m].nNbr = 0;
        for (OEIter<OEAtomBase> nbor = apr->GetAtoms(); nbor; ++nbor) {
            tas[m].nNbr++;
            if (isMust[nbor->GetIdx()] > -1) tas[m].nMustNbr++;
        }
        tas[m].score = (double) (tas[m].valence + tas[m].nattach + tas[m].nMustNbr +
                                 (tas[m].aromatic ? +1 : 0) +
                                 (apr->IsInRing() ? +1 : 0) +
                                 (apr->IsCarbon() ? -1 : 0)) - 10.0 * TCScore[m];

        tas[m].blScore = *(Tbdist[tmplID] + TmplAts[m]);
        if (dbg2) TClog << apr->GetIdx() << ">" << tas[m].score << "+" << tas[m].blScore << " ";
    }
// seeking hghest scoring must atom with two attached must atoms as root
// if none, try one must and at least one non-must  attached as root
// if none, just a polyvalent atom
    int nMust = 2;
    double bestScore;
    while (nMust >= 0) {
        bestScore = -10000.0;
        bestMusts[0] = -1;
        for (int i = 0; i < nMatch; ++i) {
            if (badRoots[i]) continue;
            if (tas[i].score - tas[i].blScore <= bestScore) continue;
            if (nMust == 2 && tas[i].nMustNbr < 2) continue;
            if (nMust == 1 && (tas[i].nMustNbr < 1 || tas[i].nNbr < 2))
                continue;
            if (nMust == 0 && tas[i].nNbr < 2) continue;
            bestScore = tas[i].score - tas[i].blScore;
            bestMusts[0] = i;
        }
        if (bestMusts[0] != -1) break;
        nMust--;
    }
    if (bestMusts[0] == -1) {
        if (dbg) TClog << "No Root atom found\n";
        return false;
    }
    *nxtRoot = bestMusts[0];
// and the other anchor atom is the next best scoring attached atom
    bestScore = -10000.0;
    bestMusts[1] = -1;

    OEAtomBase *apr = Templ.GetAtom(OEHasAtomIdx((unsigned int) TmplAts[bestMusts[0]]));
    switch (nMust) {
        case 2:
        case 1:
            for (OEIter<OEAtomBase> nbor = apr->GetAtoms(); nbor; ++nbor) {
                int midx = isMust[(int) nbor->GetIdx()];
                if (midx == -1) continue;
                double mscore = tas[midx].score;
                if (mscore > bestScore) {
                    bestScore = mscore;
                    bestMusts[1] = midx;
                }
            }
            break;
        case 0:
            for (OEIter<OEAtomBase> nbor = apr->GetAtoms(); nbor; ++nbor) {
                double mscore = (double) (nbor->GetValence() + nbor->GetHvyDegree() +
                                          //		(nbor->IsAromatic() ? -1 : 0) +
                                          //		(nbor->IsInRing() ? -1 : 0) +
                                          (nbor->IsCarbon() ? -1 : 0));
                if (mscore > bestScore) {
                    bestScore = mscore;
                    bestMusts[1] = nbor->GetIdx(); // DIFFERENT !
                }
            }
            break;
    }
    if (bestMusts[1] == -1) return false;
/*
    if (bestMusts[1] == -1 || bestMusts[0] == -1 ) {
          *nxtRoot = -1;
          if (dbg) TClog << "Root bond not found\n";
          return false;
    }
*/
// map tmpl atm IDs to cand atm IDs
    for (int i = 0; i < 2; ++i)
        for (int m = 0; m < (int) TmplAts.size(); ++m) {


            if (TmplAts[m] == TmplAts[bestMusts[i]]) {
                cAnch[i] = CandAts[m];
                break;
            }
        }
/*
    if (cAnch[0] == -1 || cAnch[1] == -1) {
          if (dbg) TClog << "Root bond not found\n";
          return false;
    }
*/
    if (dbg2) {
        TClog << "\n Anchor Atom(s):";
        for (int i = 0; i < 2; ++i) TClog << " (" << bestMusts[i] << "<=" << cAnch[i] << ")";
        TClog << endl;
    }
    return true;
}

double chkRMS(OEGraphMol molC, OEGraphMol molT, int *AMap,
              std::vector<int> Exact) {
// ****************************************************
    double Ccoo[((int) molC.NumAtoms()) * 3];
    OEGetPackedCoords(molC, Ccoo);
    double Tcoo[((int) molT.NumAtoms()) * 3];
    OEGetPackedCoords(molT, Tcoo);

    double dis2 = 0.0;
    int nwt = 0;
    int nat;
    for (nat = 0; nat < (int) molC.NumAtoms(); ++nat)
        if (AMap[nat] >= 0) {
            if (dbg2)
                TClog << "Distance: CAtm=" << nat << " Tatm=" << AMap[nat] << ": "
                      << OEGeom3DDistance(&Ccoo[nat * 3], &Tcoo[AMap[nat] * 3]) << endl;
            dis2 += OEGeom3DDistance2(&Ccoo[nat * 3], &Tcoo[AMap[nat] * 3]);
            ++nwt;
        }
    for (nat = 0; nat < (int) Exact.size(); ++nat)
        if (AMap[Exact[nat]] >= 0) {
// weighting exact matching atom distances by total of 3x
            dis2 += 2.0 * OEGeom3DDistance2(&Ccoo[Exact[nat] * 3], &Tcoo[AMap[Exact[nat]] * 3]);
            nwt += 2;
        }
    double dis = sqrt(dis2 / ((double) nwt));
    dis /= (double) Exact.size() + (double) molC.NumAtoms();
//   return dis < 10.0*TCatDisThreshold;
    return dis;
}

bool is3d(OEGraphMol molC) {
    // ***************************************************************
    double Ccoo[((int) molC.NumAtoms()) * 3];
    OEGetPackedCoords(molC, Ccoo);
    for (int i = 0; i < (int) molC.NumAtoms(); i++)
        if (Ccoo[i * 3 + 2] != 0.0) return true;
    return false;
}

bool vgOK(OEGraphMol molC) {
// ***********************************************

    double *coo = new double[3 * (int) molC.NumAtoms()];
    OEGetPackedCoords(molC, coo);

    for (OEIter<OEBondBase> ai = molC.GetBonds(); ai; ++ai) {
        double bleng = UTL_GEOM_VDIST(&coo[ai->GetBgnIdx() * 3], &coo[ai->GetEndIdx() * 3]);
        if (bleng < 1.0 || bleng > 2.4) {
            if (dbg2)
                TClog << "Bond length " << ai->GetBgnIdx() << '-' << ai->GetEndIdx() << '='
                      << bleng << "A. Conformer rejected\n";
            return false;
        }
    }
    for (OEIter<OEAtomBase> ai = molC.GetAtoms(); ai; ++ai)
        if (ai->GetExplicitDegree() > 2) {
            vector<int> atts;
            for (OEIter<OEAtomBase> nbor = ai->GetAtoms(); nbor; ++nbor)
                atts.push_back(nbor->GetIdx());
            double VA = 109.8;
            if (ai->GetDegree() < 4)
                for (OEIter<OEBondBase> bi = ai->GetBonds(); bi; ++bi) {
                    if (bi->GetOrder() == 3) {
                        VA = 180.0;
                        break;
                    }
                    if (bi->GetOrder() == 2) {
                        VA = 120.0;
                        break;
                    }
                }
            if (ai->IsAromatic()) VA = 120.0;
            if (OEAtomIsInRingSize(*ai, 3)) VA = 60.0;
            if (OEAtomIsInRingSize(*ai, 4)) VA = 90.0;

            for (int atA = 0; atA < (int) atts.size(); atA++)
                for (int atC = atA + 1; atC < (int) atts.size(); atC++) {
                    // accept rings as omega made them
                    if (Crings[atts[atA]] != 0 && Crings[atts[atC]] != 0 && Crings[ai->GetIdx()] != 0)
                        continue;
                    double valAngle = OEGeom3DAngle(&coo[3 * atts[atA]], &coo[3 * ai->GetIdx()],
                                                    &coo[3 * atts[atC]]);
                    double VA2use = VA;
                    double tolerance = 6.0;
                    if (VA < 100.0) {
                        tolerance = 20.0;
                        if (Crings[atts[atA]] != Crings[ai->GetIdx()] ||
                            Crings[atts[atC]] != Crings[ai->GetIdx()])
                            VA2use = 120.0;
                    }
                    if (ai->GetAtomicNum() == 7) tolerance = 10.0;  // amides ??
                    if (fabs(valAngle - VA2use / 57.2957795) > tolerance * TCatDisThreshold) {
                        if (dbg)
                            TClog << "Valence Angle " << atts[atA] << '-' << ai->GetIdx() << '-' << atts[atC] << '='
                                  << valAngle / 0.0174533 << " vs " << VA2use << " deg. Conformer rejected\n";
                        return false;
                    }
                }
        }

    return true;
}

void initKtmpl(const OEGraphMol mol, unsigned int tmplID, unsigned int nTmpl) {
// *********************************************
    if (tmplID == 0) tmplK = new bool **[nTmpl];
    tmplK[tmplID] = new bool *[nKTYPES];
    OESubSearch ss;
    for (int nK = 0; nK < nKTYPES; nK++) {
        tmplK[tmplID][nK] = new bool[(int) mol.NumAtoms()];
        for (int nat = 0; nat < (int) mol.NumAtoms(); nat++) tmplK[tmplID][nK][nat] = false;
        if (!ss.Init(kTypes[nK]))
            cout << "Cannot parse " << kTypes[nK] << "as SMARTS -- crashing!\n";
        for (OEIter<OEMatchBase> match = ss.Match(mol); match; ++match) {
            OEIter<OEMatchPair<OEAtomBase> > apr;
            for (apr = match->GetAtoms(); apr; ++apr) {
                int atid = apr->target->GetIdx(); //only the first SMARTS atoms
                tmplK[tmplID][nK][atid] = true;  //only the first SMARTS atoms
                break;
            }
        }
    }
}

void initKcand(const OEGraphMol mol, bool **candK) {
// *********************************************************
    for (int nK = 0; nK < nKTYPES; nK++) {
        candK[nK] = new bool[(int) mol.NumAtoms()];
        for (int nat = 0; nat < (int) mol.NumAtoms(); nat++) candK[nK][nat] = false;
        OESubSearch ss(kTypes[nK]);
        for (OEIter<OEMatchBase> match = ss.Match(mol); match; ++match) {
            OEIter<OEMatchPair<OEAtomBase> > apr;
            for (apr = match->GetAtoms(); apr; ++apr) {
                int atid = apr->target->GetIdx(); //only the first SMARTS atoms
                candK[nK][atid] = true;  //only the first SMARTS atoms
                break;
            }
        }
    }
}

int main(int argc, char **argv)
// ****************************************************************
{

    OEInterface itf(InterfaceData, argc, argv);

    TClog.open("tcfa.log");

    oemolistream ifs;
    oemolostream ofs;
    ofs.SetFormat(OEFormat::OEB);

    OEGraphMol molT;
    OEGraphMol molC;
    OEGraphMol noH;
    OEGraphMol padH;

    TCatDisThreshold = itf.Get<double>("-dis");
    int dbglvl = itf.Get<int>("-log");
    if (dbglvl > 0) dbg = true;
    if (dbglvl > 1) dbg2 = true;

    string tmplfn = itf.Get<std::string>("-tmpl");
    std::vector<OEGraphMol> Templates;
    MCSeqvLvl = itf.Get<int>("-MCSS1");
    featwt = itf.Get<double>("-featwt");
    padMolOutput = itf.Get<bool>("-pad");
    int n4syb = itf.Get<int>("-n4syb");

    if (featwt != 0.0) {
        kTypes[0] = "[!$([#6,F,Cl,Br,I,o,s,nX3,#7v5,#15v5,#16v4,#16v6,*+1,*+2,*+3])]";  //HBA
        kTypes[1] = "[!$([#6,H0,-,-2,-3])]";  //HBD
        kTypes[2] = "[*](~[*])-[R](@[R])@[R]";  // ring attachment or junction
        kTypes[3] = "[R0]=[R0]"; // acyclic double bomd
        kScore[0] = 0.8;
        kScore[1] = 0.8;
        kScore[2] = 0.35;
        kScore[3] = 0.1;
        tmplK = NULL;
    }
    if (padMolOutput) {
        OESmilesToMol(padH, "H");
        padH.SetTitle("XXX H as Pad");
    }

// read in all the templates
    if (!ifs.open(tmplfn))
        OEThrow.Fatal("Unable to open '" + tmplfn + "'");
    while (OEReadMolecule(ifs, molT)) {
// .. removing the hydrogens
        OETheFunctionFormerlyKnownAsStripSalts(molT);
        OESubsetMol(noH, molT, OEIsHeavy(), true);
        Templates.push_back(OEGraphMol(noH));
    }
    ifs.close();
    unsigned int nTmpl = (unsigned int) Templates.size();

// save template's ring systems and SSQ of bond distances for each atom pairing (to estimate centrality)
    nTrings = new unsigned int[nTmpl];
    Trings = new unsigned int *[nTmpl];
    Tbdist = new double *[nTmpl];
    for (unsigned int i = 0; i < nTmpl; ++i) {
        Trings[i] = new unsigned int[Templates[i].NumAtoms()];
        nTrings[i] = OEDetermineRingSystems(Templates[i], Trings[i]);
        int bdTmp[Templates[i].NumAtoms()];
        Tbdist[i] = new double[Templates[i].NumAtoms()];
        for (int iat = 0; iat < (int) Templates[i].NumAtoms(); ++iat) {
            getBondDist(Templates[i], (unsigned int) iat, bdTmp, NULL);
            int bdtot = 0;
            for (unsigned int jat = 0; jat < Templates[i].NumAtoms(); ++jat)
                bdtot += bdTmp[jat] * bdTmp[jat];
            *(Tbdist[i] + iat) = (double) bdtot / (10.0 * ((double) Templates[i].NumAtoms()));
        }
        if (featwt != 0.0) initKtmpl(Templates[i], i, nTmpl);
    }

// open output file
    string outfn = itf.Get<std::string>("-out");
    if (!ofs.open(outfn))
        OEThrow.Fatal("Unable to write to '" + outfn + "'");

// open candidate file
    string candfn = itf.Get<std::string>("-in");
    if (!ifs.open(candfn))
        OEThrow.Fatal("Unable to read from '" + candfn + "'");

// set up for any user-designated template atoms that must anchor alignments
    bool *hasSmarts = new bool[nTmpl];
    for (unsigned int i = 0; i < nTmpl; ++i) hasSmarts[i] = false;
    unsigned int **mustSmarts = new unsigned int *[nTmpl];
    int tid = 0;
    for (vector<OEGraphMol>::iterator Tmpl = Templates.begin();
         Tmpl != Templates.end(); ++Tmpl) {
        molT = *Tmpl;
        hasSmarts[tid] = OEHasSDData(molT, "SmartsTC");
        if (!hasSmarts[tid]) continue;
        string smarts = OEGetSDData(molT, "SmartsTC");
        OESubSearch ss(smarts.c_str());
        if (!ss) {
            cout << "Template " << tid + 1 << " has Illegal Smarts. Quitting\n" << smarts << "\n";
            return 1;
        }
        if (!ss.SingleMatch(molT)) {
            cout << "Smarts " << smarts << " not found in template " << tid + 1 << ". Will be ignored.\n";
            hasSmarts[tid] = false;
            continue;
        }
        const OEGraphMol q = ss.GetPattern();
        unsigned int nats = q.NumAtoms();
        if (nats < 3) {
            cout << "SMARTS (" << smarts << ") in Template " << tid + 1
                 << " is too small (>2 atoms needed). Will be ignored.\n";
            hasSmarts[tid] = false;
            continue;
        }
        mustSmarts[tid] = new unsigned int[nats + 1];
        mustSmarts[tid][0] = nats;
        for (OEIter<OEMatchBase> match = ss.Match(molT, true); match;
             ++match) {
            unsigned int i = 1;
            OEIter<OEMatchPair<OEAtomBase> > apr;
            for (apr = match->GetAtoms(); apr; ++apr, ++i)
                mustSmarts[tid][i] = apr->target->GetIdx();
        }
    }
//    prt("nx cand ");
    int candID = 0;
    bool **candK = NULL;
// process candidates one by one
    while (OEReadMolecule(ifs, noH)) {
        cchk(noH);
        ++candID;
        OETheFunctionFormerlyKnownAsStripSalts(noH);
        OESubsetMol(molC, noH, OEIsHeavy(), true);
        if (molC.NumAtoms() < 3) {
            TClog << "XXX Less than three heavy atoms for " << molC.GetTitle() << "\n";
            if (padMolOutput) OEWriteMolecule(ofs, padH);
            continue;
        }
        if (!is3d(molC)) cout << "Warning: All Z coords in candidate " + to_string(candID) + " = 0.0\n";
        if (dbg2) wmol(molC);
        cchk(molC);

// candidate's ring system
//        if (Crings) delete Crings;
        Crings = new unsigned int[molC.NumAtoms()];
        nCrings = OEDetermineRingSystems(molC, Crings);

// mark any user-defined additional 'match significant atoms'
        if (featwt != 0.0) {
//            if (candK != NULL) delete candK;
            candK = new bool *[nKTYPES];
            initKcand(molC, candK);
        }

// initiate best-match identification
        int cAnch[2];
        vector<int> TmplAtIDs;
        vector<int> svdExactTmplAtIDs;
        vector<int> ExactTmplAtIDs;
        vector<int> CandAtIDs;
        vector<double> TCDiffs;
        double bestScore = 0.0;
        int bestTmpl = 0;
        int nBestID = 0;
        int nMCSats = 0;
        int bestAllMCS = 0;
        int tmplID = 0;

        if (dbg) TClog << "\n\nPROCESSING CANDIDATE " << candID << " (" << OEMolToSmiles(molC) << ")\n";

// considering every template
        for (vector<OEGraphMol>::iterator Tmpl = Templates.begin();
             Tmpl != Templates.end(); ++Tmpl) {

//            prt("A");
            molT = *Tmpl;
            int mostExactMCS = 0;

// save atom IDs for specific to generic mapping
            for (OEIter<OEAtomBase> ai = molT.GetAtoms(); ai; ++ai)
                ai->SetData<OEAtomBase *>("aptr", ai);

// initiate generic mcss
            const unsigned int atomexprG = 0;
            const unsigned int bondexprG = 0;

            OEMCSSearch mcssG(molT, atomexprG, bondexprG, OEMCSType::Approximate);
            mcssG.SetMCSFunc(OEMCSMaxAtoms());

// initiate specific mcss (default case: no user-designated template atoms)
            const unsigned int atomexprS = (MCSeqvLvl == 0) ? OEExprOpts::DefaultAtoms : OEExprOpts::DefaultAtoms |
                                                                                         OEExprOpts::EqAromatic |
                                                                                         OEExprOpts::EqCAliphaticONS |
                                                                                         OEExprOpts::EqHalogen |
                                                                                         OEExprOpts::EqONS;
            const unsigned int bondexprS = (MCSeqvLvl == 0) ? OEExprOpts::DefaultBonds : OEExprOpts::DefaultBonds |
                                                                                         OEExprOpts::EqSingleDouble;

            OEMCSSearch mcssS(OEMCSType::Approximate);
            mcssS.Init(molT, atomexprS, bondexprS);
            mcssS.SetMCSFunc(OEMCSMaxAtoms());

            OEIter<OEMatchBase> miS;
/*
            if (hasSmarts)
                miS = ss.Match(molC);
            else
*/
            miS = mcssS.Match(molC, true);
            while (true) {
//                prt("b");
                if (!miS)
                    break;
                const OEMatchBase *matchS = miS;
                if (((int) matchS->NumAtoms()) < mostExactMCS) continue;
                mostExactMCS = (int) matchS->NumAtoms();

// save ids of exact matching atoms
                svdExactTmplAtIDs.clear();
                OEIter<OEMatchPair<OEAtomBase> > apr;
                for (apr = matchS->GetAtoms(); apr; ++apr) {
                    svdExactTmplAtIDs.push_back(apr->target->GetIdx());
                }
                nBestID = (int) svdExactTmplAtIDs.size();

// set up an unconstrained MCS search that must include exact matching atoms
                mcssG.ClearConstraints();

                for (OEIter<OEMatchPair<OEAtomBase> > ai = matchS->GetAtoms(); ai; ++ai) {
                    OEMatchPair<OEAtomBase> *apair = ai;
                    mcssG.AddConstraint(OEMatchPair<OEAtomBase>(apair->pattern->GetData<OEAtomBase *>("pattern_aptr"),
                                                                apair->target));
                }

// consider each result of such an unconstrained MCS 
                for (OEIter<OEMatchBase> miG = mcssG.Match(molC); miG; ++miG) {
//                    prt("C");
                    const OEMatchBase *matchG = miG;
                    int nowAllMCS = matchG->NumAtoms();
                    double nowScore = ((double) (matchG->NumAtoms() + nBestID))
                                      - TCScoreMCS(miG, -1);

                    if (featwt != 0.0)
//  modifying default MCS score? (HB, ring features, unstauration)
                        for (apr = matchG->GetAtoms(); apr; ++apr)
                            for (int nK = 0; nK < nKTYPES; nK++)
                                if (candK[nK][apr->target->GetIdx()] || tmplK[tmplID][nK][apr->pattern->GetIdx()])
                                    if (candK[nK][apr->target->GetIdx()] != tmplK[tmplID][nK][apr->pattern->GetIdx()]) {
                                        nowScore -= 10.0 * featwt * kScore[nK];
                                        if (dbg)
                                            TClog << "feat " << nK << " not matched, cand atm " << apr->target->GetIdx()
                                                  << ' ' << candK[nK][apr->target->GetIdx()] <<
                                                  ", tmpl atm " << apr->pattern->GetIdx() << ' '
                                                  << tmplK[tmplID][nK][apr->pattern->GetIdx()] << endl;
                                    }

                    if (nowScore > bestScore) {
                        if (dbg)
                            TClog << "nowMatch " << nowAllMCS << "+" << mostExactMCS << "=" << nowScore <<
                                  " > bestMatch " << bestScore << " (" << bestAllMCS << ") by  " << nBestID << " "
                                  << endl;
                        bestScore = nowScore;
                        bestAllMCS = nowAllMCS + mostExactMCS;
                        TmplAtIDs.clear();
                        CandAtIDs.clear();
                        TCDiffs.clear();

                        for (apr = miG->GetAtoms(); apr; ++apr) {
                            CandAtIDs.push_back(apr->target->GetIdx());
                            TmplAtIDs.push_back(apr->pattern->GetIdx());
                            TCDiffs.push_back(TCScoreMCS(miG, apr->target->GetIdx()));
                        }
                        nMCSats = (int) CandAtIDs.size();
// there are later exact template mappings, so save this one
                        ExactTmplAtIDs.clear();
                        for (int i = 0; i < nBestID; ++i)
                            ExactTmplAtIDs.push_back(svdExactTmplAtIDs[i]);

                        bestTmpl = tmplID;
                        if (dbg) MatchOut(matchG, candID, bestTmpl);

//                        }
                    }
                } // of generic MCSS
                ++miS;
            } // of specific MCSS
            ++tmplID;
        } // of templates

// only the best scoring MCSS result, from any template, is processed further
        if (bestScore < 0.000001) {
            TClog << "XXX No satisfactory MCS found for candidate " << molC.GetTitle() << endl;
            if (padMolOutput) OEWriteMolecule(ofs, padH);
            continue;
        }

// MCS => AMap
        int AMap[molC.NumAtoms()];
        for (unsigned int i = 0; i < molC.NumAtoms(); i++) AMap[i] = -1;
        for (int i = 0; i < nMCSats; i++) AMap[CandAtIDs[i]] = TmplAtIDs[i];
        if (dbg) {
            TClog << " Atom Mapping: ";
            for (int i = 0; i < nMCSats; i++) TClog << " " << CandAtIDs[i] << "-" << TmplAtIDs[i];
            TClog << endl;
        }

// initiate alignment with this MCS
        bool badAnchs[(unsigned int) TmplAtIDs.size()];
        for (int i = 0; i < (int) TmplAtIDs.size(); i++)
            badAnchs[i] = hasSmarts[bestTmpl] ? true : false;
        if (hasSmarts[bestTmpl])
            for (unsigned int i = 1; i < mustSmarts[bestTmpl][0]; i++) {
                int tmplAtm = mustSmarts[bestTmpl][i];
                for (int j = 0; j < nMCSats; j++) if (TmplAtIDs[j] == tmplAtm) {
                        badAnchs[j] = false;
                        break;
                }
            }

        int nbadAnch = 0;
        int rootM = -1;
        double cooInput[3 * (int) molC.NumAtoms()];
        OEGetPackedCoords(molC, cooInput);
        double cooSv[3 * (int) molC.NumAtoms()];
        bool havC = false;
        double loRMS = 1000000.0;
//        int rootSv = -1;
        bool OKanch = true;
        if (dbg2) {
            TClog << " Exact Tmpl Atoms: ";
            for (int nat = 0; nat < (int) ExactTmplAtIDs.size(); ++nat) TClog << ExactTmplAtIDs[nat] << " ";
            TClog << endl << endl;
        }
        cchk(molC);

        while (true) {
// considering each acceptable 'anchoring' for alignment
            OKanch = true;
            OESetPackedCoords(molC, cooInput);
            if (!getAnchor(cAnch, TCDiffs, CandAtIDs, TmplAtIDs, badAnchs,
                           &rootM, Templates[bestTmpl], bestTmpl)) {
                OKanch = false;
                if (rootM == -1) break;  // no possibility left
            }
            if (!atomsBonded(molC, (unsigned int) cAnch[0], (unsigned int) cAnch[1]))
                OKanch = false;

// the BIG ONE -- perform the complete alignment for this anchoring
            if (OKanch)
                if (!AlignMol(molC, Templates[bestTmpl], cAnch, AMap))
                    OKanch = false;

// valence geometries accetable?
            if (OKanch) OKanch = vgOK(molC);
            if (OKanch) {
// save alignmnt if best (has the smallest RMS over the MCS)
                double nowRMS = chkRMS(molC, Templates[bestTmpl], AMap, ExactTmplAtIDs);
                if (nowRMS < loRMS) {
                    if (dbg) TClog << "RMS distance = < " << nowRMS << " for root atm " << rootM << endl << endl;
                    havC = true;
                    OEGetPackedCoords(molC, cooSv);
                    loRMS = nowRMS;
//                    rootSv = rootM;
                }
            }
            nbadAnch++;
            if (nbadAnch >= (int) TmplAtIDs.size()) break;
            badAnchs[rootM] = true;
        }

// save the best alignment
        if (havC) {
            OESetPackedCoords(molC, cooSv);
            OEAddExplicitHydrogens(molC);
            OEAddSDData(molC, OESDDataPair("WHL_TMPL", OENumberToString(bestTmpl + 1)));
            OEAddSDData(molC, OESDDataPair("ASIM", OENumberToString(loRMS)));
            double Cfrcn = (double) CandAtIDs.size() / (double) molC.NumAtoms();
            double Tfrcn = (double) CandAtIDs.size() / (double) molT.NumAtoms();
            OEAddSDData(molC, OESDDataPair("FMAPD", OENumberToString(
                    Cfrcn > Tfrcn ? Cfrcn : Tfrcn)));
            Cfrcn = (double) ExactTmplAtIDs.size() / (double) molC.NumAtoms();
            Tfrcn = (double) ExactTmplAtIDs.size() / (double) molT.NumAtoms();
            OEAddSDData(molC, OESDDataPair("FPOSD", OENumberToString(
                    Cfrcn > Tfrcn ? Cfrcn : Tfrcn)));
            OEAddSDData(molC, OESDDataPair("AnchorAts", OENumberToString((unsigned int) cAnch[0])));
            OEWriteMolecule(ofs, molC);
        } else {
            if (padMolOutput) OEWriteMolecule(ofs, padH);
            TClog << "XXX Anchoring or Alignment failed for Candidate " << molC.GetTitle() << endl;
        }

        if (candID % 100 == 0) cout << '\n' << candID;
// cout << "finished " << molC.GetTitle() << " Cand: " << candID << " with Tmpl " << bestTmpl+1 << '\n';
        TClog << "finished Cand: " << candID << " " << molC.GetTitle() << " with Tmpl " << bestTmpl + 1 << '\n';
//        prt('\n'+ to_string(candID) +'\n');
        if (n4syb > 0 && candID > n4syb) break;

    } // of candidates

    ofs.close();

    TClog.close();

    return 0;
}
