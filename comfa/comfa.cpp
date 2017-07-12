#include "openeye.h"

#include <fstream>
#include <iostream>
#include <vector>
#include "oeplatform.h"
#include "oesystem.h"
#include "oechem.h"
#include "oedepict.h"
#include "oespicoli.h"
#include "oegrid.h"

using namespace std;
using namespace OEPlatform;
using namespace OEChem;
using namespace OESystem;
using namespace OEDepict;
using namespace OEMath;
using namespace OESpicoli;


const char *InterfaceData =
        "!PARAMETER -in\n"
                "  !TYPE string\n"
                "  !DEFAULT in.sdf\n"
                "  !BRIEF 3D from TCFA to model\n"
                "!END\n"
                "!PARAMETER -biolbl\n"
                "  !TYPE string\n"
                "  !DEFAULT LOGBIO\n"
                "  !BRIEF tag for measured affinity\n"
                "!END\n"
                "!PARAMETER -addGH\n"
                "  !TYPE bool\n"
                "  !DEFAULT true\n"
                "  !BRIEF add Gast-Huck charges\n"
                "!END\n"
                "!PARAMETER -mkmodel\n"
                "  !TYPE bool\n"
                "  !DEFAULT true\n"
                "  !BRIEF making new model\n"
                "!END\n"
                "!PARAMETER -predict\n"
                "  !TYPE bool\n"
                "  !DEFAULT true\n"
                "  !BRIEF predicting (with new or existing) model\n"
                "!END\n"
                "!PARAMETER -pred\n"
                "  !TYPE string\n"
                "  !DEFAULT pred.sdf\n"
                "  !BRIEF 3D from TC to predict\n"
                "!END\n"
                "!PARAMETER -resids\n"
                "  !TYPE bool\n"
                "  !DEFAULT false\n"
                "  !BRIEF write xv & fit model residuals to resids.dat\n"
                "!END\n"
                "!PARAMETER -ncomp\n"
                "  !TYPE int\n"
                "  !DEFAULT -1\n"
                "  !BRIEF # PLS components (-1 for SAMPLS)\n"
                "!END\n"
                "!PARAMETER -samplsall\n"
                "  !TYPE bool\n"
                "  !DEFAULT false\n"
                "  !BRIEF include all structures in SAMPLS\n"
                "!END\n"
                "!PARAMETER -model\n"
                "  !TYPE string\n"
                "  !DEFAULT model.tcfa\n"
                "  !BRIEF model file name (new or exixting)\n"
                "!END\n";


#include "pls.cpp"
#include "sampls.cpp"
#include "ghc.cpp"

//float rbAttn = 0.85f;
//ofstream TClog ;
string logbio = "LOGBIO";
bool charge = true;

void wcsv(std::string fn, double *vals, int nrow, int ncol) {
    ofstream f;
    f.open(fn);
    f.precision(5);
    for (int r = 0; r < nrow; ++r) {
        for (int c = 0; c < ncol; ++c) {
            if (c > 0) f << ',';
            f << vals[r * ncol + c];
        }
        f << endl;
    }
    f.close();
}

void wcsvF(std::string fn, double *vals, int nrow, int ncol) {
    ofstream f;

    f.open(fn);
    f.precision(5);
    double *v = vals;
    for (int r = 0; r < nrow; ++r) {
        for (int c = 0; c < ncol; ++c) {
            if (c > 0) f << ',';
            f << *v;
            v++;
        }
        f << endl;
    }
    f.close();
}

void grid2txt(string fname, double *fld, int *npts, double *lo)
// ******************************************
{

    ofstream gcout;
    gcout.open(fname);

    int xdim = npts[0];
    int ydim = npts[1];
    int zdim = npts[2];
    char buffer[80];

    int ip = 0;
    for (int iz = 0; iz < zdim; ++iz)
        for (int iy = 0; iy < ydim; ++iy)
            for (int ix = 0; ix < xdim; ++ix) {
                double x = lo[0] + 2.0f * ((double) ix);
                double y = lo[1] + 2.0f * ((double) iy);
                double z = lo[2] + 2.0f * ((double) iz);
                sprintf(buffer, "%6.2f %6.2f %6.2f %-12.6f\n", x, y, z, fld[ip]);
                gcout << buffer;
                ip++;
            }
    gcout.close();
}

void savemodel(string modelfn, int nFpt, double intcpt, double *beta,
               int *npts, double *lo, int ncomp, double sdep, double q2,
//               double r2, double s, int incomp) {
               double r2, double s, int incomp, double *ster, double *elec) {
// ********************************************
    ofstream f;
    f.open(modelfn);
    f << "Stats: #comp= " << ncomp;
    if (incomp == -1) f << " (q2= " << q2 << " SDEP= " << sdep << " ) ";
    f << " r2= " << r2 << " s= " << s << endl;
    for (int i = 0; i < 3; i++) f << npts[i] << endl;
    for (int i = 0; i < 3; i++) f << lo[i] << endl;

    f << intcpt << endl;
    for (int i = 0; i < nFpt; i++) 
	f << beta[i] << ',' << (i < nFpt/2 ? ster[i] : elec[i-nFpt/2])  << endl;
    f.close();

    float midX = (float) lo[0] + (float) npts[0];
    float midY = (float) lo[1] + (float) npts[1];
    float midZ = (float) lo[2] + (float) npts[2];
    OEScalarGrid grid(npts[0], npts[1], npts[2], midX, midY, midZ, 2.0f);

    OESurface surf;

    oeofstream sfs("ste.oesrf");
    for (int i = 0; i < nFpt/2; i++) grid[i] = (float) ster[i];
    OEMakeSurfaceFromGrid(surf,grid,0.0f,0.5f);
    OESetSurfaceColor(surf,255.0f,255.0f,0.0f,1.0f);
    OEWriteSurface(sfs, surf, OESurfaceFileType::OESurface);

    for (int i = 0; i < nFpt/2; i++) grid[i] = 0.0f - (float) ster[i];
    OEMakeSurfaceFromGrid(surf,grid,0.0f,0.5f);
    OESetSurfaceColor(surf,0.0f,255.0f,0.0f,1.0f);
    OEWriteSurface(sfs, surf, OESurfaceFileType::OESurface);

    sfs.close();

    oeofstream efs("ele.oesrf");
    for (int i = 0; i < nFpt/2; i++) grid[i] = (float) elec[i];
    OEMakeSurfaceFromGrid(surf,grid,0.0f,0.5f);
    OESetSurfaceColor(surf,255.0f,0.0f,0.0f,1.0f);
    OEWriteSurface(efs, surf, OESurfaceFileType::OESurface);

    for (int i = 0; i < nFpt/2; i++) grid[i] = 0.0f - (float) elec[i];
    OEMakeSurfaceFromGrid(surf,grid,0.0f,0.5f);
    OESetSurfaceColor(surf,0.0f,0.0f,255.0f,1.0f);
    OEWriteSurface(efs, surf, OESurfaceFileType::OESurface);

    efs.close();
}

/*
void addRBAtm( OEGraphMol mol, int nowAtm, bool *atAdded, vector<int> &NextLyr, 
	vector<int> &NextLyrFrom, bool *visited, double *wts, double wt ) {
// ******************************************

//cout << nowAtm << " add\n";
  for (OEIter<OEAtomBase>
        nbor = mol.GetAtom(OEHasAtomIdx((unsigned int) nowAtm))->GetAtoms(); nbor; ++nbor) {
// cout << ' ' << nbor->GetIdx() ;
     if (visited[nbor->GetIdx()]) continue;
     *atAdded = true;
      if (nbor->GetIdx() >= mol.NumAtoms()) cout << "xy\n";
     wts[nbor->GetIdx()] = wt;
     if (nbor->GetDegree() == 1) continue;
     OEBondBase *cBond = nbor->GetBond(mol.GetAtom(OEHasAtomIdx((unsigned int) nowAtm)));
     if (!cBond->IsInRing() && cBond->GetOrder() == 1 && 
		    !visited[(int) nbor->GetIdx()] ) {
         bool newnxt = true;
	    if ((int) NextLyr.size() > 0)
	    for (int i = 0; i < (int) NextLyr.size(); ++i)
		    if (NextLyr[i] == (int) nbor->GetIdx()) { newnxt = false; break;}
	    if (newnxt) {
	        NextLyr.push_back(nbor->GetIdx());
	        NextLyrFrom.push_back(nowAtm);
	    }
     }
     else {
         if (nbor->GetIdx() >= mol.NumAtoms()) cout << "xy\n";
	    visited[nbor->GetIdx()] = true;
	    wts[nbor->GetIdx()] = wt;
        addRBAtm(mol,nbor->GetIdx(),atAdded,NextLyr,NextLyrFrom,visited,wts,wt);
     }
  } 
}

void getRbAtWts( double *wts, OEGraphMol mol ) {
// ********************************************
  double currWt = 1.0;
  bool visited[(int) mol.NumAtoms() ];
  for (int i = 0; i < (int) mol.NumAtoms(); ++i) visited[i] = false;
  vector<int> NowLayer;
  vector<int> NextLayer;
  vector<int> NextLayerFrom;

  int rootAt = 0;
  OEStringToNumber(OEGetSDData(mol,"AnchorAts"), rootAt);

  wts[rootAt] = currWt; 
  visited[rootAt] = true;
  NowLayer.push_back(rootAt);
  while (true) {
	bool atAdded = false;
	NextLayer.clear();
	NextLayerFrom.clear();
	for (int i = 0; i < (int) NowLayer.size(); ++i) 
	    addRBAtm(mol,NowLayer[i],&atAdded,NextLayer,NextLayerFrom,visited,wts,currWt); 
	if (!atAdded) break;
	NowLayer.clear();
	for (int i = 0; i < (int) NextLayer.size(); ++i) {
		NowLayer.push_back( NextLayer[i]);
        // if ( NextLayer[i] >= mol.NumAtoms()) cout << "xy\n";
		visited[ NextLayer[i]] = true;
		wts[ NextLayer[i]] = currWt;
	}
	currWt *= rbAttn;
    NextLayer.clear();
  }
}
 */
static double scutoff[16] = {9999., 0., 2., 4., 6., 8., 10., 12., 14., 16., 18., 20., 22., 24., 26., 30.};
static double scutValues[16] = {9999., -0.1, 1., 3., 5., 7., 9., 11., 13., 15., 17., 19., 21., 23., 25., 29.};

static double ecutoff[16] = {9999., -12.9999, -10.9999, -8.9999, -6.9999, -4.9999, -2.9999, -0.9999, 1.0001, 3.0001,
                             5.0001, 7.0001, 9.0001, 11.0001, 13.0001, 15.0001};
static double ecutValues[16] = {9999., -13.9999, -11.9999, -9.9999, -7.9999, -5.9999, -3.9999, -1.9999, 0.0001, 2.0001,
                                4.0001, 6.0001, 8.0001, 10.0001, 12.0001, 14.0001};

static double LookupBinnedValue(double value, int isElectroStatic) {
// **********************************************************
    double *cutoff;
    double *cutoffValues;
    int i;

    if (isElectroStatic) {
        cutoff = ecutoff;
        cutoffValues = ecutValues;
    }
    else {
        cutoff = scutoff;
        cutoffValues = scutValues;
    }

    for (i = 1; i < 16; i++)
        if (value <= cutoff[i])
            return cutoffValues[i];

    return 0.00;
}

bool calcField(OEGraphMol mol, double *ster, double *elec, int *nstep, double *basept,
               int molID) {
// *****************************************************
#define STERIC_MAX 30.0f
#define Q2KC 332.0f
#define MIN_SQ_DISTANCE 1.0e-4
#define ELECTRO_MAX 14.0f
#define ELECTRO_MIN -14.0f

    int tnstep = nstep[0] * nstep[1] * nstep[2];
    int nats = (int) mol.NumAtoms();
    double *rbAtWts = new double[nats];
    for ( int i = 0; i < nats; i++) rbAtWts[i] = 1.0;
// getRbAtWts( rbAtWts, mol );

    double *AtChg = new double[nats];
    for ( int i = 0; i < nats; i++) AtChg[i] = 0.0;

    if (charge) {
        if (!GHCharges(mol, AtChg)) {
            cout << "Partial charge calculations failed for structure " << molID << endl;
            return false;
        }
    }
    else for (OEIter<OEAtomBase> atom = mol.GetAtoms(); atom; ++atom)
            AtChg[atom->GetIdx()] = atom->GetPartialCharge();

    double *vdwA = new double [nats];
    double *vdwB = new double [nats];

    for (OEIter<OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        double radius = 0.0, eps = 0.0;
        int anum[12] = {1, 6, 7, 8, 9, 11, 14, 15, 16, 17, 35, 53};
        double epsval[12] =
                {0.042, 0.107, 0.095, 0.116, 0.400, 0.0, 0.042, 0.314, 0.314, 0.4, 0.434, 0.6};
        double radval[12] =
                {1.5, 1.7, 1.55, 1.52, 1.47, 1.2, 2.1, 1.8, 1.8, 1.75, 1.85, 1.98};
        bool haveAW = false;
        for (int naw = 0; naw < 12; ++naw)
            if (((int) atom->GetAtomicNum()) == anum[naw]) {
                haveAW = true;
                radius = radval[naw];
                eps = epsval[naw];
            }
        if (!haveAW) {
            radius = 1.7;
            eps = .107;
        }
        radius += 1.7;
        eps = sqrt(eps * .107);
        if (atom->GetIdx() >= mol.NumAtoms()) cout << "xy\n";

        vdwA[atom->GetIdx()] = 2.0 * eps * rbAtWts[atom->GetIdx()] *
                               pow(radius, 6.0);
        vdwB[atom->GetIdx()] = eps * rbAtWts[atom->GetIdx()] *
                               pow(radius, 12.0);
    }
    double *coo = new double[3 * (int) mol.NumAtoms()];
    OEGetPackedCoords(mol, coo);

    int nstepx = nstep[0];
    int nstepy = nstep[1];
    int nstepz = nstep[2];
    double stepx = 2.0;
    double stepy = 2.0;
    double stepz = 2.0;
    double lowx = basept[0];
    double lowy = basept[1];
    double lowz = basept[2];
    /*
    nstepx = nstepy = nstepz = 2;
    stepx = stepy = stepz = 8.0;
    lowx = lowy = lowz = -5.0;
*/double *coord = coo;
    int max_steps = (int) (4.0 / stepx);
    if (max_steps <= 0 || ((double) max_steps * stepx) < 4.0) max_steps += 1;
    for (int nat = 0; nat < (int) mol.NumAtoms(); ++nat) {
        double curr_x = *coord;
        double curr_y = *(coord + 1);
        double curr_z = *(coord + 2);
        coord += 3;
        //    ofstream ofs("at1.elec");

        int iz = (int) (fabs(curr_z - lowz + 0.5) / stepz);
        int iy = (int) (fabs(curr_y - lowy + 0.5) / stepy);
        int ix = (int) (fabs(curr_x - lowx + 0.5) / stepx);

        int curr_iz = iz - max_steps;
        int curr_iy = iy - max_steps;
        int curr_ix = ix - max_steps;
        int curr_nstepsz = iz + max_steps + 1;
        int curr_nstepsy = iy + max_steps + 1;
        int curr_nstepsx = ix + max_steps + 1;
        if (curr_iz < 0) curr_iz = 0;
        if (curr_iy < 0) curr_iy = 0;
        if (curr_ix < 0) curr_ix = 0;
        if (curr_iz >= nstepz) curr_iz = nstepz - 1;
        if (curr_iy >= nstepy) curr_iy = nstepy - 1;
        if (curr_ix >= nstepx) curr_ix = nstepx - 1;
        if (curr_nstepsz > nstepz) curr_nstepsz = nstepz;
        if (curr_nstepsy > nstepy) curr_nstepsy = nstepy;
        if (curr_nstepsx > nstepx) curr_nstepsx = nstepx;
        double maxw = STERIC_MAX * rbAtWts[nat];
        double z;
        double y;
        double x;
        double *st;
        for (iz = 0, z = lowz; iz < nstepz; iz++, z += stepz) {
            double zd = z - curr_z;
            zd = zd * zd;
            for (iy = 0, y = lowy; iy < nstepy; iy++, y += stepy) {
                double yd = y - curr_y;
                yd = yd * yd;
                st = ster + (iz * nstepy * nstepx) + (iy * nstepx) + curr_ix;
                for (ix = 0, x = lowx; ix < nstepx; ix++, x += stepx) {
                    double sum_steric = *st;
                    double xd = x - curr_x;
                    double dis2 = xd * xd + yd + zd;
                    if (iz >= curr_iz && iz < curr_nstepsz &&
                        iy >= curr_iy && iy < curr_nstepsy &&
                        ix >= curr_ix && ix < curr_nstepsx) {
                        double atm_steric;
                        if (dis2 >= MIN_SQ_DISTANCE) {
                            double dis6 = dis2 * dis2 * dis2;
                            double dis12 = dis6 * dis6;
                            atm_steric = vdwB[nat] / dis12 - vdwA[nat] / dis6;
                            if (atm_steric > maxw) atm_steric = maxw;
//                          if (scrit && st == 6963 && atm_steric > 1.0) cout << nat << endl;
                        }
                        else atm_steric = maxw;
                        sum_steric += atm_steric;
                        *st = sum_steric > STERIC_MAX ? STERIC_MAX : sum_steric;
                        st++;
                    }
                    if (dis2 == 0.0) dis2 = 0.000001;
                    double *el = elec + (iz * nstepy * nstepx) + (iy * nstepx) + ix;
                    double eval = (AtChg[nat] / dis2) * rbAtWts[nat];
                    *el += eval;
                    /*
                    ofs << nat << ' ' << ix << ' ' << iy << ' ' << iz << ' ' << eval << ' ' << AtChg[nat] << ' ';
                    ofs << x << ' ' << curr_x << ' ' << y << ' ' << curr_y << ' ' << z << ' ' << curr_z << ' ' ;
                    ofs  << dis2 << endl;
                     */
                    if (nat == (int) mol.NumAtoms() - 1) {
                        *el *= Q2KC;
                        if (*el > ELECTRO_MAX) *el = ELECTRO_MAX;
                        else if (*el < ELECTRO_MIN) *el = ELECTRO_MIN;
                    }
                }  // X
            } // Y
        } // Z
//      ofs.close();
    } //atom

    for (int np = 0; np < tnstep; np++) {
        ster[np] = LookupBinnedValue(ster[np], false);
        elec[np] = LookupBinnedValue(elec[np], true);
    }

/*    if (molID <= 2 ) {
        grid2txt(to_string(molID) + "ste.fd", ster, nstep, basept);
        grid2txt(to_string(molID) + "ele.fd", elec, nstep, basept);
    }
grid2txt("ste.oe.fd",ster,true);
grid2txt("ele.oe.fd",elec,true);
*/
    return true;
}

void wresids(const double *xvresdl, const double *resdl, int nMol, const bool *OKmol) {
// *************************************************************
    ofstream fres("resids.txt");

    for (int i = 0; i < nMol; i++) {
        fres << i + 1;
        if (OKmol[i]) fres << ',' << xvresdl[i] << ',' << resdl[i * 2];
        fres << std::endl;
    }
}

void mpredict(int nFpt, string predfn, double intcpt, double *beta,
              double *ster, double *elec, int *nstep, double *basept) {
// **********************************************************************
    int nMol = 0;
    int nOKMol = 0;
    OEGraphMol molC;
    double valBio;
    oemolistream ifs;

    double sumRes = 0.0;
    double sumSsq = 0.0;
    ofstream fResid;
    fResid.open("predx.txt");

    if (!ifs.open(predfn))
        OEThrow.Fatal("Unable to read from '" + predfn + "'");
    while (OEReadMolecule(ifs, molC)) {
        ++nMol;
        if (nMol % 100 == 0) cout << "Predicting molecule " << nMol << endl;
        memset(ster, 0, (nFpt / 2) * sizeof(double));
        memset(elec, 0, (nFpt / 2) * sizeof(double));
        if (!calcField(molC, ster, elec, nstep, basept, nMol)) {
            cout << "Field calculation failed for structure " << nMol << endl;
            fResid << nMol << endl;
            continue;
        }

        double predB = intcpt;
        int npt = nFpt / 2;
        for (int i = 0; i < nFpt; i++)
            predB += beta[i] * (i < npt ? ster[i] : elec[i - npt]);

        if (OEHasSDData(molC, logbio)) {
            OEStringToNumber(OEGetSDData(molC, logbio), valBio);
            double resid = predB - valBio;
            sumRes += resid;
            sumSsq += resid * resid;
            nOKMol++;
            fResid << nMol << ',' << predB << ',' << valBio << ',' << resid;
        }
        else fResid << nMol << ',' << predB;
        /*
        if (xs) {
            double ssq = 0.0;
            for (int i = 0; i < nFpt; i++) {
                double val = fabs(xs[nMol*nFpt+i] - mean[i]);
                ssq += val * val;
            }
            double vdiff = sqrt( ssq / ((double)nFpt));
            fResid << ',' << vdiff;
        }
         */
        fResid << endl;
    }
    fResid.close();
}

void fpredict(string modelfn, string predfn) {
// ***************************************************
    ifstream mdl;
    mdl.open(modelfn);
    if (!mdl) {
        cout << "Cannot open model file '" << modelfn << "' for reading\n";
        return;
    }
    string stats;
    getline(mdl, stats);
    double lo[3];
    int npts[3];

    for (int i = 0; i < 3; i++) mdl >> npts[i];
    for (int i = 0; i < 3; i++) mdl >> lo[i];

    int nFpts = npts[0] * npts[1] * npts[2];
    double *ster = new double[nFpts];
    double *elec = new double[nFpts];

    double intcpt;
    mdl >> intcpt;
    double *beta = new double[nFpts];
    
//    for (int i = 0; i < 2 * nFpts; i++) mdl >> beta[i];
    for (int i = 0; i < 2 * nFpts; i++) {
	mdl >> stats;
	size_t comma = stats.find(",");
        if (comma == std::string::npos) beta[i] = atof(stats.c_str());
	   else beta[i] = atof(stats.substr(0,comma).c_str());
    }

    mpredict(nFpts, predfn, intcpt, beta, ster, elec, npts, lo);

}

int main(int argc, char **argv) {
// ******************************************************


    OEInterface itf(InterfaceData, argc, argv);
    string candfn = itf.Get<std::string>("-in");
    bool mkmodel = itf.Get<bool>("-mkmodel");
    charge= itf.Get<bool>("-addGH");
    int ncomp = itf.Get<int>("-ncomp");

    bool predict = itf.Get<bool>("-predict");
    bool resids = itf.Get<bool>("-resids");
    bool sampall = itf.Get<bool>("-samplsall");
    logbio = itf.Get<std::string>("-biolbl");
    string modelfn = itf.Get<std::string>("-model");
    string predfn = itf.Get<std::string>("-pred");

    if (!mkmodel && !predict) {
        cout << "Nothing to do: Neither making a model nor predicting requested\n";
        return 0;
    }

    if (mkmodel) {

        // so we are to make a new model
        unsigned int numMol = 0;
        OEGraphMol mol;
        double lo[3];
        double hi[3];
        int npts[3];
        oemolistream ifs;
        if (!ifs.open(candfn))
            OEThrow.Fatal("Unable to read (tcfa-produced) structures from '" + candfn + "'");
        bool first = true;
        while (OEReadMolecule(ifs, mol)) {
            numMol += 1;
            double *coo = new double[3 * mol.NumAtoms()];
            OEGetPackedCoords(mol, coo);
            if (first) {
                lo[0] = hi[0] = coo[0];
                lo[1] = hi[1] = coo[0];
                lo[2] = hi[2] = coo[0];
                first = false;
            }
            for (unsigned int i = 0; i < mol.NumAtoms(); i++)
                for (int j = 0; j < 3; j++) {
                    if (coo[i * 3 + j] < lo[j]) lo[j] = coo[i * 3 + j];
                    if (coo[i * 3 + j] > hi[j]) hi[j] = coo[i * 3 + j];
                }
//            delete coo;
        }
        ifs.close();

        for (int j = 0; j < 3; j++) {
            lo[j] -= 4.0f;
            hi[j] += 4.0f;
            npts[j] = ((int) ((hi[j] - lo[j]) / 2.0)) + 1;
        }

        vector<double *> fields;
        int npt = npts[0] * npts[1] * npts[2];
        int nFpt = 2 * npt;
        double *ster = new double[npt];
        double *elec = new double[npt];
        double *sumFv = new double[nFpt];
        double *ssqFv = new double[nFpt];
        double *hiFv = new double[nFpt];
        double *loFv = new double[nFpt];
        for (int i = 0; i < nFpt; ++i) {
            sumFv[i] = ssqFv[i] = 0.0f;
            hiFv[i] = -1.0e05f;
            loFv[i] = 1.0e10f;
        }
        vector<double> bioVs;
        double valBio;

        bool *OKmol = new bool[numMol];
        for (unsigned int i = 0; i < numMol; ++i) OKmol[i] = true;

        if (!ifs.open(candfn))
            OEThrow.Fatal("Unable to read from '" + candfn + "'");
        int nMol = 0;
        int nOKMol = 0;
        OEGraphMol molC;
//        ofstream ofs("14663.txt");

        while (OEReadMolecule(ifs, molC)) {
//if (nMol >= 10) break;
            ++nMol;

            if (!OEHasSDData(molC, logbio)) {
                cout << "No " << logbio << " value for structure " << nMol << endl;
                OKmol[nMol] = false;
                continue;
            }

            memset(ster, 0, npt * sizeof(double));
            memset(elec, 0, npt * sizeof(double));

            if (!calcField(molC, ster, elec, npts, lo, nMol)) {
                cout << "Field calculation failed for structure " << nMol << endl;
                OKmol[nMol] = false;
                continue;
            }

            ++nOKMol;

            // nxt 2 lines are debugging alternative to calcField
            /*
            for (int i = 0; i < npt/20; i++) ster[randint(npt)] = 30.0;
            for (int i = 0; i < npt/20; i++) elec[randint(npt)] = (i % 2 == 0) ? -15.0 : +15.0;
    */
            OEStringToNumber(OEGetSDData(molC, logbio), valBio);
            bioVs.push_back(valBio);

            double *nowFld = new double[nFpt];
            memcpy(nowFld, ster, sizeof(double) * npt);
            memcpy(nowFld + npt, elec, sizeof(double) * npt);
            fields.push_back(nowFld);
            for (int i = 0; i < nFpt; i++) {
                sumFv[i] += nowFld[i];
                ssqFv[i] += nowFld[i] * nowFld[i];
//                if (i == 14662) ofs << nowFld[i] << endl;
                if (nowFld[i] < loFv[i]) loFv[i] = nowFld[i];
                if (nowFld[i] > hiFv[i]) hiFv[i] = nowFld[i];
            }
        }
//        ofs.close();

        wcsvF("sum.fd", sumFv, nFpt, 1);
        wcsvF("ssq.fd", ssqFv, nFpt, 1);
        cout << "Fields done\n";
// identifying columns where no field variance was experienced
        bool *OKfVal = new bool[nFpt];
        for (int i = 0; i < nFpt; ++i) OKfVal[i] = false;

        int nOKfVal = 0;
        for (int i = 0; i < nFpt; ++i)
            if (ssqFv[i] != 0.0f && hiFv[i] != loFv[i]) {
                ++nOKfVal;
                OKfVal[i] = true;
            }
// dropping missing rows (structures) ..
        double *ys = new double[nOKMol];
        int iy = 0;
        for (int i = 0; i < nMol; i++)
            if (OKmol[i]) {
                ys[iy] = (double) bioVs[iy];
                ++iy;
            }
        double *xs = new double[nOKMol * nOKfVal];
        double *xsig = new double[nOKfVal + 2];
        double *xmean = new double[nOKfVal + 2];
        double *sig = new double[nOKfVal];
        double *mean = new double[nOKfVal];
        for (int i = 0; i < nOKfVal + 2; i++) xsig[i] = xmean[i] = 0.0;
        for (int i = 0; i < nOKfVal; i++) sig[i] = mean[i] = 0.0;

        int nster = 0;
        int ifd = 0;
        int ix = 0;
        for (int i = 0; i < nMol; i++)
            if (OKmol[i]) {
                double *Vals = fields[ix];
                int jfd = 0;
// .. dropping signal-less lattice points
                for (int j = 0; j < nFpt; j++)
                    if (OKfVal[j]) {
                        xs[ifd] = Vals[j];
                        // and totaling by field type
                        xmean[nOKfVal + ((j < nFpt / 2) ? 0 : 1)] += xs[ifd];
                        xsig[nOKfVal + ((j < nFpt / 2) ? 0 : 1)] += xs[ifd] * xs[ifd];
                        mean[jfd] += xs[ifd];
                        sig[jfd] += xs[ifd] * xs[ifd];
                        if (i == 0 && j < nFpt / 2) nster++;
                        ++ifd;
                        ++jfd;
                    }
                ++ix;
            }
        // CoMFA field weights are by the field, not the lattice point
        xsig[nOKfVal] = ((xsig[nOKfVal] - ((xmean[nOKfVal] * xmean[nOKfVal]) / (double) (nster * nOKMol))) /
                         ((double) ((nster * nOKMol) - 1)));
        xsig[nOKfVal] = sqrt(xsig[nOKfVal] * ((double) nster));
        xsig[nOKfVal + 1] = (
                (xsig[nOKfVal + 1] -
                 ((xmean[nOKfVal + 1] * xmean[nOKfVal + 1]) / (double) ((nOKfVal - nster) * nOKMol))) /
                ((double) ((nOKfVal - nster) * nOKMol - 1)));
        xsig[nOKfVal + 1] = sqrt(xsig[nOKfVal + 1] * ((double) (nOKfVal - nster)));
        xmean[nOKfVal] /= (double) (nOKMol * nster);
        xmean[nOKfVal + 1] /= (double) (nOKMol * (nOKfVal - nster));
        for (int j = 0; j < nOKfVal; j++) {
            xsig[j] = j < nster ? xsig[nOKfVal] : xsig[nOKfVal + 1];
            xmean[j] = j < nster ? xmean[nOKfVal] : xmean[nOKfVal + 1];
            mean[j] /= (double) nMol;
            sig[j] = sig[j] > 0.0 ? sqrt((sig[j] - ((double) nMol) * mean[j] * mean[j]) / ((double) (nMol - 1))) : 0.0;
        }
        for (int i = 0; i < nOKMol; i++)
            for (int j = 0; j < nOKfVal; j++) {
                if (i * nOKfVal + j > nOKMol * nOKfVal) cout << "xx\n";
                xs[i * nOKfVal + j] = (xs[i * nOKfVal + j] - xmean[j]) / xsig[j];
            }
/*wcsv("xavg.csv", xmean, nOKfVal, 1);
        wcsv("x.csv", xs, nOKMol, nOKfVal);
        wcsv("y.csv", ys, nOKMol, 1);
        wcsv("xsdv.csv", xsig, nOKfVal, 1);
*/
        cout << "Field values filtered and scaled\n";
        double *ywts = new double[nOKMol];
        int ierr;
        int nys = 1;
/*
    int iy = 0;
    int nOKfVal = 3000;
    int nOKMol = 260;
    int ncomp = -1;

    bool *OKfVal = new bool[nOKfVal];
    double *xs = new double[nOKMol * nOKfVal];
    double *ys = new double[nOKMol];
    double *xsig = new double[3000];
    double *xmean = new double[3000];
    */

        double *resdl = new double[2 * nOKMol];
        double *xvresdl = new double[nOKMol];
        double q2 = -0.999;
        double qsdep = -0.999;
        int max20 = 20;
        int mxncomp = (ncomp != -1) ? ncomp :
                      sampls(max20, nOKMol, ys, nOKfVal, xs, &qsdep, &q2, xvresdl, sampall);
        if (mxncomp == 0)
            OEThrow.Fatal("SAMPLS failed. No model built\n");
        cout << "SAMPLS done\n";
        int maxiter = 100;
        double eps = 1.0e-4;
        int iboot = 0;
        int icros = 0;
        int icent = 2;
        double *weytb = new double[nOKMol];
        double *xbar = new double[nOKfVal + 1];
        double *xscal = new double[nOKfVal + 1];
        double *intcpt = new double[3];
        double *coeff = new double[2 * nOKfVal];
        double *varnc = new double[mxncomp * nOKfVal];
        double *wtx = new double[mxncomp * nOKfVal];
        double *wty = new double[mxncomp];
        double *loadings = new double[mxncomp * nOKfVal];
        double *latentx = new double[mxncomp * nOKMol];
        double *latenty = new double[mxncomp * nOKMol];
        double *inner = new double[mxncomp];
        double *ssqRsdl = new double[2];
        double *r2 = new double[2];
        double *sdep = new double[nOKfVal * mxncomp];
//   double *q2 = new double[2*mxncomp];
        int *optncomp = new int[2];
        double *ypred = new double[nOKMol];
        double *sss = new double[2];
        double *ssy = new double[nOKMol];
        int *iscratch = new int[nOKMol];
        double *scratch = new double[nOKMol * (nOKfVal + 2) + (nOKfVal + 1)];
        for (iy = 0; iy < nOKMol * (nOKfVal + 2) + (nOKfVal + 1); iy++) scratch[iy] = 0.0;
        for (iy = 0; iy < 3; ++iy) intcpt[iy] = 0.0;
        for (iy = 0; iy < 2; ++iy) {
            optncomp[iy] = 0;
            r2[iy] = sss[iy] = ssqRsdl[iy] = 0.0;
        }
        for (iy = 0; iy < nOKMol; ++iy) {
            iscratch[iy] = 0;
            ssy[iy] = weytb[iy] = ypred[iy] = 0.0;
            ywts[iy] = 1.0;
        }
        for (iy = 0; iy < 2 * nOKMol; ++iy) resdl[iy] = 0.0;
        for (iy = 0; iy < 2 * nOKfVal; ++iy) coeff[iy] = 0.0;
        for (iy = 0; iy < mxncomp * nOKfVal; ++iy)
            varnc[iy] = wtx[iy] = loadings[iy] = sdep[iy] = 0.0;
        for (iy = 0; iy < mxncomp; ++iy) wty[iy] = inner[iy] = 0.0;
        for (iy = 0; iy < mxncomp * nOKMol; ++iy) latentx[iy] = latenty[iy] = 0.0;
        int plsOK = plsjer(&nOKMol, &nOKfVal, &nys, &mxncomp, &maxiter, &eps, &iboot,
                           &icros, &icent, xs, ys, ywts, weytb, xbar, xscal, intcpt, coeff,
                           varnc, wtx, wty, loadings, latentx, latenty, inner, ssqRsdl, r2,
                           NULL, NULL, optncomp, ypred, resdl, sss, ssy, &ierr, scratch, iscratch);
        if (plsOK != 1) {
            cout << "PLS error: ID =" << ierr << endl;
            return 1;
        }
        cout << "Model built\n";
        for (iy = 0; iy < nOKfVal; iy++) coeff[iy] /= xsig[iy];
        for (iy = 0; iy < nOKfVal; iy++) intcpt[0] -= xmean[iy] * coeff[iy];

// expand model; write out the contours
        double *beta = new double[nFpt];
        int nf = 0;
        for (int j = 0; j < nFpt; j++)
            if (OKfVal[j]) {
                beta[j] = coeff[nf];
                nf++;
            }
            else beta[j] = 0.0;

        double *cntrSter = new double[npt];
        double *cntrElec = new double[npt];
        for (int j = 0; j < npt; j++) {
            cntrSter[j] = 0.0;
            cntrElec[j] = 0.0;
        }
        nf = 0;
//        for (int j = 0; j < npt; j++)
        for (int j = 0; j < nFpt; j++)
            if (OKfVal[j]) {
                if (j < npt)
                    cntrSter[j] = (coeff[nf] * sig[nf]);
                else
                    cntrElec[j - npt] = (coeff[nf] * sig[nf]);
                nf++;
            }

        savemodel(modelfn, nFpt, intcpt[0], beta, npts, lo, mxncomp, qsdep, q2, r2[0], ssqRsdl[0], ncomp, cntrSter, cntrElec);

        if (resids) wresids(xvresdl, resdl, nMol, OKmol);

        if (predict)
            mpredict(nFpt, predfn, intcpt[0], beta, ster, elec, npts, lo);
    }

    if (predict && !mkmodel)
        fpredict(modelfn, predfn);

    return 0;
}
