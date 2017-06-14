#include <fstream>
#include <iostream>
#include <cmath>

double vsd(int n, double *a) {
// ***********************************
    double t = 0.0;
    for (int i = 0; i < n; ++i) t += a[i];
    t /= ((double) n);
    double s = 0.0;
    for (int i = 0; i < n; ++i) s += (a[i] - t) * (a[i] - t);
    return sqrt(s / ((double) (n - 1)));
}

void ccfit(int maxnc, int nMol, double *cc, double *beta,
           double *vcoef, double *work, double *hwork, double *ycmpnt,
           double *bio, double *y1, double *wy) {
// **********************************************
    double pac = 1.0e-08;
    double pacd = sqrt((double) (nMol - 1)) * pac;
    double pacn = vsd(nMol, bio) * pacd;

    for (int i = 0; i < nMol * maxnc; i++) vcoef[i] = 0.0;

    double wt = (double) (nMol - 1);
    double a = 0.0;
    for (int i = 0; i < nMol; ++i) a += wy[i] * bio[i];
    *y1 = a / wt;

    double a2 = 0.0;
    for (int i = 0; i < nMol; ++i) {
        work[i * 3] = bio[i] - *y1;
        a2 += wy[i] * work[i * 3] * work[i * 3];
    }
//	double y2 = a2/(wt-1.0);

    for (int nc = 0; nc < maxnc; ++nc) {
        for (int i = 0; i < nMol; ++i) {
            work[i * 3 + 2] = work[i * 3];
            vcoef[i * maxnc + nc] = work[i * 3 + 2];
        }
        double yy = 0.0;
        for (int i = 0; i < nMol; ++i) yy += work[i * 3 + 2] * work[i * 3 + 2] * wy[i];
        for (int j = 0; j < nMol; ++j) {
            a = 0.0;
            for (int i = 0; i < nMol; ++i)
                a += cc[i * nMol + j] * work[i * 3 + 2] * wy[i];
            work[j * 3 + 1] = a;
        }
// Center "s"=w2, i.e. orthogonalize to any constant-vector
        a = 0.0;
        for (int i = 0; i < nMol; ++i) a += work[i * 3 + 1] * wy[i];
        a /= wt;
        for (int i = 0; i < nMol; ++i) work[i * 3 + 1] -= a;
// Orthogonalize "s"=w2 to all previous t(jh) 
        for (int ic = 0; ic < nc; ic++) {
            a = 0.0;
            for (int i = 0; i < nMol; ++i) a += work[i * 3 + 1] * ycmpnt[i * maxnc + ic] * wy[i];
            //(use hwork(jh),the stored value of ycmpnt**2 for prev ic)
            double f = a / hwork[ic];
            for (int i = 0; i < nMol; ++i) work[i * 3 + 1] -= f * ycmpnt[i * maxnc + ic];
            for (int i = 0; i < nMol; ++i) vcoef[i * maxnc + nc] -= f * vcoef[i * maxnc + ic];
        }
        // Set (t)=(s)
        double tt = 0.0;
        double ty = 0.0;
        for (int i = 0; i < nMol; ++i) {
            double w = work[i * 3 + 1];
            tt += w * w * wy[i];
        }
        double eig = sqrt(tt / yy);
        double del = 0.0;
        for (int i = 0; i < nMol; i++) {
            double f = work[i * 3 + 1] / eig;
            del += wy[i] * fabs(f - work[i * 3 + 2]);
            work[i * 3 + 2] = f;
        }
        tt = ty = 0.0;
        for (int i = 0; i < nMol; i++) {
            double w = work[i * 3 + 1];
            tt += w * w * wy[i];
            ty += wy[i] * w * work[i * 3];
        }

        double r = (ty + pacn) / (tt + pacd);
        beta[nc] = r;
        for (int j = 0; j < nMol; j++)
            ycmpnt[j * maxnc + nc] = r * work[j * 3 + 1];
        for (int i = 0; i < nMol; ++i)
            vcoef[i * maxnc + nc] = r * vcoef[i * maxnc + nc];
        a2 = 0.0;
        for (int i = 0; i < nMol; ++i) {
            a2 += ycmpnt[i * maxnc + nc] * ycmpnt[i * maxnc + nc] * wy[i];
            work[i * 3] -= ycmpnt[i * maxnc + nc];
        }
        hwork[nc] = a2;
    }
}


void getCovar(const int nfld, double *fld, double *Covar, int nMol, int step) {
// **********************************************
    double *dists = NULL;
    if (nfld == 0)
        dists = fld;
    else {
//		ofstream ofs("12diff.txt");
        dists = new double[nMol * nMol];
        double d1;
        double d2;
        double sum;
        for (int i = 0; i < nMol; i++) {
            dists[i * nMol + i] = 0.0;  //diagonal
            for (int j = 0; j < i; j++) {
                sum = 0.0;
                for (int fv = 0; fv < nfld; fv++) {
                    d1 = fld[i * nfld * step + fv];
                    d2 = fld[j * nfld * step + fv];
                    sum += (d1 - d2) * (d1 - d2);
//			if (j==0 && i==1) ofs << fv << ' ' << d1 << ' ' << d2 << ' ' << (d1 - d2) * (d1 - d2) << ' ' << endl;
                }
                dists[i * nMol + j] = dists[j * nMol + i] = sum;
            }
        }
    }

    double s = 0.0;
    for (int i = 1; i < nMol; i++) for (int j = 0; j < i; j++) s += dists[i * nMol + j];
    double s2 = s;
    double s1 = s2 / ((double) nMol);
    for (int i = 0; i < nMol; i++) {
        double t = -s1;
        for (int j = 0; j < nMol; j++) t += dists[i * nMol + j];
        Covar[i * nMol + i] = t / ((double) nMol);
    }
    for (int i = 1; i < nMol; i++)
        for (int j = 0; j < i; j++) {
            double t = 0.5 * (Covar[i * nMol + i] + Covar[j * nMol + j] - dists[i * nMol + j]);
            Covar[i * nMol + j] = Covar[j * nMol + i] = t;
        }
    delete dists;
}

int sampls(int maxnc, int manyMol, double *manybio, int nfld,
           double *fld, double *sdep, double *q2, double *res, bool doall) {
// **************************************************

//   cout << "into sampls\n";
    int nMol = 1000;
    int step = manyMol / nMol;
    if (manyMol < 1000 || doall) {
        nMol = manyMol;
        step = 1;
    }
    double *Covar = new double[nMol * nMol];
    getCovar(nfld, fld, Covar, nMol, step);

    double *bio = new double[nMol];
    for (int i = 0; i < nMol; i++) bio[i] = manybio[i * step];

    double *ycmpnt = new double[nMol * maxnc];
    double wy[nMol];
    double *yrwork = new double[nMol * maxnc];
    double *beta = new double[maxnc];
    double *zwork = new double[nMol * maxnc];
    double *work = new double[3 * nMol];
    double *hwork = new double[maxnc];
    double *vcoeff = new double[nMol * maxnc];
    for (int i = 0; i < nMol * maxnc; i++) vcoeff[i] = 0.0;
    double y1work = 0.0;

    double avgBio = 0.0;
    for (int nm = 0; nm < nMol; nm++) avgBio += bio[nm];
    avgBio /= (double) nMol;

    for (int lxv = 0; lxv < nMol; ++lxv) {
        if ((lxv + 1) % 100 == 0) std::cout << " Processing compound " << lxv + 1 << std::endl;
        for (int i = 0; i < nMol; ++i) wy[i] = 1.0;
        wy[lxv] = 0.0;

        ccfit(maxnc, nMol, Covar, beta, zwork, work, hwork, yrwork, bio, &y1work, wy);

        for (int ic = 0; ic < maxnc; ++ic) {
            ycmpnt[lxv * maxnc + ic] = yrwork[lxv * maxnc + ic];
            if (ic == 0) ycmpnt[lxv * maxnc] += (y1work - avgBio);
            for (int i = 0; i < nMol; i++)
                if (i != lxv)
                    vcoeff[i * maxnc + ic] += zwork[i * maxnc + ic] / ((double) nMol);
        }
    }

    double a = 0.0;
    for (int nm = 0; nm < nMol; nm++) a += bio[nm];
    double y1 = a / ((double) nMol);
    double y2 = 0.0;
    for (int nm = 0; nm < nMol; nm++) y2 += (bio[nm] - y1) * (bio[nm] - y1);

    double *pwork = new double[2 * maxnc];
    for (int nc = 0; nc < maxnc; nc++) pwork[nc * 2 + 1] = 0.0;
    for (int nm = 0; nm < nMol; nm++) {
        double dy = bio[nm] - y1;
        for (int nc = 0; nc < maxnc; ++nc) {
            dy -= ycmpnt[nm * maxnc + nc];
            pwork[nc * 2] = dy;
            pwork[nc * 2 + 1] += dy * dy;
        }
    }
    int optnc = 0;
    double optsd = sqrt(y2 / ((double) (nMol - 1)));
    double varBio = optsd;
    for (int nc = 0; nc < maxnc; ++nc) {
        int nDegF = nMol - nc - 2;
        nDegF = nDegF > 0 ? nDegF : 1;
        double nowsd = sqrt(pwork[nc * 2 + 1] / ((double) (nDegF)));
        if (nowsd > optsd) break;
        if (nowsd < 0.2 && nowsd / varBio < 0.2) break;
        optsd = nowsd;
        optnc = nc + 1;
    }
    *sdep = optsd;
    *q2 = 1.0 - optsd * optsd / (varBio * varBio);

    for (int nm = 0; nm < nMol; nm++) {
        double biodiff = bio[nm] - avgBio;
        for (int nc = 0; nc < optnc; ++nc) {
            biodiff -= ycmpnt[nm * maxnc + nc];
            pwork[nc] = biodiff;
        }
        res[nm] = pwork[optnc - 1];
    }
    return optnc;
}

/*
int main() {

  double bio[14] = {-2.15,-1.28,-1.19,-1.0,-.75,-.63,-.42,-.40,-.15,
      -.05,.02,.04,.35,.40};
      double fld[42] = {-.61,4.4,9.63, -.61,6.0,11.0, 0,3.39,9.67, 0,3.3,8.34,
      -.05,2.243,10.321, .42,3.197,8.836, -.17,3.95,8.8, -.24,4.233,9.417,
      -.4,5.849,11.216, -.2,5.119,12.085, .23,4.21,9.113, .23,4.06,8.828,
      .03,4.655,10.07, .54,4.27,8.847};

      double *res = new double[14];
      double *xsig = new double[3];
      double *xmean = new double[3];
      for (int i = 0; i < 3; i++) {
          double sum = 0.0;
          double ssq = 0.0;
          for (int j = 0; j < 14; j++) {
              sum += fld[j*3 + i];
              ssq += fld[j*3 + i] * fld[j*3 + i];
          }
          xsig[i] = sqrt((ssq - (sum*sum/14.0) ) / 13.0);
      }

    double sdep;
    double q2;
    std::ifstream fbio{"fort.1"};
    std::ifstream fdist{"fort.3"};
    std::string s;
    int ncomp = 0;
    int nmol = 0;
    fbio >> s >> ncomp;
    fbio >> s >> nmol;
    double *bio = new double[nmol];
    double *res = new double[nmol];
    for (int i = 0; i < nmol + 8; i++) fbio >> s;
    for (int ct = 0; ct < nmol; ct++)
        fbio >> s >> bio[ct];

    double *dists = new double[nmol * nmol];
    for (int i = 0; i < nmol * nmol; i++) dists[i] = 0.0;
    for (int i = 0; i < 5; i++) fdist >> s;
    int i;
    int j;
    double d;
    while (true) {
        fdist >> i >> j >> d;
        i--;
        j--;
        d *= d;
        dists[i * nmol + j] = d;
        dists[j * nmol + i] = d;
        if (i == nmol - 1 && j == i - 1) break;
    }
    std::cout << sampls(ncomp, nmol, bio, 0, dists, &sdep, &q2, res, true);

//    std::cout << sampls(4, 14, bio, 3, fld, &sdep, &q2, res, xsig );
}
 */


