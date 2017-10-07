// SolverFFTW_N3.cpp
// created by Kuangdai on 21-Apr-2016 
// perform FFT using fftw

#include "SolverFFTW_N3.h"

int SolverFFTW_N3::sNmax = 0;
std::vector<PlanFFTW> SolverFFTW_N3::sR2CPlans;
std::vector<PlanFFTW> SolverFFTW_N3::sC2RPlans;
std::vector<RMatXN3> SolverFFTW_N3::sR2C_RMats;
std::vector<CMatXN3> SolverFFTW_N3::sR2C_CMats;
std::vector<RMatXN3> SolverFFTW_N3::sC2R_RMats;
std::vector<CMatXN3> SolverFFTW_N3::sC2R_CMats;

void SolverFFTW_N3::initialize(int Nmax) {
    int ndim = 3;
    int xx = nPntElem * ndim;
    sNmax = Nmax;
    sR2CPlans.reserve(Nmax);
    sC2RPlans.reserve(Nmax);
    sR2C_RMats.reserve(Nmax);
    sR2C_CMats.reserve(Nmax);
    sC2R_RMats.reserve(Nmax);
    sC2R_CMats.reserve(Nmax);
    for (int NR = 1; NR <= Nmax; NR++) {
        int NC = NR / 2 + 1;
        int n[] = {NR};
        sR2C_RMats.push_back(RMatXN3(NR, xx));
        sR2C_CMats.push_back(CMatXN3(NC, xx));
        sC2R_RMats.push_back(RMatXN3(NR, xx));
        sC2R_CMats.push_back(CMatXN3(NC, xx));
        Real *r2c_r = &(sR2C_RMats[NR - 1](0, 0));
        Complex *r2c_c = &(sR2C_CMats[NR - 1](0, 0));
        sR2CPlans.push_back(planR2CFFTW(1, n, xx, r2c_r, n, 1, NR, complexFFTW(r2c_c), n, 1, NC, SolverFFTW::mWisdomLearnOption));   
        Real *c2r_r = &(sC2R_RMats[NR - 1](0, 0));
        Complex *c2r_c = &(sC2R_CMats[NR - 1](0, 0));
        sC2RPlans.push_back(planC2RFFTW(1, n, xx, complexFFTW(c2r_c), n, 1, NC, c2r_r, n, 1, NR, SolverFFTW::mWisdomLearnOption)); 
    }
}

void SolverFFTW_N3::finalize() {
    for (int i = 0; i < sNmax; i++) {
        distroyFFTW(sR2CPlans[i]);
        distroyFFTW(sC2RPlans[i]);
    }
    sNmax = 0;
}

void SolverFFTW_N3::computeR2C(int nr) {
    execFFTW(sR2CPlans[nr - 1]);
    Real inv_nr = one / (Real)nr;
    sR2C_CMats[nr - 1] *= inv_nr;
}

void SolverFFTW_N3::computeC2R(int nr) {
    execFFTW(sC2RPlans[nr - 1]);
}
