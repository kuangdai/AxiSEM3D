// PreloopFFTW.cpp
// created by Kuangdai on 21-Sep-2016 
// perform FFT using fftw

#include "PreloopFFTW.h"

int PreloopFFTW::sNmax = 0;
std::vector<fftw_plan> PreloopFFTW::sR2CPlans;
std::vector<fftw_plan> PreloopFFTW::sC2RPlans;
std::vector<RDColX> PreloopFFTW::sR2C_RMats;
std::vector<CDColX> PreloopFFTW::sR2C_CMats;
std::vector<RDColX> PreloopFFTW::sC2R_RMats;
std::vector<CDColX> PreloopFFTW::sC2R_CMats;

void PreloopFFTW::initialize(int Nmax) {
    int xx = 1;
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
        sR2C_RMats.push_back(RDColX(NR, 1));
        sR2C_CMats.push_back(CDColX(NC, 1));
        sC2R_RMats.push_back(RDColX(NR, 1));
        sC2R_CMats.push_back(CDColX(NC, 1));
        double *r2c_r = &(sR2C_RMats[NR - 1](0, 0));
        ComplexD *r2c_c = &(sR2C_CMats[NR - 1](0, 0));
        sR2CPlans.push_back(fftw_plan_many_dft_r2c(
            1, n, xx, r2c_r, n, 1, NR, reinterpret_cast<fftw_complex*>(r2c_c), n, 1, NC, FFTW_ESTIMATE));   
        double *c2r_r = &(sC2R_RMats[NR - 1](0, 0));
        ComplexD *c2r_c = &(sC2R_CMats[NR - 1](0, 0));
        sC2RPlans.push_back(fftw_plan_many_dft_c2r(
            1, n, xx, reinterpret_cast<fftw_complex*>(c2r_c), n, 1, NC, c2r_r, n, 1, NR, FFTW_ESTIMATE)); 
    }
}

void PreloopFFTW::finalize() {
    for (int i = 0; i < sNmax; i++) {
        fftw_destroy_plan(sR2CPlans[i]);
        fftw_destroy_plan(sC2RPlans[i]);
    }
    sNmax = 0;
}

void PreloopFFTW::computeR2C(int nr) {
    fftw_execute(sR2CPlans[nr - 1]);
    double inv_nr = one / (double)nr;
    sR2C_CMats[nr - 1] *= inv_nr;
}

void PreloopFFTW::computeC2R(int nr) {
    fftw_execute(sC2RPlans[nr - 1]);
}
