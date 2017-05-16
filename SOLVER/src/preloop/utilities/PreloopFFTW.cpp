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

void PreloopFFTW::checkAndInit(int nr) {
    if (nr <= sNmax) {
        return;
    }
    int xx = 1;
    for (int NR = sNmax + 1; NR <= nr; NR++) {
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
    sNmax = nr;
}

void PreloopFFTW::finalize() {
    for (int i = 0; i < sNmax; i++) {
        fftw_destroy_plan(sR2CPlans[i]);
        fftw_destroy_plan(sC2RPlans[i]);
    }
    sNmax = 0;
}

void PreloopFFTW::computeR2C(int nr) {
    checkAndInit(nr); 
    fftw_execute(sR2CPlans[nr - 1]);
    double inv_nr = one / (double)nr;
    sR2C_CMats[nr - 1] *= inv_nr;
}

void PreloopFFTW::computeC2R(int nr) {
    checkAndInit(nr); 
    fftw_execute(sC2RPlans[nr - 1]);
}

bool PreloopFFTW::isLuckyNumber(int n, bool forceOdd)
{
    int num = n;
    
    // We always hope to use even numbers that are generally faster,
    // but the Nyquist frequency sometimes causes trouble.
    // force odd
    if (forceOdd && num % 2 == 0) {
        return false;
    }
    
    // use even when n > 10
    if (!forceOdd && num % 2 != 0 && num > 10) {
        return false;
    }
    
    for (int i = 2; i <= num; i++) {  
        while(num % i == 0) {
            num /= i;
            if (i > 13) {
                return false;
            }
        }
    }
    num = n;
    int e = 0;
    while(num % 11 == 0) {
        num /= 11;
        e++;
    }
    num = n;
    int f = 0;
    while(num % 13 == 0) {
        num /= 13;
        f++;
    }
    if (e + f > 1) {
        return false;
    }
    return true;
}

int PreloopFFTW::nextLuckyNumber(int n, bool forceOdd)
{
    while(true) {
        if (isLuckyNumber(n, forceOdd)) {
            return n;
        }
        n++;
    }
}
