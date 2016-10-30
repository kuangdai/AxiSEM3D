// PreloopFFTW_N.h
// created by Kuangdai on 21-Sep-2016 
// perform FFT using fftw

#pragma once

#include <fftw3.h>
#include <vector>
#include "eigenp.h"

class PreloopFFTW {
public:
    // initialize plans
    static void initialize(int Nmax);
    // finalize plans
    static void finalize();
    
    // get input and output
    static RDColX &getR2C_RMat(int nr) {return sR2C_RMats[nr - 1];};
    static CDColX &getR2C_CMat(int nr) {return sR2C_CMats[nr - 1];};
    static RDColX &getC2R_RMat(int nr) {return sC2R_RMats[nr - 1];};    
    static CDColX &getC2R_CMat(int nr) {return sC2R_CMats[nr - 1];};
     
    // forward, real => complex
    static void computeR2C(int nr);
    // backward, complex => real
    static void computeC2R(int nr);
    
private:
    static int sNmax;
    static std::vector<fftw_plan> sR2CPlans;
    static std::vector<fftw_plan> sC2RPlans;
    static std::vector<RDColX> sR2C_RMats;
    static std::vector<CDColX> sR2C_CMats;
    static std::vector<RDColX> sC2R_RMats;
    static std::vector<CDColX> sC2R_CMats;
};
