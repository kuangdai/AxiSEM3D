// SolverFFTW_N3.h
// created by Kuangdai on 21-Apr-2016 
// perform FFT using fftw

#pragma once

#include "SolverFFTW.h"

class SolverFFTW_N3 {
public:
    // initialize plans
    static void initialize(int Nmax);
    // finalize plans
    static void finalize();
    
    // get input and output
    static RMatXN3 &getR2C_RMat(int nr) {return sR2C_RMats[nr - 1];};
    static CMatXN3 &getR2C_CMat(int nr) {return sR2C_CMats[nr - 1];};
    static RMatXN3 &getC2R_RMat(int nr) {return sC2R_RMats[nr - 1];};    
    static CMatXN3 &getC2R_CMat(int nr) {return sC2R_CMats[nr - 1];};
     
    // forward, real => complex
    static void computeR2C(int nr);
    // backward, complex => real
    static void computeC2R(int nr);
        
private:
    static int sNmax;
    static std::vector<PlanFFTW> sR2CPlans;
    static std::vector<PlanFFTW> sC2RPlans;
    static std::vector<RMatXN3> sR2C_RMats;
    static std::vector<CMatXN3> sR2C_CMats;
    static std::vector<RMatXN3> sC2R_RMats;
    static std::vector<CMatXN3> sC2R_CMats;
};
