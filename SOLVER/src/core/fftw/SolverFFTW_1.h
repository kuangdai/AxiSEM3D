// SolverFFTW_1.h
// created by Kuangdai on 21-Apr-2016 
// perform FFT using fftw

#pragma once

#include "SolverFFTW.h"

class SolverFFTW_1 {
public:
    // initialize plans
    static void initialize(int Nmax);
    // finalize plans
    static void finalize();
    
    // get input and output
    static RColX &getR2C_RMat(int nr) {return sR2C_RMats[nr - 1];};
    static CColX &getR2C_CMat(int nr) {return sR2C_CMats[nr - 1];};
    static RColX &getC2R_RMat(int nr) {return sC2R_RMats[nr - 1];};    
    static CColX &getC2R_CMat(int nr) {return sC2R_CMats[nr - 1];};
     
    // forward, real => complex
    static void computeR2C(int nr);
    // backward, complex => real
    static void computeC2R(int nr);
        
private:
    
    static int sNmax;
    static std::vector<PlanFFTW> sR2CPlans;
    static std::vector<PlanFFTW> sC2RPlans;
    static std::vector<RColX> sR2C_RMats;
    static std::vector<CColX> sR2C_CMats;
    static std::vector<RColX> sC2R_RMats;
    static std::vector<CColX> sC2R_CMats;
};
