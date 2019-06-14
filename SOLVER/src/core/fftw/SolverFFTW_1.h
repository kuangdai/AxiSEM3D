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
    static RColX &getR2C_RMat() {return sR2C_RMat;};
    static CColX &getR2C_CMat() {return sR2C_CMat;};
    static RColX &getC2R_RMat() {return sC2R_RMat;};    
    static CColX &getC2R_CMat() {return sC2R_CMat;};
     
    // forward, real => complex
    static void computeR2C(int nr);
    // backward, complex => real
    static void computeC2R(int nr);
        
private:
    static int sNmax;
    static std::vector<PlanFFTW> sR2CPlans;
    static std::vector<PlanFFTW> sC2RPlans;
    static RColX sR2C_RMat;
    static CColX sR2C_CMat;
    static RColX sC2R_RMat;
    static CColX sC2R_CMat;
};
