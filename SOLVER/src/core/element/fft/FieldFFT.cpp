// FieldFFT.cpp
// created by Kuangdai on 19-May-2017 
// FFT between Fourier and physical spaces


#include "FieldFFT.h"
#include "SolverFFTW_N3.h"
#include "SolverFFTW_N6.h"
#include "SolverFFTW_N9.h"

void FieldFFT::transformF2P(const vec_ar3_CMatPP &uc, int Nr) {
    int Nu = Nr / 2;
    CMatXN3 &ucf = SolverFFTW_N3::getC2R_CMat(Nr);
    makeFlat<vec_ar3_CMatPP, CMatXN3>(uc, ucf, Nu);
    SolverFFTW_N3::computeC2R(Nr);
    // output to SolverFFTW_N3::getC2R_RMat(Nr);
}

void FieldFFT::transformF2P(const vec_ar6_CMatPP &uc, int Nr) {
    int Nu = Nr / 2;
    CMatXN6 &ucf = SolverFFTW_N6::getC2R_CMat(Nr);
    makeFlat<vec_ar6_CMatPP, CMatXN6>(uc, ucf, Nu);
    SolverFFTW_N6::computeC2R(Nr);
    // output to SolverFFTW_N6::getC2R_RMat(Nr);
}

void FieldFFT::transformF2P(const vec_ar9_CMatPP &uc, int Nr) {
    int Nu = Nr / 2;
    CMatXN9 &ucf = SolverFFTW_N9::getC2R_CMat(Nr);
    makeFlat<vec_ar9_CMatPP, CMatXN9>(uc, ucf, Nu);
    SolverFFTW_N9::computeC2R(Nr);
    // output to SolverFFTW_N9::getC2R_RMat(Nr);
}

void FieldFFT::transformP2F(vec_ar3_CMatPP &uc, int Nr) {
    int Nu = Nr / 2;
    // input from SolverFFTW_N3::getR2C_RMat(Nr);
    SolverFFTW_N3::computeR2C(Nr);
    CMatXN3 &ucf = SolverFFTW_N3::getR2C_CMat(Nr);
    makeStruct<vec_ar3_CMatPP, CMatXN3>(uc, ucf, Nu);
}

void FieldFFT::transformP2F(vec_ar6_CMatPP &uc, int Nr) {
    int Nu = Nr / 2;
    // input from SolverFFTW_N6::getR2C_RMat(Nr);
    SolverFFTW_N6::computeR2C(Nr);
    CMatXN6 &ucf = SolverFFTW_N6::getR2C_CMat(Nr);
    makeStruct<vec_ar6_CMatPP, CMatXN6>(uc, ucf, Nu);
}

void FieldFFT::transformP2F(vec_ar9_CMatPP &uc, int Nr) {
    int Nu = Nr / 2;
    // input from SolverFFTW_N9::getR2C_RMat(Nr);
    SolverFFTW_N9::computeR2C(Nr);
    CMatXN9 &ucf = SolverFFTW_N9::getR2C_CMat(Nr);
    makeStruct<vec_ar9_CMatPP, CMatXN9>(uc, ucf, Nu);
}



