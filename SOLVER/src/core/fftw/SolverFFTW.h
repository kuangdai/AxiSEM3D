// SolverFFTW.h
// created by Kuangdai on 21-Apr-2016 
// perform FFT using fftw

#pragma once

#include <fftw3.h>
#include <vector>
#include "eigenc.h"

#ifdef _USE_DOUBLE
    typedef fftw_plan PlanFFTW;
    #define complexFFTW reinterpret_cast<fftw_complex*>
    #define planR2CFFTW fftw_plan_many_dft_r2c
    #define planC2RFFTW fftw_plan_many_dft_c2r
    #define distroyFFTW fftw_destroy_plan
    #define execFFTW fftw_execute
#else
    typedef fftwf_plan PlanFFTW;
    #define complexFFTW reinterpret_cast<fftwf_complex*>
    #define planR2CFFTW fftwf_plan_many_dft_r2c
    #define planC2RFFTW fftwf_plan_many_dft_c2r
    #define distroyFFTW fftwf_destroy_plan
    #define execFFTW fftwf_execute
#endif

class SolverFFTW {
public:
    static void importWisdom(bool disableWisdom);
    static void exportWisdom();    
    static unsigned mWisdomLearnOption;
    
private:    
    static bool mDisableWisdom; 
};
