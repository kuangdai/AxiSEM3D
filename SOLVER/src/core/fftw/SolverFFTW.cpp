// SolverFFTW.h
// created by Kuangdai on 21-Apr-2016 
// perform FFT using fftw

#include "SolverFFTW.h"
#include "XMPI.h"

#include <fstream>
#include <sstream>

unsigned SolverFFTW::mWisdomLearnOption = FFTW_PATIENT;
bool SolverFFTW::mDisableWisdom = false; 

void SolverFFTW::importWisdom(bool disableWisdom) {
    mDisableWisdom = disableWisdom;
    if (mDisableWisdom) {
        mWisdomLearnOption = FFTW_ESTIMATE;
        return;
    }
    mWisdomLearnOption = FFTW_PATIENT;
    
    std::string wisdomstr = "";
    if (XMPI::root()) {
        #ifdef _USE_DOUBLE
            std::string fname = fftwWisdomDirectory + "/fftw_wisdom.double";
        #else
            std::string fname = fftwWisdomDirectory + "/fftw_wisdom.float";
        #endif
        std::ifstream fs(fname);
        std::stringstream buffer;
        if (fs) {
            buffer << fs.rdbuf();
            fs.close();
        }
        wisdomstr = buffer.str();
    }
    XMPI::bcast(wisdomstr);
    if (wisdomstr.length() > 0) {
        #ifdef _USE_DOUBLE
            fftw_import_wisdom_from_string(wisdomstr.c_str());
        #else
            fftwf_import_wisdom_from_string(wisdomstr.c_str());
        #endif
    }
}

void SolverFFTW::exportWisdom() {
    if (mDisableWisdom) {
        return;
    }
    
    if (XMPI::root()) {
        XMPI::mkdir(fftwWisdomDirectory);
        if (!XMPI::dirExists(fftwWisdomDirectory)) {
            throw std::runtime_error("SolverFFTW::exportWisdom || "
                "Error creating FFTW wisdom directory: ||" + fftwWisdomDirectory);    
        }
        #ifdef _USE_DOUBLE
            fftw_export_wisdom_to_filename((fftwWisdomDirectory + "/fftw_wisdom.double").c_str());
        #else
            fftwf_export_wisdom_to_filename((fftwWisdomDirectory + "/fftw_wisdom.float").c_str());
        #endif
    }
}

