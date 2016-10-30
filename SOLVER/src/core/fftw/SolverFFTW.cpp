// SolverFFTW.h
// created by Kuangdai on 21-Apr-2016 
// perform FFT using fftw

#include "SolverFFTW.h"
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include "XMPI.h"

void SolverFFTW::importWisdom() {
    std::string wisdomstr;
    if (XMPI::root()) {
        #ifdef _USE_DOUBLE
            fftw_import_wisdom_from_filename((fftwWisdomDirectory + "/fftw_wisdom.double").c_str());
            char *wisdom = fftw_export_wisdom_to_string();
        #else
            fftwf_import_wisdom_from_filename((fftwWisdomDirectory + "/fftw_wisdom.float").c_str());
            char *wisdom = fftwf_export_wisdom_to_string();
        #endif
        wisdomstr = std::string(wisdom);
        free(wisdom);
    }
    XMPI::bcast(wisdomstr);
    if (!XMPI::root()) {
        #ifdef _USE_DOUBLE
            fftw_import_wisdom_from_string(wisdomstr.c_str());
        #else
            fftwf_import_wisdom_from_string(wisdomstr.c_str());
        #endif
    }
}

void SolverFFTW::exportWisdom() {
    if (XMPI::root()) {
        if (!boost::filesystem::exists(fftwWisdomDirectory)) 
            boost::filesystem::create_directory(fftwWisdomDirectory);
        if (!boost::filesystem::exists(fftwWisdomDirectory)) 
            throw std::runtime_error("SolverFFTW::exportWisdom || Error creating FFTW wisdom directory: ||" + fftwWisdomDirectory);    
        #ifdef _USE_DOUBLE
            fftw_export_wisdom_to_filename((fftwWisdomDirectory + "/fftw_wisdom.double").c_str());
        #else
            fftwf_export_wisdom_to_filename((fftwWisdomDirectory + "/fftw_wisdom.float").c_str());
        #endif
    }
}