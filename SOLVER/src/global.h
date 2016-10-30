// global.h
// created by Kuangdai on 27-Mar-2016 
// global constants and macros

#pragma once

// tiny number
const float  tinySingle = 1e-5f;
const double tinyDouble = 1e-10;

// solver precision
#ifdef _USE_DOUBLE
    typedef double Real;
    const Real tinyReal = tinyDouble;
#else
    typedef float Real;
    const Real tinyReal = tinySingle;
#endif

// complex
#include <complex>
typedef std::complex<Real>   Complex;
typedef std::complex<double> ComplexD;

// polynomial order
#ifndef _NPOL
    #define _NPOL 4
#endif

const int nPol = _NPOL;
const int nPntEdge = nPol + 1;
const int nPntElem = nPntEdge * nPntEdge;
const int nPE = nPntElem;

// constants
const Real zero     = (Real)0.0;
const Real one      = (Real)1.0;
const Real two      = (Real)2.0;
const Real three    = (Real)3.0;
const Real half     = (Real)0.5;
const Real third    = (Real)(1.0 / 3.0);
const Complex ii    = Complex(zero, one);
const Complex czero = Complex(zero, zero);
const double pi     = 3.141592653589793238463;
const double degree = pi / 180.;
const ComplexD iid  = ComplexD(0., 1.);

// shorten cwiseProduct 
#define schur cwiseProduct

// top-level source dir
const std::string projectDirectory = _PROJECT_DIR;
const std::string fftwWisdomDirectory = _FFTW_WISDOM_DIR;
