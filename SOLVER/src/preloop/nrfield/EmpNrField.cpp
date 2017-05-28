// EmpNrField.cpp
// created by Kuangdai on 13-May-2016 
// constant nr integer field

#include "EmpNrField.h"
#include "Geodesy.h"
#include <sstream>

EmpNrField::EmpNrField(bool useLucky, int nu_ref, int nu_min, 
    bool scaleS, bool scaleT, bool scaleD, 
    double powS,
    double factPI, double startT, double powT,
    double factD0, double startD, double endD): 
NrField(useLucky), mNuRef(nu_ref), mNuMin(nu_min),
mScaleS(scaleS), mScaleT(scaleT), mScaleD(scaleD),
mROuter(Geodesy::getROuter()), mPowS(powS), 
mFactPI(factPI), mStartT(startT), mPowT(powT),
mFactD0(factD0), mStartD(startD), mEndD(endD) {
    // nothing
}

int EmpNrField::getNrAtPoint(const RDCol2 &coords) const {
    double s = coords(0);
    double z = coords(1);
    double r, theta;
    Geodesy::rtheta(coords, r, theta);
    double d = mROuter - r;

    // reference value
    double nu = (double)mNuRef;
    // distance to axis
    if (mScaleS) {
        nu *= pow(s / mROuter, mPowS);
    }
    // epicentral distance
    if (mScaleT && theta > mStartT) {
        nu *= 1. + (mFactPI - 1.) * pow((theta - mStartT) / (pi - mStartT), mPowT);
    }
    // surface wave 
    if (mScaleD && d <= mEndD) {
        if (d <= mStartD) {
            nu *= mFactD0;
        } else {
            nu *= 1. + (mFactD0 - 1.) / (mStartD - mEndD) * (d - mEndD);
        }
    }
    // minimum value
    int nr = 2 * std::max(mNuMin, (int)ceil(nu)) + 1;
    if (nr <= 0) {
        throw std::runtime_error("EmpNrField::getNrAtPoint || Non-positive Nr.");
    }
    return nr;
}

std::string EmpNrField::verbose() const {
    std::stringstream ss;
    ss << "\n================= Fourier Expansion Order ==================" << std::endl;
    ss << "  Type                       =   Empirical" << std::endl;
    ss << "  Reference Order            =   " << mNuRef << std::endl;
    ss << "  Minimum Order              =   " << mNuMin << std::endl;
    ss << "  Scaled by" << std::endl;
    ss << "    a) Distance to Axis      =   " << (mScaleS ? "YES" : "NO") << std::endl;
    ss << "    b) Epicentral Distance   =   " << (mScaleT ? "YES" : "NO") << std::endl;
    ss << "    c) Surface Enhancement   =   " << (mScaleD ? "YES" : "NO") << std::endl;
    ss << "  Use FFTW Lucky Numbers     =   " << (mUseLuckyNumber ? "YES" : "NO") << std::endl;
    ss << "================= Fourier Expansion Order ==================\n" << std::endl;
    return ss.str();
}

