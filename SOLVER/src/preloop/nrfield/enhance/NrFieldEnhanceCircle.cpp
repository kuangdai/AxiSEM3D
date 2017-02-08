// NrFieldEnhanceCircle.cpp
// created by Kuangdai on 2-Feb-2017 
// circle-shaped enhanced nr integer field

#include "NrFieldEnhanceCircle.h"
#include <sstream>

NrFieldEnhanceCircle::NrFieldEnhanceCircle(int ref, bool decrease, 
    double r, double theta, double diameter, double hwhm, double value):
NrFieldEnhance(ref, decrease), mR(r), mTheta(theta), 
mDiameter(diameter), mHWHM(hwhm), mValue(value) {
    // nothing    
}

std::string NrFieldEnhanceCircle::verbose() const {
    std::stringstream ss;
    ss << "\n=================== Local Nu Enhancement ===================" << std::endl;
    ss << "  Type                    =   Circle" << std::endl;
    ss << "  Location (r, theta)     =   " <<  "(" << mR / 1e3 << ", " << mTheta / degree << ")" << std::endl;
    ss << "  Diameter / km           =   " << mDiameter / 1e3 << std::endl;
    ss << "  HWHM / km               =   " << mHWHM / 1e3 << std::endl;
    if (mReference == 0) {
        ss << "  Nu (factor) at center   =   " << (int)(mValue / 2) << std::endl;
        ss << "  Reference Type          =   " << "Absolute" << std::endl;
    } else {
        ss << "  Nu (factor) at center   =   " << mValue << std::endl;
        ss << "  Reference Type          =   " << (mReference == 1 ? "Reference Base" : "Reference Current") << std::endl;
    }
    ss << "  Decrease Original Nu    =   " << (mDecrease ? "YES" : "NO") << std::endl;
    ss << "=================== Local Nu Enhancement ===================\n" << std::endl;
    return ss.str();
}

double NrFieldEnhanceCircle::getValue(const RDCol2 &sz_target) const {
    RDCol2 sz_center;
    sz_center(0) = mR * sin(mTheta);
    sz_center(1) = mR * cos(mTheta);
    double distance = (sz_center - sz_target).norm();
    
    // treat as center if inside bubble
    distance -= mDiameter * .5; 
    if (distance < 0.) distance = 0.;
    
    // outside range
    if (distance > 4. * mHWHM) return 0.;
    
    // compute Gaussian
    double stddev = mHWHM / sqrt(2. * log(2.));
    return mValue * exp(-distance * distance / (stddev * stddev * 2.));
}

