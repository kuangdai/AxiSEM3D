// Volumetric3D_bubble.cpp
// created by Kuangdai on 19-Oct-2016 
// a bubble-shaped heterogeneity

#include "Volumetric3D_bubble.h"
#include "XMath.h"
#include <sstream>

void Volumetric3D_bubble::initialize(const std::vector<double> &params) {
    // need at least 6 parameters to make a bubble
    if (params.size() < 6) throw std::runtime_error("Volumetric3D_bubble::initialize || "
        "Not enough parameters to initialize a Volumetric3D_bubble object.");
    
    // initialize location    
    mRadius = params[0] * 1e3;
    mTheta = params[1] * degree;
    mPhi = params[2] * degree;
    
    // initialize Gaussian parameters
    mMax = params[3];
    mHWHM = params[4] * 1e3;
    
    // initialize reference type
    if (params[5] < 0.5) 
        mReferenceType = ReferenceTypes::Absolute;
    else if (params[5] < 1.5) 
        mReferenceType = ReferenceTypes::Reference1D;
    else if (params[5] < 2.5)
        mReferenceType = ReferenceTypes::Reference3D;
    else 
        mReferenceType = ReferenceTypes::ReferenceDiff;    
        
    // initialize optional flags    
    if (params.size() >= 7) mChangeVp = (params[6] > 0.);
    if (params.size() >= 8) mChangeVs = (params[7] > 0.);
    if (params.size() >= 9) mChangeRho = (params[8] > 0.);
}

bool Volumetric3D_bubble::get3dProperties(double r, double theta, double phi, double rElemCenter,
    double &dvpv, double &dvph, double &dvsv, double &dvsh, double &drho) const {
    
    // find the distance between bubble center and target point
    RDCol3 rtpBubble, rtpTarget;
    rtpBubble(0) = mRadius;
    rtpBubble(1) = mTheta;
    rtpBubble(2) = mPhi;
    const RDCol3 &xyzBubble = XMath::toCartesian(rtpBubble);
    rtpTarget(0) = r;
    rtpTarget(1) = theta;
    rtpTarget(2) = phi;
    const RDCol3 &xyzTarget = XMath::toCartesian(rtpTarget);
    double distance = (xyzBubble - xyzTarget).norm();
    
    // zero results
    dvpv = dvph = dvsv = dvsh = drho = 0.;
    
    // outside range
    if (distance > 4. * mHWHM) return false;
    
    // compute Gaussian
    double stddev = mHWHM / sqrt(2. * log(2.));
    double gaussian = mMax * exp(-distance * distance / (stddev * stddev * 2.));
    
    // set perturbations    
    if (mChangeVp) dvpv = dvph = gaussian;
    if (mChangeVs) dvsv = dvsh = gaussian;
    if (mChangeRho) drho = gaussian;
    
    return true;
}

std::string Volumetric3D_bubble::verbose() const {
    std::stringstream ss;
    ss << "\n======================= 3D Volumetric ======================" << std::endl;
    ss << "  Model Name       =   bubble" << std::endl;
    ss << "  Radius / km      =   " << mRadius / 1e3 << std::endl;
    ss << "  Theta / degree   =   " << mTheta / degree << std::endl;
    ss << "  Phi / degree     =   " << mPhi / degree << std::endl;
    ss << "  Maximum          =   " << mMax << std::endl;
    ss << "  HWHM / km        =   " << mHWHM / 1e3 << std::endl;
    ss << "  Reference Type   =   " << ReferenceTypesString[mReferenceType] << std::endl;
    ss << "  Affect VP        =   " << (mChangeVp ? "YES" : "NO") << std::endl;
    ss << "  Affect VS        =   " << (mChangeVs ? "YES" : "NO") << std::endl;
    ss << "  Affect Density   =   " << (mChangeRho ? "YES" : "NO") << std::endl;
    ss << "======================= 3D Volumetric ======================\n" << std::endl;
    return ss.str();
}

