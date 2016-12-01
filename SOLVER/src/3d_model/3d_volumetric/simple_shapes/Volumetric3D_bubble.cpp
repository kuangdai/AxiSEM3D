// Volumetric3D_bubble.cpp
// created by Kuangdai on 19-Oct-2016 
// a bubble-shaped heterogeneity

#include "Volumetric3D_bubble.h"
#include "XMath.h"
#include <sstream>

void Volumetric3D_bubble::initialize(const std::vector<std::string> &params) {
    // need at least 6 parameters to make a bubble
    if (params.size() < 7) throw std::runtime_error("Volumetric3D_bubble::initialize || "
        "Not enough parameters to initialize a Volumetric3D_bubble object.");
    
    const std::string source = "Volumetric3D_bubble::initialize";
    
    // initialize location 
    XMath::castValue(mDepth, params[0], source); mDepth *= 1e3;
    XMath::castValue(mLat, params[1], source);
    XMath::castValue(mLon, params[2], source);
    
    // initialize Gaussian parameters
    XMath::castValue(mRadius, params[3], source); mRadius *= 1e3;
    XMath::castValue(mHWHM, params[4], source); mHWHM *= 1e3;
    XMath::castValue(mMax, params[5], source);
    
    // initialize reference type
    if (boost::iequals(params[6], "Absolute") || boost::iequals(params[6], "Abs")) 
        mReferenceType = ReferenceTypes::Absolute;
    else if (boost::iequals(params[6], "Reference1D") || boost::iequals(params[6], "Ref1D") || boost::iequals(params[6], "1D")) 
        mReferenceType = ReferenceTypes::Reference1D;
    else if (boost::iequals(params[6], "Reference3D") || boost::iequals(params[6], "Ref3D") || boost::iequals(params[6], "3D")) 
        mReferenceType = ReferenceTypes::Reference3D;
    else if (boost::iequals(params[6], "ReferencePerturb") || boost::iequals(params[6], "RefPerturb") || boost::iequals(params[6], "Perturb")) 
        mReferenceType = ReferenceTypes::ReferencePerturb;        
    else 
        throw std::runtime_error("Volumetric3D_bubble::initialize || Unknown reference type: " + params[6] + ".");
        
    try {
        int ipar = 7;
        XMath::castValue(mChangeVp, params[ipar++], source);
        XMath::castValue(mChangeVs, params[ipar++], source);
        XMath::castValue(mChangeRho, params[ipar++], source);
    } catch (std::out_of_range) {
        // nothing
    }    
}

bool Volumetric3D_bubble::get3dProperties(double r, double theta, double phi, double rElemCenter,
    double &dvpv, double &dvph, double &dvsv, double &dvsh, double &drho) const {
    
    // find the distance between bubble center and target point
    RDCol3 rtpBubble, rtpTarget;
    rtpBubble(0) = XMath::getROuter() - mDepth;
    rtpBubble(1) = XMath::lat2Theta(mLat, mDepth);
    rtpBubble(2) = XMath::lon2Phi(mLon);
    const RDCol3 &xyzBubble = XMath::toCartesian(rtpBubble);
    rtpTarget(0) = r;
    rtpTarget(1) = theta;
    rtpTarget(2) = phi;
    const RDCol3 &xyzTarget = XMath::toCartesian(rtpTarget);
    double distance = (xyzBubble - xyzTarget).norm();
    
    // zero results
    dvpv = dvph = dvsv = dvsh = drho = 0.;
    
    // treat as center if inside bubble
    distance -= mRadius; 
    if (distance < 0.) distance = 0.;
    
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
    ss << "  Model Name          =   bubble" << std::endl;
    ss << "  Depth / km          =   " << mDepth / 1e3 << std::endl;
    ss << "  Lat / degree        =   " << mLat << std::endl;
    ss << "  Lon / degree        =   " << mLon << std::endl;
    ss << "  Bubble Radius / km  =   " << mRadius / 1e3 << std::endl;
    ss << "  HWHM / km           =   " << mHWHM / 1e3 << std::endl;
    ss << "  Maximum at Center   =   " << mMax << std::endl;
    ss << "  Reference Type      =   " << ReferenceTypesString[mReferenceType] << std::endl;
    ss << "  Affect VP           =   " << (mChangeVp ? "YES" : "NO") << std::endl;
    ss << "  Affect VS           =   " << (mChangeVs ? "YES" : "NO") << std::endl;
    ss << "  Affect Density      =   " << (mChangeRho ? "YES" : "NO") << std::endl;
    ss << "======================= 3D Volumetric ======================\n" << std::endl;
    return ss.str();
}

