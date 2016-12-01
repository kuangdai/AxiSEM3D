// Volumetric3D_cylinder.cpp
// created by Kuangdai on 19-Oct-2016 
// a cylinder-shaped heterogeneity

#include "Volumetric3D_cylinder.h"
#include "XMath.h"
#include <sstream>

void Volumetric3D_cylinder::initialize(const std::vector<std::string> &params) {
    // need at least 6 parameters to make a cylinder
    if (params.size() < 11) throw std::runtime_error("Volumetric3D_cylinder::initialize || "
        "Not enough parameters to initialize a Volumetric3D_cylinder object.");
        
    const std::string source = "Volumetric3D_cylinder::initialize";
    
    // initialize location    
    XMath::castValue(mD1, params[0], source); mD1 *= 1e3;
    XMath::castValue(mLat1, params[1], source);
    XMath::castValue(mLon1, params[2], source);

    XMath::castValue(mD2, params[3], source); mD2 *= 1e3;
    XMath::castValue(mLat2, params[4], source);
    XMath::castValue(mLon2, params[5], source);
    
    // initialize Gaussian parameters
    XMath::castValue(mRadius, params[6], source); mRadius *= 1e3;
    XMath::castValue(mHWHM_lateral, params[7], source); mHWHM_lateral *= 1e3;
    XMath::castValue(mHWHM_top_bot, params[8], source); mHWHM_top_bot *= 1e3;
    XMath::castValue(mMaxAxis, params[9], source);
    
    // initialize reference type
    if (boost::iequals(params[10], "Absolute") || boost::iequals(params[10], "Abs")) 
        mReferenceType = ReferenceTypes::Absolute;
    else if (boost::iequals(params[10], "Reference1D") || boost::iequals(params[10], "Ref1D") || boost::iequals(params[10], "1D")) 
        mReferenceType = ReferenceTypes::Reference1D;
    else if (boost::iequals(params[10], "Reference3D") || boost::iequals(params[10], "Ref3D") || boost::iequals(params[10], "3D")) 
        mReferenceType = ReferenceTypes::Reference3D;
    else if (boost::iequals(params[10], "ReferencePerturb") || boost::iequals(params[10], "RefPerturb") || boost::iequals(params[10], "Perturb")) 
        mReferenceType = ReferenceTypes::ReferencePerturb;        
    else 
        throw std::runtime_error("Volumetric3D_bubble::initialize || Unknown reference type: " + params[10] + ".");
        
    try {
        int ipar = 11;
        XMath::castValue(mChangeVp, params[ipar++], source);
        XMath::castValue(mChangeVs, params[ipar++], source);
        XMath::castValue(mChangeRho, params[ipar++], source);
    } catch (std::out_of_range) {
        // nothing
    }    
}

bool Volumetric3D_cylinder::get3dProperties(double r, double theta, double phi, double rElemCenter,
    double &dvpv, double &dvph, double &dvsv, double &dvsh, double &drho) const {
    // zero results
    dvpv = dvph = dvsv = dvsh = drho = 0.;
    
    // distance from point to axis
    RDCol3 rtpPoint1, rtpPoint2, rtpTarget;
    rtpPoint1(0) = XMath::getROuter() - mD1;
    rtpPoint1(1) = XMath::lat2Theta(mLat1, mD1);
    rtpPoint1(2) = XMath::lon2Phi(mLon1);
    rtpPoint2(0) = XMath::getROuter() - mD2;
    rtpPoint2(1) = XMath::lat2Theta(mLat2, mD2);
    rtpPoint2(2) = XMath::lon2Phi(mLon2);
    rtpTarget(0) = r;
    rtpTarget(1) = theta;
    rtpTarget(2) = phi;
    const RDCol3 &xyzPoint1 = XMath::toCartesian(rtpPoint1);
    const RDCol3 &xyzPoint2 = XMath::toCartesian(rtpPoint2);
    const RDCol3 &xyzTarget = XMath::toCartesian(rtpTarget);
    const RDCol3 &xyzDiff1 = xyzTarget - xyzPoint1;
    const RDCol3 &xyzDiff2 = xyzTarget - xyzPoint2; 
    double length = (xyzPoint1 - xyzPoint2).norm();
    double distToLine = xyzDiff1.cross(xyzDiff2).norm() / length;
    
    // outside range
    if (distToLine > mRadius + 4. * mHWHM_lateral) return false;
    
    // compute Gaussian lateral
    double distToSurf = distToLine - mRadius;
    if (distToSurf < 0.) distToSurf = 0.;
    double stddev = mHWHM_lateral / sqrt(2. * log(2.));
    double gaussian_lateral = mMaxAxis * exp(-distToSurf * distToSurf / (stddev * stddev * 2.));
    
    // distance from point to ends
    double d1 = xyzDiff1.norm();
    double d2 = xyzDiff2.norm();
    double dmax = std::max(d1, d2);
    double distToTopBot = sqrt(dmax * dmax - distToLine * distToLine) - length;
    
    // outside range
    if (distToTopBot > 4. * mHWHM_top_bot) return false;
    
    double gaussian_topbot = gaussian_lateral;
    if (distToTopBot > 0.) {
        // compute Gaussian top-bottom
        stddev = mHWHM_top_bot / sqrt(2. * log(2.));
        gaussian_topbot = gaussian_lateral * exp(-gaussian_topbot * gaussian_topbot / (stddev * stddev * 2.)); 
    }
    
    // set perturbations    
    if (mChangeVp) dvpv = dvph = gaussian_topbot;
    if (mChangeVs) dvsv = dvsh = gaussian_topbot;
    if (mChangeRho) drho = gaussian_topbot;
    
    return true;
}

std::string Volumetric3D_cylinder::verbose() const {
    std::stringstream ss;
    ss << "\n======================= 3D Volumetric ======================" << std::endl;
    ss << "  Model Name             =   cylinder" << std::endl;
    ss << "  Depth_1 / km           =   " << mD1 / 1e3 << std::endl;
    ss << "  Lat_1 / degree         =   " << mLat1 << std::endl;
    ss << "  Lon_1 / degree         =   " << mLon1 << std::endl;
    ss << "  Depth_2 / km           =   " << mD2 / 1e3 << std::endl;
    ss << "  Lat_2 / degree         =   " << mLat2 << std::endl;
    ss << "  Lon_2 / degree         =   " << mLon2 << std::endl;
    ss << "  Cylinder Radius / km   =   " << mRadius / 1e3 << std::endl;
    ss << "  HWHM lateral / km      =   " << mHWHM_lateral / 1e3 << std::endl;
    ss << "  HWHM top & bot / km    =   " << mHWHM_top_bot / 1e3 << std::endl;
    ss << "  Maximum on Axis        =   " << mMaxAxis << std::endl;
    ss << "  Reference Type         =   " << ReferenceTypesString[mReferenceType] << std::endl;
    ss << "  Affect VP              =   " << (mChangeVp ? "YES" : "NO") << std::endl;
    ss << "  Affect VS              =   " << (mChangeVs ? "YES" : "NO") << std::endl;
    ss << "  Affect Density         =   " << (mChangeRho ? "YES" : "NO") << std::endl;
    ss << "======================= 3D Volumetric ======================\n" << std::endl;
    return ss.str();
}

