// Volumetric3D_cylinder.cpp
// created by Kuangdai on 19-Oct-2016 
// a cylinder-shaped heterogeneity

#include "Volumetric3D_cylinder.h"
#include "Parameters.h"
#include "Geodesy.h"
#include <boost/algorithm/string.hpp>
#include <sstream>

void Volumetric3D_cylinder::initialize(const std::vector<std::string> &params) {
    // need at least 10 parameters to make a cylinder
    if (params.size() < 10) {
        throw std::runtime_error("Volumetric3D_cylinder::initialize || "
            "Not enough parameters for a cylinder-shaped heterogeneity. Need 10 at least.");
    }
        
    const std::string source = "Volumetric3D_cylinder::initialize";
    
    // property name
    bool found = false;
    for (int i = 0; i < Volumetric3D::MaterialPropertyString.size(); i++) {
        if (boost::iequals(params[0], Volumetric3D::MaterialPropertyString[i])) {
            mMaterialProp = Volumetric3D::MaterialProperty(i);
            found = true;
            break;
        }
    }
    if (!found) {
        throw std::runtime_error("Volumetric3D_cylinder::initialize || "
            "Unknown material property, name = " + params[0]);
    }
    
    // reference type
    found = false;
    for (int i = 0; i < Volumetric3D::MaterialRefTypeString.size(); i++) {
        if (boost::iequals(params[1], Volumetric3D::MaterialRefTypeString[i]) ||
            boost::iequals(params[1], Volumetric3D::MaterialRefTypeStringShort[i])) {
            mReferenceType = Volumetric3D::MaterialRefType(i);
            found = true;
            break;
        }
    }
    if (!found) {
        throw std::runtime_error("Volumetric3D_cylinder::initialize || "
            "Unknown material reference type, type = " + params[1]);
    }
    
    // value inside
    Parameters::castValue(mValueInside, params[2], source);
    
    // radius
    Parameters::castValue(mRadius, params[3], source); mRadius *= 1e3;
    
    // location    
    Parameters::castValue(mD1, params[4], source); mD1 *= 1e3;
    Parameters::castValue(mLat1, params[5], source);
    Parameters::castValue(mLon1, params[6], source);
    
    Parameters::castValue(mD2, params[7], source); mD2 *= 1e3;
    Parameters::castValue(mLat2, params[8], source);
    Parameters::castValue(mLon2, params[9], source);
    
    // optional
    try {
        int ipar = 10;
        Parameters::castValue(mSourceCentered, params.at(ipar++), source);
        Parameters::castValue(mFluid, params.at(ipar++), source);
        Parameters::castValue(mHWHM_lateral, params.at(ipar++), source); mHWHM_lateral *= 1e3;
        Parameters::castValue(mHWHM_top_bot, params.at(ipar++), source); mHWHM_top_bot *= 1e3;
    } catch (std::out_of_range) {
        // nothing
    }    
    
    // compute xyz of endpoints and length
    RDCol3 rtpPoint1, rtpPoint2;
    if (mSourceCentered) {
        RDCol3 rtpPoint1Src, rtpPoint2Src;
        rtpPoint1Src(0) = Geodesy::getROuter() - mD1;
        rtpPoint1Src(1) = mLat1 * degree;
        rtpPoint1Src(2) = mLon1 * degree;
        rtpPoint1 = Geodesy::rotateSrc2Glob(rtpPoint1Src, mSrcLat, mSrcLon, mSrcDep);
        rtpPoint2Src(0) = Geodesy::getROuter() - mD2;
        rtpPoint2Src(1) = mLat2 * degree;
        rtpPoint2Src(2) = mLon2 * degree;
        rtpPoint2 = Geodesy::rotateSrc2Glob(rtpPoint2Src, mSrcLat, mSrcLon, mSrcDep);
    } else {
        rtpPoint1(0) = Geodesy::getROuter() - mD1;
        rtpPoint1(1) = Geodesy::lat2Theta_d(mLat1, mD1);
        rtpPoint1(2) = Geodesy::lon2Phi(mLon1);
        rtpPoint2(0) = Geodesy::getROuter() - mD2;
        rtpPoint2(1) = Geodesy::lat2Theta_d(mLat2, mD2);
        rtpPoint2(2) = Geodesy::lon2Phi(mLon2);
    }
    mXyzPoint1 = Geodesy::toCartesian(rtpPoint1);
    mXyzPoint2 = Geodesy::toCartesian(rtpPoint2);
    mLength = (mXyzPoint1 - mXyzPoint2).norm();
    
    // use 20% of radius for lateral HWHM if not specified
    if (mHWHM_lateral < 0.) {
        mHWHM_lateral = mRadius * .2;
    }
    
    // use 10% of cylinder length for top-bot HWHM if not specified
    if (mHWHM_top_bot < 0.) {
        mHWHM_top_bot = mLength * .1;
    }
    
    // for Absolute models
    if (mReferenceType == Volumetric3D::MaterialRefType::Absolute) {
        // decay is not allowed
        mHWHM_lateral = mHWHM_top_bot = 0.;
        // convert to SI
        mValueInside *= MaterialPropertyAbsSI[mMaterialProp];
    }
}

bool Volumetric3D_cylinder::get3dProperties(double r, double theta, double phi, double rElemCenter,
    std::vector<MaterialProperty> &properties, 
    std::vector<MaterialRefType> &refTypes,
    std::vector<double> &values) const {
    
    // header
    properties = std::vector<MaterialProperty>(1, mMaterialProp);
    refTypes = std::vector<MaterialRefType>(1, mReferenceType);
    values = std::vector<double>(1, 0.);
    
    // distance from point to axis
    RDCol3 rtpTarget;
    rtpTarget(0) = r;
    rtpTarget(1) = theta;
    rtpTarget(2) = phi;
    const RDCol3 &xyzTarget = Geodesy::toCartesian(rtpTarget);
    const RDCol3 &xyzDiff1 = xyzTarget - mXyzPoint1;
    const RDCol3 &xyzDiff2 = xyzTarget - mXyzPoint2; 
    double distToLine = xyzDiff1.cross(xyzDiff2).norm() / mLength;
    
    // outside range
    if (distToLine > mRadius + 4. * mHWHM_lateral) {
        return false;
    }
    
    // compute Gaussian lateral
    double distToSurf = distToLine - mRadius;
    if (distToSurf < 0.) {
        distToSurf = 0.;
    }
    double stddev = mHWHM_lateral / sqrt(2. * log(2.));
    double gaussian_lateral = mValueInside * exp(-distToSurf * distToSurf / (stddev * stddev * 2.));
    
    // distance from point to ends
    double d1 = xyzDiff1.norm();
    double d2 = xyzDiff2.norm();
    double dmax = std::max(d1, d2);
    double distToTopBot = sqrt(dmax * dmax - distToLine * distToLine) - mLength;
    
    // outside range
    if (distToTopBot > 4. * mHWHM_top_bot) {
        return false;
    }
    
    double gaussian_topbot = gaussian_lateral;
    if (distToTopBot > 0.) {
        // compute Gaussian top-bottom
        stddev = mHWHM_top_bot / sqrt(2. * log(2.));
        gaussian_topbot = gaussian_lateral * exp(-gaussian_topbot * gaussian_topbot / (stddev * stddev * 2.)); 
    }
    
    // set perturbations    
    values[0] = gaussian_topbot;
    return true;    
}

std::string Volumetric3D_cylinder::verbose() const {
    std::stringstream ss;
    ss << "\n======================= 3D Volumetric ======================" << std::endl;
    ss << "  Model Name               =   cylinder" << std::endl;
    ss << "  Material Property        =   " << MaterialPropertyString[mMaterialProp] << std::endl;
    ss << "  Reference Type           =   " << MaterialRefTypeString[mReferenceType] << std::endl;
    if (mReferenceType == Volumetric3D::MaterialRefType::Absolute) {
        ss << "  Value Inside             =   " << mValueInside / MaterialPropertyAbsSI[mMaterialProp] << std::endl;
    } else {
        ss << "  Value Inside             =   " << mValueInside << std::endl;
    }
    ss << "  Cylinder Radius / km     =   " << mRadius / 1e3 << std::endl;
    ss << "  Depth_1 / km             =   " << mD1 / 1e3 << std::endl;
    ss << "  Lat_1 or Theta_1 / deg   =   " << mLat1 << std::endl;
    ss << "  Lon_1 or Phi_1 / deg     =   " << mLon1 << std::endl;
    ss << "  Depth_2 / km             =   " << mD2 / 1e3 << std::endl;
    ss << "  Lat_2 or Theta_2 / deg   =   " << mLat2 << std::endl;
    ss << "  Lon_2 or Phi_2 / deg     =   " << mLon2 << std::endl;
    ss << "  Source-centered          =   " << (mSourceCentered ? "YES" : "NO") << std::endl;
    ss << "  HWHM lateral / km        =   " << mHWHM_lateral / 1e3 << std::endl;
    ss << "  HWHM top-bot / km        =   " << mHWHM_top_bot / 1e3 << std::endl;
    ss << "======================= 3D Volumetric ======================\n" << std::endl;
    return ss.str();
}

