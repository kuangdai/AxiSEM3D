// Volumetric3D_bubble.cpp
// created by Kuangdai on 19-Oct-2016 
// a bubble-shaped heterogeneity

#include "Volumetric3D_bubble.h"
#include "Parameters.h"
#include "Geodesy.h"
#include <boost/algorithm/string.hpp>
#include <sstream>

void Volumetric3D_bubble::initialize(const std::vector<std::string> &params) {
    // need at least 7 parameters to make a bubble
    if (params.size() < 7) {
        throw std::runtime_error("Volumetric3D_bubble::initialize || "
            "Not enough parameters for a bubble-shaped heterogeneity. Need 7 at least.");
    }
        
    const std::string source = "Volumetric3D_bubble::initialize";
    
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
        throw std::runtime_error("Volumetric3D_bubble::initialize || "
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
        throw std::runtime_error("Volumetric3D_bubble::initialize || "
            "Unknown material reference type, type = " + params[1]);
    }
    
    // value inside
    Parameters::castValue(mValueInside, params[2], source);
    
    // radius
    Parameters::castValue(mRadius, params[3], source); mRadius *= 1e3;
    
    // location    
    Parameters::castValue(mDepth, params[4], source); mDepth *= 1e3;
    Parameters::castValue(mLat, params[5], source);
    Parameters::castValue(mLon, params[6], source);
    
    // optional
    try {
        int ipar = 7;
        Parameters::castValue(mSourceCentered, params.at(ipar++), source);
        Parameters::castValue(mFluid, params.at(ipar++), source);
        Parameters::castValue(mHWHM, params.at(ipar++), source); mHWHM *= 1e3;
    } catch (std::out_of_range) {
        // nothing
    }    
    
    // compute xyz of endpoints and length
    RDCol3 rtpBubble;
    if (mSourceCentered) {
        RDCol3 rtpBubbleSrc;
        rtpBubbleSrc(0) = Geodesy::getROuter() - mDepth;
        rtpBubbleSrc(1) = mLat * degree;
        rtpBubbleSrc(2) = mLon * degree;
        rtpBubble = Geodesy::rotateSrc2Glob(rtpBubbleSrc, mSrcLat, mSrcLon, mSrcDep);
    } else {
        rtpBubble(0) = Geodesy::getROuter() - mDepth;
        rtpBubble(1) = Geodesy::lat2Theta_d(mLat, mDepth);
        rtpBubble(2) = Geodesy::lon2Phi(mLon);
    }
    mXyzBubble = Geodesy::toCartesian(rtpBubble);
    
    // use 20% of radius for HWHM if not specified
    if (mHWHM < 0.) {
        mHWHM = mRadius * .2;
    }
    
    // for Absolute models
    if (mReferenceType == Volumetric3D::MaterialRefType::Absolute) {
        // decay is not allowed
        mHWHM = 0.;
        // convert to SI
        mValueInside *= MaterialPropertyAbsSI[mMaterialProp];
    }
}

bool Volumetric3D_bubble::get3dProperties(double r, double theta, double phi, double rElemCenter,
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
    double distance = (mXyzBubble - xyzTarget).norm();
    
    // treat as center if inside bubble
    distance -= mRadius; 
    if (distance < 0.) {
        distance = 0.;
    }
    
    // outside range
    if (distance > 4. * mHWHM) {
        return false;
    }
    
    // compute Gaussian
    double stddev = mHWHM / sqrt(2. * log(2.));
    double gaussian = mValueInside * exp(-distance * distance / (stddev * stddev * 2.));
    
    // set perturbations    
    values[0] = gaussian;
    return true;    
}

std::string Volumetric3D_bubble::verbose() const {
    std::stringstream ss;
    ss << "\n======================= 3D Volumetric ======================" << std::endl;
    ss << "  Model Name           =   bubble" << std::endl;
    ss << "  Material Property    =   " << MaterialPropertyString[mMaterialProp] << std::endl;
    ss << "  Reference Type       =   " << MaterialRefTypeString[mReferenceType] << std::endl;
    if (mReferenceType == Volumetric3D::MaterialRefType::Absolute) {
        ss << "  Value Inside         =   " << mValueInside / MaterialPropertyAbsSI[mMaterialProp] << std::endl;
    } else {
        ss << "  Value Inside         =   " << mValueInside << std::endl;
    }
    ss << "  Bubble Radius / km   =   " << mRadius / 1e3 << std::endl;
    ss << "  Depth / km           =   " << mDepth / 1e3 << std::endl;
    ss << "  Lat or Theta / deg   =   " << mLat << std::endl;
    ss << "  Lon or Phi / deg     =   " << mLon << std::endl;
    ss << "  Source-centered      =   " << (mSourceCentered ? "YES" : "NO") << std::endl;
    ss << "  HWHM / km            =   " << mHWHM / 1e3 << std::endl;
    ss << "======================= 3D Volumetric ======================\n" << std::endl;
    return ss.str();
}

