// Volumetric3D_EMC.cpp
// created by Kuangdai on 16-May-2017 
// genetral Volumetric3D model with IRIS-EMC format

#include "Volumetric3D_EMC.h"
#include <sstream>
#include "XMPI.h"
#include "XMath.h"
#include "Geodesy.h"
#include "Parameters.h"
#include "NetCDF_Reader.h"
#include "NetCDF_ReaderAscii.h"

#include <dirent.h>
#include <sys/types.h>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <algorithm>

bool compareFunc(const std::pair<std::string, float> &a, 
    const std::pair<std::string, float> &b) {
    return a.second < b.second;
}

void Volumetric3D_EMC::initialize() {
    // meta data
    Eigen::Matrix<float, Eigen::Dynamic, 1> fdata, fdep, flat, flon;
    if (XMPI::root()) {
        std::vector<size_t> dims;
        std::string fname = Parameters::sInputDirectory + "/" + mFileName;
        if (mOneFilePerDepth) {
            // get all files
            DIR *dir = opendir(fname.c_str());
            if (!dir) {
                throw std::runtime_error("Volumetric3D_EMC::initialize || " 
                    "Error opening directory of data files || Directory = " + fname);    
            }
            struct dirent *entry;
            std::vector<std::pair<std::string, float>> fileDepth;
            while ((entry = readdir(dir)) != NULL) {
                // file name
                std::string fn(entry->d_name);
                if (fn.find(".grd") == std::string::npos) {
                    continue;
                }
                // depth
                std::string depthStr(fn);
                boost::replace_first(depthStr, ".grd", "");
                std::size_t found = depthStr.find_last_of("_");
                if (found == std::string::npos) {
                    throw std::runtime_error("Volumetric3D_EMC::initialize || " 
                        "Error processing file name || File = " + fn);    
                }
                depthStr = depthStr.substr(found + 1);
                float depth = 0.;
                try {
                    depth = boost::lexical_cast<float>(depthStr);
                } catch (std::exception) {
                    throw std::runtime_error("Volumetric3D_EMC::initialize || " 
                        "Error processing file name || File = " + fn);
                }
                // add pair
                fileDepth.push_back(std::make_pair(fname + "/" + fn, depth));
            }
            closedir(dir);
            
            // sort depth
            std::sort(fileDepth.begin(), fileDepth.end(), compareFunc);
            
            // read lat and lon
            Eigen::Matrix<double, Eigen::Dynamic, 1> dlat, dlon;
            if (NetCDF_Reader::checkNetCDF_isAscii(fileDepth[0].first)) {
                NetCDF_ReaderAscii reader;
                reader.open(fileDepth[0].first);
                reader.read1D("lat", dlat);
                reader.read1D("lon", dlon);
                reader.close();
            } else {
                NetCDF_Reader reader;
                reader.open(fileDepth[0].first);
                reader.read1D("lat", dlat);
                reader.read1D("lon", dlon);
                reader.close();
            }
            flat = dlat.cast<float>();
            flon = dlon.cast<float>();
            
            // depths and data
            size_t depthLen = flat.size() * flon.size();   
            size_t totalLen = depthLen * fileDepth.size();
            fdep.resize(fileDepth.size());
            fdata.resize(totalLen);
            for (int i = 0; i < fileDepth.size(); i++) {
                fdep(i) = fileDepth[i].second;
                Eigen::Matrix<float, Eigen::Dynamic, 1> temp;
                if (NetCDF_Reader::checkNetCDF_isAscii(fileDepth[i].first)) {
                    NetCDF_ReaderAscii reader;
                    reader.open(fileDepth[i].first);
                    reader.readMetaData(mVarName, temp, dims);
                    reader.close();
                } else {
                    NetCDF_Reader reader;
                    reader.open(fileDepth[i].first);
                    reader.readMetaData(mVarName, temp, dims);
                    reader.close();
                }
                fdata.block(i * depthLen, 0, depthLen, 1) = temp;
            }
            dims.insert(dims.begin(), fdep.size());
            
            // perturbations are given in percentage
            fdata *= .01f;
        } else {
            if (NetCDF_Reader::checkNetCDF_isAscii(fname)) {
                NetCDF_ReaderAscii reader;
                reader.open(fname);
                reader.read1D("depth", fdep);
                reader.read1D("latitude", flat);
                reader.read1D("longitude", flon);
                reader.readMetaData(mVarName, fdata, dims);
                reader.close();
            } else {
                NetCDF_Reader reader;
                reader.open(fname);
                reader.read1D("depth", fdep);
                reader.read1D("latitude", flat);
                reader.read1D("longitude", flon);
                reader.readMetaData(mVarName, fdata, dims);
                reader.close();
            }
        }
        
        // check dimensions
        if (dims.size() != 3) {
            throw std::runtime_error("Volumetric3D_EMC::initialize || Inconsistent data dimensions || "
                "File/Directory = " + fname);
        }
        if (dims[0] != fdep.size() || dims[1] != flat.size() || dims[2] != flon.size()) {
            throw std::runtime_error("Volumetric3D_EMC::initialize || Inconsistent data dimensions || "
                "File/Directory = " + fname);
        }
        if (!XMath::sortedAscending(fdep) || !XMath::sortedAscending(flat) || !XMath::sortedAscending(flon)) {
            throw std::runtime_error("Volumetric3D_EMC::initialize || Grid coordinates are not sorted ascendingly || "
                "File/Directory = " + fname);
        }
    }
    XMPI::bcastEigen(fdep);
    XMPI::bcastEigen(flat);
    XMPI::bcastEigen(flon);
    XMPI::bcastEigen(fdata);
    
    mGridDep = fdep.cast<double>();
    mGridLat = flat.cast<double>();
    mGridLon = flon.cast<double>();
    RDColX data = fdata.cast<double>();
    
    // SI
    mGridDep *= 1e3;
    if (mReferenceType == Volumetric3D::MaterialRefType::Absolute) {
        // convert to SI
        data *= MaterialPropertyAbsSI[mMaterialProp];
    }
    
    // apply factor
    data *= mFactor;
    
    // special flag
    if (!boost::iequals(mModelFlag, "none")) {
        if (boost::iequals(mModelFlag, "abs")) {
            data.array() = data.array().abs();
            data *= mModelFlagFactor;
        } else if (boost::iequals(mModelFlag, "pow")) {
            double absmax = data.array().abs().maxCoeff();
            double scalefact = std::abs(absmax / std::pow(absmax, mModelFlagFactor));
            data.array() = data.array().sign() * (data.array().abs().pow(mModelFlagFactor) * scalefact);
        } else {
            throw std::runtime_error("Volumetric3D_EMC::initialize || "
                "Unknown special model flag, flag = " + mModelFlag);
        }
    }
    
    // reshape data
    int pos = 0;
    for (int i = 0; i < mGridDep.size(); i++) {
        RDMatXX mat(mGridLat.size(), mGridLon.size());
        for (int j = 0; j < mGridLat.size(); j++) {
            for (int k = 0; k < mGridLon.size(); k++) {
                mat(j, k) = data(pos++);
            }
        }
        mGridData.push_back(mat);
    }
}

void Volumetric3D_EMC::initialize(const std::vector<std::string> &params) {
    if (params.size() < 4) throw std::runtime_error("Volumetric3D_EMC::initialize || "
        "Not enough parameters to initialize a Volumetric3D_EMC object, at least 4 needed.");
    
    const std::string source = "Volumetric3D_EMC::initialize";
        
    // initialize data
    Parameters::castValue(mFileName, params[0], source);
    Parameters::castValue(mVarName, params[1], source);
    
    // property name
    bool found = false;
    for (int i = 0; i < Volumetric3D::MaterialPropertyString.size(); i++) {
        if (boost::iequals(params[2], Volumetric3D::MaterialPropertyString[i])) {
            mMaterialProp = Volumetric3D::MaterialProperty(i);
            found = true;
            break;
        }
    }
    if (!found) {
        throw std::runtime_error("Volumetric3D_EMC::initialize || "
            "Unknown material property, name = " + params[2]);
    }
    
    // reference type
    found = false;
    for (int i = 0; i < Volumetric3D::MaterialRefTypeString.size(); i++) {
        if (boost::iequals(params[3], Volumetric3D::MaterialRefTypeString[i]) ||
            boost::iequals(params[3], Volumetric3D::MaterialRefTypeStringShort[i])) {
            mReferenceType = Volumetric3D::MaterialRefType(i);
            found = true;
            break;
        }
    }
    if (!found) {
        throw std::runtime_error("Volumetric3D_EMC::initialize || "
            "Unknown material reference type, type = " + params[3]);
    }
    
    try {
        int ipar = 4;
        Parameters::castValue(mFactor, params.at(ipar++), source);
        Parameters::castValue(mGeographic, params.at(ipar++), source);
        Parameters::castValue(mOneFilePerDepth, params.at(ipar++), source);
        Parameters::castValue(mModelFlag, params.at(ipar++), source);
        Parameters::castValue(mModelFlagFactor, params.at(ipar++), source);
    } catch (std::out_of_range) {
        // nothing
    }
    
    if (!boost::iequals(mModelFlag, "none")) {
        if (mReferenceType == MaterialRefType::Absolute) {
            throw std::runtime_error("Volumetric3D_EMC::initialize || "
                "Imposing special model flag on an absolute model.");
        }
    }
    initialize();
}

bool Volumetric3D_EMC::get3dProperties(double r, double theta, double phi, double rElemCenter,
    std::vector<MaterialProperty> &properties, 
    std::vector<MaterialRefType> &refTypes,
    std::vector<double> &values) const {
    
    // header
    properties = std::vector<MaterialProperty>(1, mMaterialProp);
    refTypes = std::vector<MaterialRefType>(1, mReferenceType);
    values = std::vector<double>(1, 0.);
        
    // to geocentric
    if (mGeographic) {
        // which radius to use?
        theta = pi / 2. - Geodesy::theta2Lat_d(theta, 0.) * degree;
    }
    
    // regularise
    double dep = Geodesy::getROuter() - r;
    double lat = 90. - theta / degree;
    double lon = phi / degree;
    XMath::checkLimits(dep, 0., Geodesy::getROuter());
    XMath::checkLimits(lat, -90., 90.);
    if (mGridLon[0] < 0.) {
        // lon starts from -180.
        if (lon > 180.) {
            lon -= 360.;
        }
        XMath::checkLimits(lon, -180., 180.);
    } else {
        // lon starts from 0.
        XMath::checkLimits(lon, 0., 360.);
    }
    
    // check center
    double dmin = mGridDep[0];
    double dmax = mGridDep[mGridDep.size() - 1];
    double dcenter = Geodesy::getROuter() - rElemCenter;
    if (dcenter < dmin || dcenter > dmax) {
        return false;
    }
    if (dep < dmin && dep > dmin * 0.999999) {
        dep = dmin;
    }
    if (dep > dmax && dep < dmax * 1.000001) {
        dep = dmax;
    }
    
    // interpolation
    int ldep0, llat0, llon0, ldep1, llat1, llon1;
    double wdep0, wlat0, wlon0, wdep1, wlat1, wlon1;
    XMath::interpLinear(dep, mGridDep, ldep0, wdep0);
    XMath::interpLinear(lat, mGridLat, llat0, wlat0);
    XMath::interpLinear(lon, mGridLon, llon0, wlon0);    
    if (ldep0 < 0 || llat0 < 0 || llon0 < 0) {
        return false;
    }
    
    ldep1 = ldep0 + 1;
    llat1 = llat0 + 1;
    llon1 = llon0 + 1;
    wdep1 = 1. - wdep0;
    wlat1 = 1. - wlat0;
    wlon1 = 1. - wlon0;
    
    values[0] += mGridData[ldep0](llat0, llon0) * wdep0 * wlat0 * wlon0;
    values[0] += mGridData[ldep0](llat1, llon0) * wdep0 * wlat1 * wlon0;
    values[0] += mGridData[ldep0](llat0, llon1) * wdep0 * wlat0 * wlon1;
    values[0] += mGridData[ldep0](llat1, llon1) * wdep0 * wlat1 * wlon1;
    values[0] += mGridData[ldep1](llat0, llon0) * wdep1 * wlat0 * wlon0;
    values[0] += mGridData[ldep1](llat1, llon0) * wdep1 * wlat1 * wlon0;
    values[0] += mGridData[ldep1](llat0, llon1) * wdep1 * wlat0 * wlon1;
    values[0] += mGridData[ldep1](llat1, llon1) * wdep1 * wlat1 * wlon1;
    return true;
}

std::string Volumetric3D_EMC::verbose() const {
    std::stringstream ss;
    ss << "\n======================= 3D Volumetric =======================" << std::endl;
    ss << "  Model Name           =   EMC" << std::endl;
    ss << "  Data File            =   " << mFileName << std::endl;
    ss << "  Variable Name        =   " << mVarName << std::endl;
    ss << "  Material Property    =   " << MaterialPropertyString[mMaterialProp] << std::endl;
    ss << "  Reference Type       =   " << MaterialRefTypeString[mReferenceType] << std::endl;
    ss << "  Num. Depths          =   " << mGridDep.size() << std::endl;
    ss << "  Num. Latitudes       =   " << mGridLat.size() << std::endl;
    ss << "  Num. Longitudes      =   " << mGridLon.size() << std::endl;
    ss << "  Depth Range          =   [" << mGridDep.minCoeff() << ", " << mGridDep.maxCoeff() << "]" << std::endl;
    ss << "  Latitude Range       =   [" << mGridLat.minCoeff() << ", " << mGridLat.maxCoeff() << "]" << std::endl;
    ss << "  Longitude Range      =   [" << mGridLon.minCoeff() << ", " << mGridLon.maxCoeff() << "]" << std::endl;
    ss << "  Factor               =   " << mFactor << std::endl;
    ss << "  Use Geographic       =   " << (mGeographic ? "YES" : "NO") << std::endl;
    ss << "  One File per Depth   =   " << (mOneFilePerDepth ? "YES" : "NO") << std::endl;
    if (!boost::iequals(mModelFlag, "none")) {
        ss << "  Special Model Flag   =   " << mModelFlag << std::endl;
        ss << "  Model Flag Factor    =   " << mModelFlagFactor << std::endl;
    }
    ss << "======================= 3D Volumetric =======================\n" << std::endl;
    return ss.str();
}



