// Geometric3D_EMC.cpp
// created by Kuangdai on 16-May-2017 
// topography on a boundary at any depth, with IRIS-EMC format

#include "Geometric3D_EMC.h"
#include <sstream>
#include "XMPI.h"
#include "XMath.h"
#include "Geodesy.h"
#include "Parameters.h"
#include "NetCDF_Reader.h"
#include "NetCDF_ReaderAscii.h"

void Geometric3D_EMC::initialize() {
    // EMC is in float
    Eigen::Matrix<float, Eigen::Dynamic, 1> flat, flon;
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> fdata;
    
    // read file
    if (XMPI::root()) {
        std::string fname = Parameters::sInputDirectory + "/" + mFileName;
        if (NetCDF_Reader::checkNetCDF_isAscii(fname)) {
            NetCDF_ReaderAscii reader;
            reader.open(fname);
            reader.read1D("latitude", flat);
            reader.read1D("longitude", flon);
            reader.read2D(mVarName, fdata);
            reader.close();
        } else {
            NetCDF_Reader reader;
            reader.open(fname);
            reader.read1D("latitude", flat);
            reader.read1D("longitude", flon);
            reader.read2D(mVarName, fdata);
            reader.close();
        }
        if (fdata.rows() != flat.size() || fdata.cols() != flon.size()) {
            throw std::runtime_error("Geometric3D_EMC::initialize || "
                "Inconsistent data dimensions || File = " + fname);
        }
        if (!XMath::sortedAscending(flat) || !XMath::sortedAscending(flon)) {
            throw std::runtime_error("Geometric3D_EMC::initialize || "
                "Grid coordinates are not sorted ascendingly || File = " + fname);
        }
    }
    
    // broadcast
    XMPI::bcastEigen(flat);
    XMPI::bcastEigen(flon);
    XMPI::bcastEigen(fdata);
    
    // to double
    mGridLat = flat.cast<double>();
    mGridLon = flon.cast<double>();
    mGridData = fdata.cast<double>();
    
    // to SI
    mGridData *= 1e3;
    
    // apply factor
    mGridData *= mFactor;

    ////////// plot computed data ////////////  
    // std::fstream fsdr;
    // fsdr.open("/Users/kuangdai/Desktop/crust1/cmb.txt", std::fstream::out);
    // double r = mRLayer; // double r = mRSurf;
    // int intGrid = 1;
    // for (int i = 0; i <= 180 * intGrid; i++) {
    //     double theta = i * degree / intGrid;
    //     for (int j = -180 * intGrid; j <= 180 * intGrid; j++) {
    //         double phi = j * degree / intGrid;
    //         if (phi < 0) phi += 2. * pi;
    //         fsdr << getDeltaR(r, theta, phi, r) << " ";
    //     }
    //     fsdr << std::endl;
    // }   
    // fsdr.close();
    // exit(0);
    ////////// plot computed data ////////////  
}

void Geometric3D_EMC::initialize(const std::vector<std::string> &params) {
    if (params.size() < 5) throw std::runtime_error("Geometric3D_EMC::initialize || "
        "Not enough parameters to initialize a Geometric3D_EMC object, at least 5 needed.");
    
    const std::string source = "Geometric3D_EMC::initialize";
        
    // initialize location 
    Parameters::castValue(mRLayer, params[0], source); mRLayer *= 1e3;
    Parameters::castValue(mRLower, params[1], source); mRLower *= 1e3;
    Parameters::castValue(mRUpper, params[2], source); mRUpper *= 1e3;
    
    // initialize data
    Parameters::castValue(mFileName, params[3], source);
    Parameters::castValue(mVarName, params[4], source);
    
    try {
        int ipar = 5;
        Parameters::castValue(mFactor, params.at(ipar++), source);
        Parameters::castValue(mGeographic, params.at(ipar++), source);
    } catch (std::out_of_range) {
        // nothing
    }
    initialize();
}

double Geometric3D_EMC::getDeltaR(double r, double theta, double phi, double rElemCenter) const {
    if (r > mRUpper || r < mRLower) { 
        return 0.;
    }
    
    // to geocentric
    if (mGeographic) {
        // which radius to use?
        theta = pi / 2. - Geodesy::theta2Lat_d(theta, 0.) * degree;
    }
    
    // regularise
    double lat = 90. - theta / degree;
    double lon = phi / degree;
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
    
    // interpolation on sphere
    int llat0, llon0, llat1, llon1;
    double wlat0, wlon0, wlat1, wlon1;
    XMath::interpLinear(lat, mGridLat, llat0, wlat0);
    XMath::interpLinear(lon, mGridLon, llon0, wlon0);    
    if (llat0 < 0 || llon0 < 0) {
        return 0.;
    }
    
    llat1 = llat0 + 1;
    llon1 = llon0 + 1;
    wlat1 = 1. - wlat0;
    wlon1 = 1. - wlon0;
    
    double dr = 0.;
    dr += mGridData(llat0, llon0) * wlat0 * wlon0;
    dr += mGridData(llat1, llon0) * wlat1 * wlon0;
    dr += mGridData(llat0, llon1) * wlat0 * wlon1;
    dr += mGridData(llat1, llon1) * wlat1 * wlon1;
    
    // interpolation along radius    
    if (rElemCenter < mRLayer) {
        return dr / (mRLayer - mRLower) * (r - mRLower);
    } else {
        return dr / (mRUpper - mRLayer) * (mRUpper - r);
    }
}

std::string Geometric3D_EMC::verbose() const {
    std::stringstream ss;
    ss << "\n======================= 3D Geometric =======================" << std::endl;
    ss << "  Model Name          =   EMC" << std::endl;
    ss << "  Layer radius (km)   =   " << mRLayer / 1e3 << std::endl;
    ss << "  Radius range (km)   =   [" << mRLower / 1e3 << ", " << mRUpper / 1e3 << "]" << std::endl;
    ss << "  Latitude Range      =   [" << mGridLat.minCoeff() << ", " << mGridLat.maxCoeff() << "]" << std::endl;
    ss << "  Longitude Range     =   [" << mGridLon.minCoeff() << ", " << mGridLon.maxCoeff() << "]" << std::endl;
    ss << "  Data Range (km)     =   [" << mGridData.minCoeff() / 1e3 << ", " << mGridData.maxCoeff() / 1e3 << "]" << std::endl;
    ss << "  Data File           =   " << mFileName << std::endl;
    ss << "  Variable Name       =   " << mVarName << std::endl;
    ss << "  Num. Latitudes      =   " << mGridLat.size() << std::endl;
    ss << "  Num. Longitudes     =   " << mGridLon.size() << std::endl;
    ss << "  Factor              =   " << mFactor << std::endl;
    ss << "  Use Geographic      =   " << (mGeographic ? "YES" : "NO") << std::endl;
    ss << "======================= 3D Geometric =======================\n" << std::endl;
    return ss.str();
}



