// Geometric3D_Internal.cpp
// created by Kuangdai on 19-Jan-2017 
// topography on a general internal boundary

#include "Geometric3D_Internal.h"
#include <sstream>
#include <fstream>
#include "XMPI.h"
#include "XMath.h"
#include "Parameters.h"

void Geometric3D_Internal::initialize() {
    // dimension
    int nlon = mN360;
    int nlat = mN360 / 2 + 1;
    
    // read raw data
    mData = RDMatXX::Zero(nlat, nlon);
    if (XMPI::root()) {
        std::string fname = Parameters::sInputDirectory + "/" + mFileName;
        std::fstream fs(fname, std::fstream::in);
        if (!fs) throw std::runtime_error("Geometric3D_Internal::initialize || "
            "Error opening internal topography data file: ||" + fname);
        for (int i = 0; i < nlat; i++)
            for (int j = 0; j < nlon; j++) 
                fs >> mData(i, j);    
        fs.close();
    }
    mData *= 1e3;
    
    // broadcast
    XMPI::bcastEigen(mData);
    
    // Gaussian smoothing
    IColX orderRow = IColX::Constant(nlat, mGaussianOrder);
    IColX orderCol = IColX::Constant(nlon, mGaussianOrder);
    RDColX devRow = RDColX::Constant(nlat, mGaussianDev);
    RDColX devCol = RDColX::Constant(nlon, mGaussianDev);
    XMath::gaussianSmoothing(mData, orderRow, devRow, true, orderCol, devCol, false);
    
    // apply factor
    mData *= mFactor;

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

void Geometric3D_Internal::initialize(const std::vector<std::string> &params) {
    if (params.size() < 5) throw std::runtime_error("Geometric3D_Internal::initialize || "
        "Not enough parameters to initialize a Geometric3D_Internal object.");
    
    const std::string source = "Geometric3D_Internal::initialize";
        
    // initialize location 
    XMath::castValue(mRLayer, params[0], source); mRLayer *= 1e3;
    XMath::castValue(mRLower, params[1], source); mRLower *= 1e3;
    XMath::castValue(mRUpper, params[2], source); mRUpper *= 1e3;
    
    // initialize data
    XMath::castValue(mN360, params[3], source); 
    XMath::castValue(mFileName, params[4], source);
    
    try {
        int ipar = 5;
        XMath::castValue(mGeographic, params.at(ipar++), source);
        XMath::castValue(mFactor, params.at(ipar++), source);
        XMath::castValue(mGaussianOrder, params.at(ipar++), source);
        XMath::castValue(mGaussianDev, params.at(ipar++), source);
    } catch (std::out_of_range) {
        // nothing
    }
    initialize();
}

double Geometric3D_Internal::getDeltaR(double r, double theta, double phi, double rElemCenter) const {
    if (rElemCenter > mRUpper || rElemCenter < mRLower) { 
        return 0.;
    }
    
    // to geocentric
    if (mGeographic) theta = pi / 2. - XMath::theta2Lat(theta, 0.) * degree;
    
    // regularise
    if (theta < tinyDouble) theta = tinyDouble;
    if (theta > pi - tinyDouble) theta = pi - tinyDouble;
    if (phi < tinyDouble) phi = tinyDouble;
    if (phi > 2. * pi - tinyDouble) phi = 2. * pi - tinyDouble;
    
    // interpolation on sphere
    double resolution = 2. * pi / mN360;
    
    int ilat0 = (int)(theta / resolution);
    int ilat1 = ilat0 + 1;
    double flat0 = 1. - 1. / resolution * (theta - ilat0 * resolution);
    double flat1 = 1. - flat0;
    
    int ilon0 = (int)(phi / resolution);
    int ilon1 = ilon0 + 1;
    if (ilon1 == mN360) ilon1 = 0;
    double flon0 = 1. - 1. / resolution * (phi - ilon0 * resolution);
    double flon1 = 1. - flon0;
    
    double dr = 0.;
    dr += mData(ilat0, ilon0) * flat0 * flon0;
    dr += mData(ilat0, ilon1) * flat0 * flon1;
    dr += mData(ilat1, ilon0) * flat1 * flon0;
    dr += mData(ilat1, ilon1) * flat1 * flon1;
    
    // interpolation along radius    
    if (rElemCenter < mRLayer)
        return dr / (mRLayer - mRLower) * (r - mRLower);
    else 
        return dr / (mRUpper - mRLayer) * (mRUpper - r);
}

std::string Geometric3D_Internal::verbose() const {
    std::stringstream ss;
    ss << "\n======================= 3D Geometric =======================" << std::endl;
    ss << "  Model Name            =   Internal" << std::endl;
    ss << "  Depth (km)            =   " << mRLayer / 1e3 << std::endl;
    ss << "  Scope (km)            =   [" << mRLower / 1e3 << ", " << mRUpper / 1e3 << "]" << std::endl;
    ss << "  Range (km)            =   [" << mData.minCoeff() / 1e3 << ", " << mData.maxCoeff() / 1e3 << "]" << std::endl;
    ss << "  Points (Lat x Lon)    =   " << mN360 / 2 + 1 << " x " << mN360 << std::endl;
    ss << "  Data File             =   " << mFileName << std::endl;
    ss << "  Use Geographic        =   " << (mGeographic ? "YES" : "NO") << std::endl;
    ss << "  Factor                =   " << mFactor << std::endl;
    ss << "  Smoothing Order       =   " << mGaussianOrder << std::endl;
    ss << "  Smoothing Intensity   =   " << mGaussianDev << std::endl;
    ss << "======================= 3D Geometric =======================\n" << std::endl;
    return ss.str();
}



