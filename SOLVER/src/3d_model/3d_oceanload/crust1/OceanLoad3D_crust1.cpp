// OceanLoad3D_crust1.h
// created by Kuangdai on 8-Oct-2016 
// crustal model CRUST 1.0 
// http://igppweb.ucsd.edu/~gabi/crust1.html

#include "OceanLoad3D_crust1.h"
#include "Volumetric3D_crust1.h"
#include <sstream>
#include <fstream>
#include "XMPI.h"
#include "XMath.h"

const size_t OceanLoad3D_crust1::sNLayer = 9;
const size_t OceanLoad3D_crust1::sNLat = 180;
const size_t OceanLoad3D_crust1::sNLon = 360;

void OceanLoad3D_crust1::initialize() {
    // read raw data
    size_t nrow = sNLat * sNLon;
    RDMatXX elevation = RDMatXX::Zero(nrow, sNLayer);
    if (XMPI::root()) {
        std::string fname = projectDirectory + "/src/3d_model/3d_volumetric/crust1/data/crust1.bnds";
        std::fstream fs(fname, std::fstream::in);
        if (!fs) {
            throw std::runtime_error("OceanLoad3D_crust1::initialize || "
                "Error opening crust1.0 data file: ||" + fname);
        }
        for (size_t i = 0; i < nrow; i++) {
            for (size_t j = 0; j < sNLayer; j++) {
                fs >> elevation(i, j);
            }
        }
        fs.close();
    }
    // broadcast
    XMPI::bcastEigen(elevation);
    
    // water depth
    size_t colWaterBot = 1;
    if (mIncludeIceAsWater) {
        colWaterBot = 2;
    }
    RDColX depthVec = (elevation.col(0) - elevation.col(colWaterBot)) * 1e3;
    
    // determine ocean depth ONLY by elevation for SPECFEM benchmark
    if (mBenchmarkSPECFEM) {
        depthVec = elevation.col(2) * 1e3;
        mGaussianOrder = 0;
        mNPointInterp = 2;
        mIncludeIceAsWater = false;
        mGeographic = false;
    }
    
    RDMatXX depth(sNLat, sNLon);
    for (size_t i = 0; i < sNLat; i++) {
        depth.row(i) = depthVec.block(i * sNLon, 0, sNLon, 1).transpose();
    } 
        
    //////////// plot raw data ////////////  
    // std::fstream fs;
    // fs.open("/Users/kuangdai/Desktop/crust1/water.txt", std::fstream::out);
    // fs << depth << std::endl;
    // fs.close();
    //////////// plot raw data //////////// 
    
    // Gaussian smoothing
    IColX orderRow = IColX::Constant(sNLat, mGaussianOrder);
    IColX orderCol = IColX::Constant(sNLon, mGaussianOrder);
    RDColX devRow = RDColX::Constant(sNLat, mGaussianDev);
    RDColX devCol = RDColX::Constant(sNLon, mGaussianDev);
    XMath::gaussianSmoothing(depth, orderRow, devRow, true, orderCol, devCol, false);
    
    // cast to integer theta with unique polar values
    mDepth = RDMatXX::Zero(sNLat + 1, sNLon);
    // fill north and south pole
    mDepth.row(0).fill(depth.row(0).sum() / sNLon);
    mDepth.row(sNLat).fill(depth.row(sNLat - 1).sum() / sNLon);
    // interp at integer theta
    for (size_t i = 1; i < sNLat; i++) {
        mDepth.row(i) = (depth.row(i - 1) + depth.row(i)) * .5;
    }
    
    //////////// plot raw data ////////////  
    // std::fstream fs;
    // fs.open("/Users/kuangdai/Desktop/crust1/smwater.txt", std::fstream::out);
    // fs << mDepth << std::endl;
    // fs.close();
    //////////// plot raw data //////////// 
    
}

void OceanLoad3D_crust1::initialize(const std::vector<std::string> &params) {
    try {
        size_t ipar = 0;
        const std::string source = "OceanLoad3D_crust1::initialize";
        XMath::castValue(mGaussianOrder, params.at(ipar++), source);
        XMath::castValue(mGaussianDev, params.at(ipar++), source);
        XMath::castValue(mNPointInterp, params.at(ipar++), source);
        XMath::castValue(mGeographic, params.at(ipar++), source);
        XMath::castValue(mIncludeIceAsWater, params.at(ipar++), source);
        XMath::castValue(mBenchmarkSPECFEM, params.at(ipar++), source);
    } catch (std::out_of_range) {
        // nothing
    }
    initialize();
}

double OceanLoad3D_crust1::getOceanDepth(double theta, double phi) const {
    // convert theta to co-latitude  
    if (mGeographic) {
        theta = pi / 2. - XMath::theta2Lat(theta, 0.) * degree;
    }
    
    // interpolation on sphere
    double depth = 0.;
    std::vector<size_t> ilat, ilon;
    std::vector<double> wlat, wlon;
    Volumetric3D_crust1::interpThetaPhi(theta, phi, mNPointInterp, ilat, ilon, wlat, wlon);
    for (size_t i = 0; i < mNPointInterp; i++) {
        for (size_t j = 0; j < mNPointInterp; j++) {
            double weight = wlat[i] * wlon[j];
            depth += weight * mDepth(ilat[i], ilon[j]);
        }
    }
    if (mBenchmarkSPECFEM) {
        double elevation = depth;
        double MINIMUM_THICKNESS_3D_OCEANS = 50.;
        if (elevation >= - MINIMUM_THICKNESS_3D_OCEANS) {
            depth = 0.;
        } else {
            depth = std::abs(elevation);
        }
    }
    return depth;
}

std::string OceanLoad3D_crust1::verbose() const {
    std::stringstream ss;
    ss << "\n======================= 3D OceanLoad =======================" << std::endl;
    ss << "  Model Name            =   Crust 1.0" << std::endl;
    ss << "  Smoothing Order       =   " << mGaussianOrder << std::endl;
    ss << "  Smoothing Intensity   =   " << mGaussianDev << std::endl;
    ss << "  Num. Interp. Points   =   " << mNPointInterp << std::endl;
    ss << "  Use Geographic        =   " << (mGeographic ? "YES" : "NO") << std::endl;
    ss << "  Add Ice as Water      =   " << (mIncludeIceAsWater ? "YES" : "NO") << std::endl;
    if (mBenchmarkSPECFEM) {
        ss << "  Using SPECFEM method, computing ocean depth by elevation." << std::endl;
    }
    ss << "======================= 3D OceanLoad =======================\n" << std::endl;
    return ss.str();
}


