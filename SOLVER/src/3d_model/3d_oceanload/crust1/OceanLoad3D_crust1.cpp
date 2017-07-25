// OceanLoad3D_crust1.h
// created by Kuangdai on 8-Oct-2016 
// crustal model CRUST 1.0 
// http://igppweb.ucsd.edu/~gabi/crust1.html

#include "OceanLoad3D_crust1.h"
#include <sstream>
#include <fstream>
#include "XMPI.h"
#include "XMath.h"
#include "Geodesy.h"
#include "Parameters.h"

const int OceanLoad3D_crust1::sNLayer = 9;
const int OceanLoad3D_crust1::sNLat = 180;
const int OceanLoad3D_crust1::sNLon = 360;

void OceanLoad3D_crust1::initialize() {
    // read raw data
    int nrow = sNLat * sNLon;
    RDMatXX elevation = RDMatXX::Zero(nrow, sNLayer);
    if (XMPI::root()) {
        std::string fname = projectDirectory + "/src/3d_model/3d_volumetric/crust1/data/crust1.bnds";
        std::fstream fs(fname, std::fstream::in);
        if (!fs) {
            throw std::runtime_error("OceanLoad3D_crust1::initialize || "
                "Error opening crust1.0 data file: ||" + fname);
        }
        for (int i = 0; i < nrow; i++) {
            for (int j = 0; j < sNLayer; j++) {
                fs >> elevation(i, j);
            }
        }
        fs.close();
    }
    // broadcast
    XMPI::bcastEigen(elevation);
    
    // water depth
    int colWaterBot = 1;
    if (mIncludeIceAsWater) {
        colWaterBot = 2;
    }
    RDColX depthVec = (elevation.col(0) - elevation.col(colWaterBot)) * 1e3;
    
    // determine ocean depth ONLY by elevation for SPECFEM benchmark
    if (mBenchmarkSPECFEM) {
        depthVec = elevation.col(2) * 1e3;
        mGaussianOrder = 0;
        mIncludeIceAsWater = false;
        mGeographic = false;
    }
    
    RDMatXX depth(sNLat, sNLon);
    for (int i = 0; i < sNLat; i++) {
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
    for (int i = 1; i < sNLat; i++) {
        mDepth.row(i) = (depth.row(i - 1) + depth.row(i)) * .5;
    }
    // reverse south to north
    mDepth = mDepth.colwise().reverse().eval();
    
    // grid lat and lon
    mGridLat = RDColX(sNLat + 1);
    mGridLon = RDColX(sNLon + 1); // one bigger than data
    for (int i = 0; i < sNLat + 1; i++) {
        mGridLat[i] = i * 1. - 90.;
    }
    for (int i = 0; i < sNLon + 1; i++) {
        mGridLon[i] = i * 1. - 179.5;
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
        int ipar = 0;
        const std::string source = "OceanLoad3D_crust1::initialize";
        Parameters::castValue(mIncludeIceAsWater, params.at(ipar++), source);
        Parameters::castValue(mGaussianOrder, params.at(ipar++), source);
        Parameters::castValue(mGaussianDev, params.at(ipar++), source);
        Parameters::castValue(mGeographic, params.at(ipar++), source);
        Parameters::castValue(mBenchmarkSPECFEM, params.at(ipar++), source);
    } catch (std::out_of_range) {
        // nothing
    }
    initialize();
}

double OceanLoad3D_crust1::getOceanDepth(double theta, double phi) const {
    // convert theta to co-latitude  
    if (mGeographic) {
        theta = pi / 2. - Geodesy::theta2Lat_d(theta, 0.) * degree;
    }
    
    // regularise
    double lat = 90. - theta / degree;
    double lon = phi / degree;
    if (lon > 180.) {
        lon -= 360.;
    }
    XMath::checkLimits(lat, -90., 90.);
    XMath::checkLimits(lon, -180., 180.);
    if (lon < -179.5) {
        lon += 360.; 
    }
    
    // interpolation on sphere
    int llat0, llon0, llat1, llon1;
    double wlat0, wlon0, wlat1, wlon1;
    XMath::interpLinear(lat, mGridLat, llat0, wlat0);
    XMath::interpLinear(lon, mGridLon, llon0, wlon0);    
    llat1 = llat0 + 1;
    llon1 = llon0 + 1;
    wlat1 = 1. - wlat0;
    wlon1 = 1. - wlon0;
    if (llon1 == sNLon) {
        llon1 = 0;
    }
    
    double depth = 0.;
    depth += mDepth(llat0, llon0) * wlat0 * wlon0;
    depth += mDepth(llat1, llon0) * wlat1 * wlon0;
    depth += mDepth(llat0, llon1) * wlat0 * wlon1;
    depth += mDepth(llat1, llon1) * wlat1 * wlon1;
    
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
    ss << "  Add Ice as Water      =   " << (mIncludeIceAsWater ? "YES" : "NO") << std::endl;
    ss << "  Smoothing Order       =   " << mGaussianOrder << std::endl;
    ss << "  Smoothing Intensity   =   " << mGaussianDev << std::endl;
    ss << "  Use Geographic        =   " << (mGeographic ? "YES" : "NO") << std::endl;
    if (mBenchmarkSPECFEM) {
        ss << "  Using SPECFEM method, computing ocean depth by elevation." << std::endl;
    }
    ss << "======================= 3D OceanLoad =======================\n" << std::endl;
    return ss.str();
}


