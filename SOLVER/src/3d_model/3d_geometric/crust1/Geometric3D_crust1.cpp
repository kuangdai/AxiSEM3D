// Geometric3D_crust1.h
// created by Kuangdai on 17-Jun-2016 
// crustal model CRUST 1.0 
// http://igppweb.ucsd.edu/~gabi/crust1.html

#include "Geometric3D_crust1.h"
#include <sstream>
#include <fstream>
#include "XMPI.h"
#include "XMath.h"
#include "Geodesy.h"
#include "Parameters.h"

const int Geometric3D_crust1::sNLayer = 9;
const int Geometric3D_crust1::sNLat = 180;
const int Geometric3D_crust1::sNLon = 360;

void Geometric3D_crust1::initialize() {
    // read raw data
    int nrow = sNLat * sNLon;
    RDMatXX elevation = RDMatXX::Zero(nrow, sNLayer);
    if (XMPI::root()) {
        std::string fname = projectDirectory + "/src/3d_model/3d_volumetric/crust1/data/crust1.bnds";
        std::fstream fs(fname, std::fstream::in);
        if (!fs) {
            throw std::runtime_error("Geometric3D_crust1::initialize || "
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
    
    // surface and moho undulations
    // NOTE: ellipticity should not be considered here because all geometric 
    //       models are defined independently w.r.t. reference sphere 
    int colSurf = 5; // no ice, no sediment
    if (mIncludeIce) {
        colSurf = 1; // ice
        mIncludeSediment = true;
    } else if (mIncludeSediment) {
        colSurf = 2; // sediment
    }
    int colMoho = 8;
    RDColX deltaRSurfVec = elevation.col(colSurf) * 1e3;
    RDColX deltaRMohoVec = elevation.col(colMoho) * 1e3;
    deltaRMohoVec += RDColX::Constant(nrow, mRSurf - mRMoho);
    // cast to matrix
    RDMatXX deltaRSurf(sNLat, sNLon);
    RDMatXX deltaRMoho(sNLat, sNLon);
    for (int i = 0; i < sNLat; i++) {
        deltaRSurf.row(i) = deltaRSurfVec.block(i * sNLon, 0, sNLon, 1).transpose();
        deltaRMoho.row(i) = deltaRMohoVec.block(i * sNLon, 0, sNLon, 1).transpose();
    }

    //////////// plot raw data ////////////  
    // std::fstream fs;
    // fs.open("/Users/kuangdai/Desktop/crust1/raw.txt", std::fstream::out);
    // fs << deltaRMoho << std::endl;
    // fs.close();
    //////////// plot raw data //////////// 
    
    // Gaussian smoothing
    IColX orderRow = IColX::Constant(sNLat, mGaussianOrder);
    IColX orderCol = IColX::Constant(sNLon, mGaussianOrder);
    RDColX devRow = RDColX::Constant(sNLat, mGaussianDev);
    RDColX devCol = RDColX::Constant(sNLon, mGaussianDev);
    // smooth poles more
    // int npolar = mGaussianOrder + 1;
    // int opolar = sNLon;
    // for (int i = 0; i < npolar; i++) {
    //     double order = (double)(opolar - mGaussianOrder) / npolar * (npolar - i) + mGaussianOrder;
    //     orderRow(i) = orderRow(sNLat - 1 - i) = round(order);
    // }
    XMath::gaussianSmoothing(deltaRSurf, orderRow, devRow, true, orderCol, devCol, false);
    XMath::gaussianSmoothing(deltaRMoho, orderRow, devRow, true, orderCol, devCol, false);
    
    // cast to integer theta with unique polar values
    mDeltaRSurf = RDMatXX::Zero(sNLat + 1, sNLon);
    mDeltaRMoho = RDMatXX::Zero(sNLat + 1, sNLon);
    // fill north and south pole
    mDeltaRSurf.row(0).fill(deltaRSurf.row(0).sum() / sNLon);
    mDeltaRMoho.row(0).fill(deltaRMoho.row(0).sum() / sNLon);
    mDeltaRSurf.row(sNLat).fill(deltaRSurf.row(sNLat - 1).sum() / sNLon);
    mDeltaRMoho.row(sNLat).fill(deltaRMoho.row(sNLat - 1).sum() / sNLon);
    // interp at integer theta
    for (int i = 1; i < sNLat; i++) {
        mDeltaRSurf.row(i) = (deltaRSurf.row(i - 1) + deltaRSurf.row(i)) * .5;
        mDeltaRMoho.row(i) = (deltaRMoho.row(i - 1) + deltaRMoho.row(i)) * .5; 
    }
    // reverse south to north
    mDeltaRSurf = mDeltaRSurf.colwise().reverse().eval();
    mDeltaRMoho = mDeltaRMoho.colwise().reverse().eval();
    
    // apply factor
    mDeltaRSurf *= mSurfFactor;
    mDeltaRMoho *= mMohoFactor;
    
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
    // fs.open("/Users/kuangdai/Desktop/crust1/smt.txt", std::fstream::out);
    // fs << deltaRMoho << std::endl;
    // fs.close();
    //////////// plot raw data //////////// 
    
    //////////// plot computed data ////////////  
    // std::fstream fsdr;
    // fsdr.open("/Users/kuangdai/Desktop/crust1/drmoho1.txt", std::fstream::out);
    // double r = mRMoho; // double r = mRSurf;
    // int intGrid = 4;
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
    //////////// plot computed data ////////////  
    
    //////////// generate SPECFEM topo_bathy ////////////  
    // std::cout << verbose() << std::endl;
    // mGeographic = false;
    // std::fstream fsdr;
    // fsdr.open("/Users/kuangdai/Desktop/crust1/specfem_topo_bathy.txt", std::fstream::out);
    // double r = mRSurf; 
    // int intGrid = 15;
    // for (int i = 0; i <= 180 * intGrid; i++) {
    //     double theta = i * degree / intGrid;
    //     for (int j = 0; j < 360 * intGrid; j++) {
    //         double phi = j * degree / intGrid;
    //         fsdr << round(getDeltaR(r, theta, phi, r)) << std::endl;
    //     }
    // }   
    // fsdr.close();
    // exit(0);
    //////////// generate SPECFEM topo_bathy ////////////  
}

void Geometric3D_crust1::initialize(const std::vector<std::string> &params) {
    try {
        int ipar = 0;
        const std::string source = "Geometric3D_crust1::initialize";
        Parameters::castValue(mIncludeSediment, params.at(ipar++), source);
        Parameters::castValue(mSurfFactor, params.at(ipar++), source);
        Parameters::castValue(mMohoFactor, params.at(ipar++), source);
        Parameters::castValue(mGaussianOrder, params.at(ipar++), source);
        Parameters::castValue(mGaussianDev, params.at(ipar++), source);
        Parameters::castValue(mGeographic, params.at(ipar++), source);
        Parameters::castValue(mRBase, params.at(ipar++), source); mRBase *= 1e3;
        Parameters::castValue(mRMoho, params.at(ipar++), source); mRMoho *= 1e3;
        Parameters::castValue(mIncludeIce, params.at(ipar++), source);
    } catch (std::out_of_range) {
        // nothing
    }
    mRSurf = Geodesy::getROuter();
    initialize();
}

double Geometric3D_crust1::getDeltaR(double r, double theta, double phi, double rElemCenter) const {
    if (rElemCenter > mRSurf || rElemCenter < mRBase) { 
        return 0.;
    }    
    
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
    
    double drSurf = 0.;
    drSurf += mDeltaRSurf(llat0, llon0) * wlat0 * wlon0;
    drSurf += mDeltaRSurf(llat1, llon0) * wlat1 * wlon0;
    drSurf += mDeltaRSurf(llat0, llon1) * wlat0 * wlon1;
    drSurf += mDeltaRSurf(llat1, llon1) * wlat1 * wlon1;
    
    double drMoho = 0.;
    drMoho += mDeltaRMoho(llat0, llon0) * wlat0 * wlon0;
    drMoho += mDeltaRMoho(llat1, llon0) * wlat1 * wlon0;
    drMoho += mDeltaRMoho(llat0, llon1) * wlat0 * wlon1;
    drMoho += mDeltaRMoho(llat1, llon1) * wlat1 * wlon1;

    // interpolation along radius    
    if (rElemCenter < mRMoho) {
        return drMoho / (mRMoho - mRBase) * (r - mRBase);
    } else {
        return (drSurf - drMoho) / (mRSurf - mRMoho) * (r - mRMoho) + drMoho;
    } 
}

std::string Geometric3D_crust1::verbose() const {
    std::stringstream ss;
    ss << "\n======================= 3D Geometric =======================" << std::endl;
    ss << "  Model Name            =   Crust 1.0" << std::endl;
    ss << "  Scope                 =   crust" << std::endl;
    ss << "  Radii (km)            =   [" << mRBase / 1e3 << ", " << mRSurf / 1e3 << "]" << std::endl;
    ss << "  RMoho (km)            =   " << mRMoho / 1e3 << std::endl;
    ss << "  Include Ice           =   " << (mIncludeIce ? "YES" : "NO") << std::endl;
    ss << "  Include Sediment      =   " << (mIncludeSediment ? "YES" : "NO") << std::endl;
    ss << "  Surface Factor        =   " << mSurfFactor << std::endl;
    ss << "  Moho Factor           =   " << mMohoFactor << std::endl;
    ss << "  Smoothing Order       =   " << mGaussianOrder << std::endl;
    ss << "  Smoothing Intensity   =   " << mGaussianDev << std::endl;
    ss << "  Use Geographic        =   " << (mGeographic ? "YES" : "NO") << std::endl;
    ss << "======================= 3D Geometric =======================\n" << std::endl;
    return ss.str();
}

