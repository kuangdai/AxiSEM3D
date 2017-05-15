// Geometric3D_crust1.h
// created by Kuangdai on 17-Jun-2016 
// crustal model CRUST 1.0 
// http://igppweb.ucsd.edu/~gabi/crust1.html

#include "Geometric3D_crust1.h"
#include "Volumetric3D_crust1.h"
#include <sstream>
#include <fstream>
#include "XMPI.h"
#include "XMath.h"

const size_t Geometric3D_crust1::sNLayer = 9;
const size_t Geometric3D_crust1::sNLat = 180;
const size_t Geometric3D_crust1::sNLon = 360;

void Geometric3D_crust1::initialize() {
    // read raw data
    size_t nrow = sNLat * sNLon;
    RDMatXX elevation = RDMatXX::Zero(nrow, sNLayer);
    if (XMPI::root()) {
        std::string fname = projectDirectory + "/src/3d_model/3d_volumetric/crust1/data/crust1.bnds";
        std::fstream fs(fname, std::fstream::in);
        if (!fs) {
            throw std::runtime_error("Geometric3D_crust1::initialize || "
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
    
    // surface and moho undulations
    // NOTE: ellipticity should not be considered here because all geometric 
    //       models are defined independently w.r.t. reference sphere 
    size_t colSurf = 5; // no ice, no sediment
    if (mIncludeIce) {
        colSurf = 1; // ice
        mIncludeSediment = true;
    } else if (mIncludeSediment) {
        colSurf = 2; // sediment
    }
    size_t colMoho = 8;
    RDColX deltaRSurfVec = elevation.col(colSurf) * 1e3;
    RDColX deltaRMohoVec = elevation.col(colMoho) * 1e3;
    deltaRMohoVec += RDColX::Constant(nrow, mRSurf - mRMoho);
    // cast to matrix
    RDMatXX deltaRSurf(sNLat, sNLon);
    RDMatXX deltaRMoho(sNLat, sNLon);
    for (size_t i = 0; i < sNLat; i++) {
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
    // size_t npolar = mGaussianOrder + 1;
    // size_t opolar = sNLon;
    // for (size_t i = 0; i < npolar; i++) {
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
    for (size_t i = 1; i < sNLat; i++) {
        mDeltaRSurf.row(i) = (deltaRSurf.row(i - 1) + deltaRSurf.row(i)) * .5;
        mDeltaRMoho.row(i) = (deltaRMoho.row(i - 1) + deltaRMoho.row(i)) * .5; 
    }
    
    // apply factor
    mDeltaRSurf *= mSurfFactor;
    mDeltaRMoho *= mMohoFactor;
    
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
        size_t ipar = 0;
        const std::string source = "Geometric3D_crust1::initialize";
        XMath::castValue(mIncludeSediment, params.at(ipar++), source);
        XMath::castValue(mSurfFactor, params.at(ipar++), source);
        XMath::castValue(mMohoFactor, params.at(ipar++), source);
        XMath::castValue(mGaussianOrder, params.at(ipar++), source);
        XMath::castValue(mGaussianDev, params.at(ipar++), source);
        XMath::castValue(mNPointInterp, params.at(ipar++), source);
        XMath::castValue(mGeographic, params.at(ipar++), source);
        XMath::castValue(mRBase, params.at(ipar++), source); mRBase *= 1e3;
        XMath::castValue(mRMoho, params.at(ipar++), source); mRMoho *= 1e3;
        XMath::castValue(mIncludeIce, params.at(ipar++), source);
    } catch (std::out_of_range) {
        // nothing
    }
    initialize();
}

double Geometric3D_crust1::getDeltaR(double r, double theta, double phi, double rElemCenter) const {
    if (rElemCenter > mRSurf || rElemCenter < mRBase) { 
        return 0.;
    }    
    
    // convert theta to co-latitude 
    if (mGeographic) {
        theta = pi / 2. - XMath::theta2Lat(theta, mRSurf - r) * degree;
    }
    
    // interpolation on sphere
    double drSurf = 0.;
    double drMoho = 0.;
    std::vector<size_t> ilat, ilon;
    std::vector<double> wlat, wlon;
    Volumetric3D_crust1::interpThetaPhi(theta, phi, mNPointInterp, ilat, ilon, wlat, wlon);
    for (size_t i = 0; i < mNPointInterp; i++) {
        for (size_t j = 0; j < mNPointInterp; j++) {
            double weight = wlat[i] * wlon[j];
            drSurf += weight * mDeltaRSurf(ilat[i], ilon[j]);
            drMoho += weight * mDeltaRMoho(ilat[i], ilon[j]);
        }
    }

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
    ss << "  Num. Interp. Points   =   " << mNPointInterp << std::endl;
    ss << "  Use Geographic        =   " << (mGeographic ? "YES" : "NO") << std::endl;
    ss << "======================= 3D Geometric =======================\n" << std::endl;
    return ss.str();
}

