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
        if (!fs) throw std::runtime_error("Geometric3D_crust1::initialize || "
            "Error opening crust1.0 data file: ||" + fname);
        for (int i = 0; i < nrow; i++)
            for (int j = 0; j < sNLayer; j++) fs >> elevation(i, j);    
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
    
    // apply factor
    mDeltaRSurf *= mSurfFactor;
    mDeltaRMoho *= mMohoFactor;
    
    // compute polar values in Cartesian
    // computePolar();
    
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
    // 
    // std::fstream fsdr;
    // fsdr.open("/Users/kuangdai/Desktop/crust1/specfem_topo_bathy.txt", std::fstream::out);
    // double r = mRSurf; 
    // int intGrid = 15;
    // for (int i = 0; i < 180 * intGrid; i++) {
    //     double theta = (i + .5) * degree / intGrid;
    //     for (int j = 0; j < 360 * intGrid; j++) {
    //         double phi = j * degree / intGrid;
    //         fsdr << round(getDeltaR(r, theta, phi, r)) << std::endl;
    //     }
    // }   
    // fsdr.close();
    // exit(0);
    //////////// generate SPECFEM topo_bathy ////////////  
}

void Geometric3D_crust1::initialize(const std::vector<double> &params) {
    try {
        int ipar = 0;
        mIncludeSediment = (params.at(ipar++) > tinyDouble);
        mSurfFactor = params.at(ipar++);
        mMohoFactor = params.at(ipar++);
        mGaussianOrder = round(params.at(ipar++));
        mGaussianDev = params.at(ipar++);
        mNPointInterp = round(params.at(ipar++));
        mGeographic = (params.at(ipar++) > tinyDouble);
        mRBase = params.at(ipar++) * 1e3;
        mRMoho = params.at(ipar++) * 1e3;
        mRSurf = params.at(ipar++) * 1e3;
        mIncludeIce = (params.at(ipar++) > tinyDouble);
    } catch (std::out_of_range) {
        // nothing
    }
    initialize();
}

double Geometric3D_crust1::getDeltaR(double r, double theta, double phi, double rElemCenter) const {
    if (rElemCenter > mRSurf || rElemCenter < mRBase) { 
        return 0.;
    }    
    
    // interpolation on sphere
    double drSurf = 0.;
    double drMoho = 0.;
    std::vector<int> ilat, ilon;
    std::vector<double> wlat, wlon;
    if (mGeographic) theta = pi / 2. - XMath::theta2Lat(theta, mRSurf - r) * degree;
    Volumetric3D_crust1::interpThetaPhi(theta, phi, mNPointInterp, ilat, ilon, wlat, wlon);
    for (int i = 0; i < mNPointInterp; i++) {
        for (int j = 0; j < mNPointInterp; j++) {
            double weight = wlat[i] * wlon[j];
            drSurf += weight * mDeltaRSurf(ilat[i], ilon[j]);
            drMoho += weight * mDeltaRMoho(ilat[i], ilon[j]);
        }
    }

    // interpolation along radius    
    if (rElemCenter <= mRMoho)
        return drMoho / (mRMoho - mRBase) * (r - mRBase);
    else 
        return (drSurf - drMoho) / (mRSurf - mRMoho) * (r - mRMoho) + drMoho;
}

// bool Geometric3D_crust1::getNablaDeltaR(double r, double theta, double phi, double rElemCenter,
//     double &deltaR_r, double &deltaR_theta, double &deltaR_phi) const {
//     double small = 0.01 * degree;
//     if (theta >= small && theta <= pi - small) {
//         return Geometric3D::getNablaDeltaR(r, theta, phi, rElemCenter, 
//             deltaR_r, deltaR_theta, deltaR_phi);
//     }
//     
//     if (rElemCenter < mRBase) { // no need to check r > mRSurf
//         deltaR_r = deltaR_theta = deltaR_phi = 0.;
//         return false;
//     }
//        
//     bool north = theta < small;
//     double drSurf = 0.;
//     double drMoho = 0.;
//     double drdtSurf = 0., drdpSurf = 0.;
//     double drdtMoho = 0., drdpMoho = 0.;
//     if (north) {
//         drSurf = mDeltaRSurf(0, 0);
//         drMoho = mDeltaRMoho(0, 0);
//         drdtSurf =  mDrDxNorthSurf(0) * cos(phi) + mDrDxNorthSurf(1) * sin(phi);
//         drdpSurf = -mDrDxNorthSurf(0) * sin(phi) + mDrDxNorthSurf(1) * cos(phi);
//         drdtMoho =  mDrDxNorthMoho(0) * cos(phi) + mDrDxNorthMoho(1) * sin(phi);
//         drdpMoho = -mDrDxNorthMoho(0) * sin(phi) + mDrDxNorthMoho(1) * cos(phi);
//     } else {
//         drSurf = mDeltaRSurf(sNLat, 0);
//         drMoho = mDeltaRMoho(sNLat, 0);
//         drdtSurf = -mDrDxSouthSurf(0) * cos(phi) - mDrDxSouthSurf(1) * sin(phi);
//         drdpSurf = -mDrDxSouthSurf(0) * sin(phi) + mDrDxSouthSurf(1) * cos(phi);
//         drdtMoho = -mDrDxSouthMoho(0) * cos(phi) - mDrDxSouthMoho(1) * sin(phi);
//         drdpMoho = -mDrDxSouthMoho(0) * sin(phi) + mDrDxSouthMoho(1) * cos(phi);
//     }
//     
//     if (rElemCenter <= mRMoho) {
//         deltaR_r = drMoho / (mRMoho - mRBase);
//         deltaR_theta = drdtMoho / (mRMoho - mRBase) * (r - mRBase);
//         deltaR_phi = drdpMoho / (mRMoho - mRBase) * (r - mRBase);
//     } else {
//         deltaR_r = (drSurf - drMoho) / (mRSurf - mRMoho);
//         deltaR_theta = (drdtSurf - drdtMoho) / (mRSurf - mRMoho) * (r - mRMoho) + drdtMoho;
//         deltaR_phi = (drdpSurf - drdpMoho) / (mRSurf - mRMoho) * (r - mRMoho) + drdpMoho;    
//     }
//     return true;
// }

std::string Geometric3D_crust1::verbose() const {
    std::stringstream ss;
    ss << "\n======================= 3D Geometric ======================" << std::endl;
    ss << "  Model Name            =   Crust 1.0" << std::endl;
    ss << "  Scope                 =   crust" << std::endl;
    ss << "  Radii (km)            =   [" << mRBase / 1e3 << ", " << mRSurf / 1e3 << "]" << std::endl;
    ss << "  RMoho (km)            =   " << mRMoho << std::endl;
    ss << "  Include Ice           =   " << (mIncludeIce ? "YES" : "NO") << std::endl;
    ss << "  Include Sediment      =   " << (mIncludeSediment ? "YES" : "NO") << std::endl;
    ss << "  Surface Factor        =   " << mSurfFactor << std::endl;
    ss << "  Moho Factor           =   " << mMohoFactor << std::endl;
    ss << "  Smoothing Order       =   " << mGaussianOrder << std::endl;
    ss << "  Smoothing Intensity   =   " << mGaussianDev << std::endl;
    ss << "  Num. Interp. Points   =   " << mNPointInterp << std::endl;
    ss << "  Use Geographic        =   " << (mGeographic ? "YES" : "NO") << std::endl;
    ss << "======================= 3D Geometric ======================\n" << std::endl;
    return ss.str();
}

// void Geometric3D_crust1::computePolar() {
//     mDrDxNorthSurf.setZero();
//     mDrDxSouthSurf.setZero();
//     mDrDxNorthMoho.setZero();
//     mDrDxSouthMoho.setZero();
//     double small = 0.01 * degree;
//     double theta_north = small;
//     double theta_south = pi - small;
//     for (int iphi = 0; iphi < 360; iphi++) {
//         double phi = iphi * degree;
//         double dr;
//         double dt_north_surf, dp_north_surf;
//         double dt_north_moho, dp_north_moho;
//         double dt_south_surf, dp_south_surf;
//         double dt_south_moho, dp_south_moho;
//         Geometric3D::getNablaDeltaR(mRSurf, theta_north, phi, mRSurf - 1., dr, dt_north_surf, dp_north_surf);
//         Geometric3D::getNablaDeltaR(mRMoho, theta_north, phi, mRMoho - 1., dr, dt_north_moho, dp_north_moho);
//         Geometric3D::getNablaDeltaR(mRSurf, theta_south, phi, mRSurf - 1., dr, dt_south_surf, dp_south_surf);
//         Geometric3D::getNablaDeltaR(mRMoho, theta_south, phi, mRMoho - 1., dr, dt_south_moho, dp_south_moho);
//         mDrDxNorthSurf(0) +=  dt_north_surf * cos(phi) - dp_north_surf * sin(phi);
//         mDrDxNorthSurf(1) +=  dt_north_surf * sin(phi) + dp_north_surf * cos(phi);
//         mDrDxNorthMoho(0) +=  dt_north_moho * cos(phi) - dp_north_moho * sin(phi);
//         mDrDxNorthMoho(1) +=  dt_north_moho * sin(phi) + dp_north_moho * cos(phi);
//         mDrDxSouthSurf(0) += -dt_south_surf * cos(phi) - dp_south_surf * sin(phi);
//         mDrDxSouthSurf(1) += -dt_south_surf * sin(phi) + dp_south_surf * cos(phi);
//         mDrDxSouthMoho(0) += -dt_south_moho * cos(phi) - dp_south_moho * sin(phi);
//         mDrDxSouthMoho(1) += -dt_south_moho * sin(phi) + dp_south_moho * cos(phi);
//     }
//     mDrDxNorthSurf /= 360.;
//     mDrDxSouthSurf /= 360.;
//     mDrDxNorthMoho /= 360.;
//     mDrDxSouthMoho /= 360.;
// }

