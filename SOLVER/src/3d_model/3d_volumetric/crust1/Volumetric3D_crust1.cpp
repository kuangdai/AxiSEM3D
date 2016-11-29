// Volumetric3D_crust1.h
// created by Kuangdai on 16-May-2016 
// crustal model CRUST 1.0 
// http://igppweb.ucsd.edu/~gabi/crust1.html

#include "Volumetric3D_crust1.h"
#include <sstream>
#include <fstream>
#include "XMPI.h"
#include "XMath.h"

const int Volumetric3D_crust1::sNLayer = 9;
const int Volumetric3D_crust1::sNLat = 180;
const int Volumetric3D_crust1::sNLon = 360;

void Volumetric3D_crust1::initialize() {
    // read raw data
    int nrow = sNLat * sNLon;
    RDMatXX bnd, vp_, vs_, rho;
    bnd = vp_ = vs_ = rho = RDMatXX::Zero(nrow, sNLayer);
    if (XMPI::root()) {
        std::string path = projectDirectory + "/src/3d_model/3d_volumetric/crust1/data";
        std::fstream fsbnd, fsvp_, fsvs_, fsrho;
        fsbnd.open(path + "/crust1.bnds", std::fstream::in);
        fsvp_.open(path + "/crust1.vp", std::fstream::in);
        fsvs_.open(path + "/crust1.vs", std::fstream::in);
        fsrho.open(path + "/crust1.rho", std::fstream::in);
        if (!fsbnd || !fsvp_ || !fsvs_ || !fsrho) 
            throw std::runtime_error("Volumetric3D_crust1::initialize || "
                "Error opening crust1.0 data files at directory: ||" + path);
        for (int i = 0; i < nrow; i++) {
            for (int j = 0; j < sNLayer; j++) {
                fsbnd >> bnd(i, j);
                fsvp_ >> vp_(i, j);
                fsvs_ >> vs_(i, j);
                fsrho >> rho(i, j);
            }
        }
        fsbnd.close();
        fsvp_.close();
        fsvs_.close();
        fsrho.close();    
    }
    // broadcast
    XMPI::bcastEigen(bnd);
    XMPI::bcastEigen(vp_);
    XMPI::bcastEigen(vs_);
    XMPI::bcastEigen(rho);
    
    // cast to integer theta
    mRl = mVp = mVs = mRh = RDMatXX::Zero(nrow + sNLon, sNLayer);
    for (int col = 0; col < sNLayer; col++) {
        mRl.block(0, col, sNLon, 1).fill(bnd.block(0, col, sNLon, 1).sum() / sNLon);
        mVp.block(0, col, sNLon, 1).fill(vp_.block(0, col, sNLon, 1).sum() / sNLon);
        mVs.block(0, col, sNLon, 1).fill(vs_.block(0, col, sNLon, 1).sum() / sNLon);
        mRh.block(0, col, sNLon, 1).fill(rho.block(0, col, sNLon, 1).sum() / sNLon);
        mRl.block(nrow, col, sNLon, 1).fill(bnd.block(nrow - sNLon, col, sNLon, 1).sum() / sNLon);
        mVp.block(nrow, col, sNLon, 1).fill(vp_.block(nrow - sNLon, col, sNLon, 1).sum() / sNLon);
        mVs.block(nrow, col, sNLon, 1).fill(vs_.block(nrow - sNLon, col, sNLon, 1).sum() / sNLon);
        mRh.block(nrow, col, sNLon, 1).fill(rho.block(nrow - sNLon, col, sNLon, 1).sum() / sNLon);
    }
    for (int i = 1; i < sNLat; i++) {
        mRl.block(i * sNLon, 0, sNLon, sNLayer) = (bnd.block(i * sNLon, 0, sNLon, sNLayer) + bnd.block((i - 1) * sNLon, 0, sNLon, sNLayer)) * .5;
        mVp.block(i * sNLon, 0, sNLon, sNLayer) = (vp_.block(i * sNLon, 0, sNLon, sNLayer) + vp_.block((i - 1) * sNLon, 0, sNLon, sNLayer)) * .5;
        mVs.block(i * sNLon, 0, sNLon, sNLayer) = (vs_.block(i * sNLon, 0, sNLon, sNLayer) + vs_.block((i - 1) * sNLon, 0, sNLon, sNLayer)) * .5;
        mRh.block(i * sNLon, 0, sNLon, sNLayer) = (rho.block(i * sNLon, 0, sNLon, sNLayer) + rho.block((i - 1) * sNLon, 0, sNLon, sNLayer)) * .5;
    } 
    
    // convert to SI
    mVp *= 1e3;
    mVs *= 1e3;
    mRh *= 1e3;
    mRl *= 1e3;
    
    // check sediment thickness to avoid too thin sediment layers
    for (int row = 0; row < mRl.rows(); row++) {
        // 3 sediments + ice, check from bottom to up
        for (int ised = 4; ised >= 1; ised--) {
            double sed = mRl(row, ised) - mRl(row, ised + 1);
            if (sed > 1. && sed < mMinimumSedimentThickness) {
                mVp(row, ised) = mVp(row, ised + 1);
                mVs(row, ised) = mVs(row, ised + 1);
                mRh(row, ised) = mRh(row, ised + 1);
            }
        }
    }

    // linear mapping to sphere
    int colSurf = 5; // no ice, no sediment
    if (mIncludeIce) {
        colSurf = 1; // ice
        mIncludeSediment = true;
    } else if (mIncludeSediment) {
        colSurf = 2; // sediment
    }
    int colMoho = 8;
    const RDColX &rmoho = RDColX::Constant(nrow + sNLon, mRMoho);
    const RDColX &rdiff = (RDColX::Constant(nrow + sNLon, mRSurf - mRMoho).array() 
        / (mRl.col(colSurf) - mRl.col(colMoho)).array()).matrix();
    RDMatXX copyRl = mRl;    
    for (int i = 0; i < sNLayer; i++) {
        mRl.col(i).array() = rdiff.array() * (copyRl.col(i) - copyRl.col(colMoho)).array() + rmoho.array();
    }
}

void Volumetric3D_crust1::initialize(const std::vector<double> &params) {
    try {
        int ipar = 0;
        mIncludeSediment = (params.at(ipar++) > tinyDouble);
        mMinimumSedimentThickness = params.at(ipar++) * 1e3;
        mNPointInterp = round(params.at(ipar++));
        mGeographic = (params.at(ipar++) > tinyDouble);
        mRMoho = params.at(ipar++) * 1e3;
        mRSurf = params.at(ipar++) * 1e3;
        mIncludeIce = (params.at(ipar++) > tinyDouble);
    } catch (std::out_of_range) {
        // nothing
    }
    initialize();
}

bool Volumetric3D_crust1::get3dProperties(double r, double theta, double phi, double rElemCenter,
    double &vpv, double &vph, double &vsv, double &vsh, double &rho) const {
    // not in crust 
    if (rElemCenter > mRSurf || rElemCenter < mRMoho) {
        vpv = vph = vsv = vsh = rho = 0.;
        return false;
    }    
    // interpolation
    std::vector<int> ilat, ilon;
    std::vector<double> wlat, wlon;
    if (mGeographic) theta = pi / 2. - XMath::theta2Lat(theta, mRSurf - r) * degree;
    interpThetaPhi(theta, phi, mNPointInterp, ilat, ilon, wlat, wlon);
    
    // element inside but point slightly outside
    if (r >= mRSurf * 0.999999) r = mRSurf * 0.999999;
    if (r <= mRMoho * 1.000001) r = mRMoho * 1.000001;
    
    vpv = 0.;
    vsv = 0.;
    rho = 0.;
    int istart = mIncludeSediment ? 3 : 6; 
    for (int i = 0; i < mNPointInterp; i++) {
        for (int j = 0; j < mNPointInterp; j++) {
            double weight = wlat[i] * wlon[j];
            int rowdata = ilat[i] * sNLon + ilon[j];
            bool found = false;
            for (int iLayer = istart; iLayer <= 8; iLayer++) {
                if (r >= mRl(rowdata, iLayer)) {
                    double vp_read = mVp(rowdata, iLayer - 1);
                    double vs_read = mVs(rowdata, iLayer - 1);
                    double rh_read = mRh(rowdata, iLayer - 1);
                    // make sure we don't read a zero thickness layer
                    if (vp_read < 1. || vs_read < 1. || rh_read < 1.)
                        throw std::runtime_error("Volumetric3D_crust1::get3dProperties || Invalid radius.");
                    vpv += vp_read * weight;
                    vsv += vs_read * weight;
                    rho += rh_read * weight;
                    found = true;
                    break;
                }
            } 
            if (!found) throw std::runtime_error("Volumetric3D_crust1::get3dProperties || Invalid radius.");
        }
    }
    // when sediment is considered, Vs at some spots  
    // could be too small to keep the simulation stable
    vsv = std::max(vsv, 500.);
    vpv = std::max(vpv, sqrt(2) * vsv);
    // isotropic
    vph = vpv;
    vsh = vsv;
    return true;
}

std::string Volumetric3D_crust1::verbose() const {
    std::stringstream ss;
    ss << "\n======================= 3D Volumetric ======================" << std::endl;
    ss << "  Model Name            =   Crust 1.0" << std::endl;
    ss << "  Scope                 =   Crust" << std::endl;
    ss << "  Radii (km)            =   [" << mRMoho / 1e3 << ", " << mRSurf / 1e3 << "]" << std::endl;
    ss << "  Reference Type        =   Absolute" << std::endl;
    ss << "  Max. Fourier Order    =   180" << std::endl;
    ss << "  Anisotropic           =   NO" << std::endl;
    ss << "  Include Ice           =   " << (mIncludeIce ? "YES" : "NO") << std::endl;
    ss << "  Include Sediment      =   " << (mIncludeSediment ? "YES" : "NO") << std::endl;
    ss << "  Min. Sediment / km    =   " << mMinimumSedimentThickness / 1e3 << std::endl;
    ss << "  Num. Interp. Points   =   " << mNPointInterp << std::endl;
    ss << "  Use Geographic        =   " << (mGeographic ? "YES" : "NO") << std::endl;
    ss << "======================= 3D Volumetric ======================\n" << std::endl;
    return ss.str();
}

void Volumetric3D_crust1::interpThetaPhi(double theta, double phi, int np, 
    std::vector<int> &ilat, std::vector<int> &ilon, 
    std::vector<double> &wlat, std::vector<double> &wlon) {
    
    // target colat and lon
    double lat = theta / degree; 
    if (lat < 0.) lat = 0.;
    if (lat > 180.) lat = 180.;
    double lon = phi / degree;
    if (lon < 0.) lon = 0.;
    if (lon > 360.) lon = 360.;
    if (lon > 180.) lon -= 360.;
    
    // populate model sampling points
    // two cycles are needed for longitude
    double lat_all[sNLat + 1], lon_all[sNLon * 2]; 
    int ilon_all[sNLon * 2];
    for (int i = 0; i <= sNLat; i++) lat_all[i] = i * 1.0; 
    for (int i = 0; i < sNLon * 2; i++) {
        lon_all[i] = i * 1.0 - 359.5;
        if (i < sNLon / 2) ilon_all[i] = i + sNLon / 2;
        else if (i < sNLon / 2 * 3) ilon_all[i] = i - sNLon / 2;
        else ilon_all[i] = i - sNLon / 2 * 3; 
    }

    // init output
    ilat.reserve(np);
    ilon.reserve(np);
    wlat.reserve(np);
    wlon.reserve(np);
    
    // find closest lat
    std::vector<double> closest_lats;
    for (int ip = 0; ip < np; ip++) {
        double distmin = 1e100;
        int closest_ilat = -1;
        for (int i = 0; i < sNLat + 1; i++) {
            double dist = std::abs(lat - lat_all[i]);
            if (dist < distmin) {
                distmin = dist;
                closest_ilat = i;
            }
        }
        ilat.push_back(closest_ilat);
        closest_lats.push_back(lat_all[closest_ilat]);
        lat_all[closest_ilat] = 1e100; // set to crazy value
    }
    XMath::interpLagrange(lat, np, closest_lats.data(), wlat.data());
    
    // find closest lon
    std::vector<double> closest_lons;
    for (int ip = 0; ip < np; ip++) {
        double distmin = 1e100;
        int closest_ilon = -1;
        for (int i = 0; i < sNLon * 2; i++) {
            double dist = std::abs(lon - lon_all[i]);
            if (dist < distmin) {
                distmin = dist;
                closest_ilon = i;
            }
        }
        ilon.push_back(ilon_all[closest_ilon]);
        closest_lons.push_back(lon_all[closest_ilon]);
        lon_all[closest_ilon] = 1e100; // set to crazy value
    }
    XMath::interpLagrange(lon, np, closest_lons.data(), wlon.data());
}


