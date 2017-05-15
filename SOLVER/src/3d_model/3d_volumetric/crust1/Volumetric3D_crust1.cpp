// Volumetric3D_crust1.h
// created by Kuangdai on 16-May-2016 
// crustal model CRUST 1.0 
// http://igppweb.ucsd.edu/~gabi/crust1.html

#include "Volumetric3D_crust1.h"
#include <sstream>
#include <fstream>
#include "XMPI.h"
#include "XMath.h"
#include "SpectralConstants.h"
#include <algorithm>
#include "ExodusModel.h"

const size_t Volumetric3D_crust1::sNLayer = 9;
const size_t Volumetric3D_crust1::sNLat = 180;
const size_t Volumetric3D_crust1::sNLon = 360;

void Volumetric3D_crust1::initialize() {
    // read raw data
    size_t nrow = sNLat * sNLon;
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
        for (size_t i = 0; i < nrow; i++) {
            for (size_t j = 0; j < sNLayer; j++) {
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
    for (size_t col = 0; col < sNLayer; col++) {
        mRl.block(0, col, sNLon, 1).fill(bnd.block(0, col, sNLon, 1).sum() / sNLon);
        mVp.block(0, col, sNLon, 1).fill(vp_.block(0, col, sNLon, 1).sum() / sNLon);
        mVs.block(0, col, sNLon, 1).fill(vs_.block(0, col, sNLon, 1).sum() / sNLon);
        mRh.block(0, col, sNLon, 1).fill(rho.block(0, col, sNLon, 1).sum() / sNLon);
        mRl.block(nrow, col, sNLon, 1).fill(bnd.block(nrow - sNLon, col, sNLon, 1).sum() / sNLon);
        mVp.block(nrow, col, sNLon, 1).fill(vp_.block(nrow - sNLon, col, sNLon, 1).sum() / sNLon);
        mVs.block(nrow, col, sNLon, 1).fill(vs_.block(nrow - sNLon, col, sNLon, 1).sum() / sNLon);
        mRh.block(nrow, col, sNLon, 1).fill(rho.block(nrow - sNLon, col, sNLon, 1).sum() / sNLon);
    }
    for (size_t i = 1; i < sNLat; i++) {
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
    
    // layers
    if (mIncludeIce) mIncludeSediment = true;
    size_t colSurf = columnSurf();
    size_t colMoho = 8;

    // linear mapping to sphere
    const RDColX &rmoho = RDColX::Constant(nrow + sNLon, mRMoho);
    const RDColX &rdiff = (RDColX::Constant(nrow + sNLon, mRSurf - mRMoho).array() 
        / (mRl.col(colSurf) - mRl.col(colMoho)).array()).matrix();
    RDMatXX copyRl = mRl;    
    for (size_t i = 0; i < sNLayer; i++) {
        mRl.col(i).array() = rdiff.array() * (copyRl.col(i) - copyRl.col(colMoho)).array() + rmoho.array();
    }
    
    // above: physical layers
    ///////////////////////////////////////////////////////
    // blow: elemental layers
    
    // form gll boundaries
    size_t nEleCrust = mElementBoundaries.size() - 1;
    const RDColP &eta = SpectralConstants::getP_GLL();
    size_t numGll = 1 + nPol * nEleCrust;
    mRlGLL = RDColX(numGll);
    mRlGLL(0) = mRSurf; // surface
    size_t index = 1;
    for (size_t iele = 0; iele < nEleCrust; iele++) {
        double eBoundTop = mElementBoundaries[iele];
        double eBoundBot = mElementBoundaries[iele + 1];
        double eHeight = eBoundTop - eBoundBot;
        for (size_t jpol = 1; jpol <= nPol; jpol++) {
            double z = eBoundTop - (eta(jpol) - eta(0)) / (eta(nPol) - eta(0)) * eHeight;
            mRlGLL(index++) = z;
        }
    }
    
    // zero properties
    mVpGLL = mVsGLL = mRhGLL = RDMatXX::Zero(nrow + sNLon, numGll);
    
    // form values at GLL boundaries
    for (size_t igll = 0; igll < numGll; igll++) {
        size_t itop = igll - 1 ? igll >= 1 : 0;
        size_t imid = igll;
        size_t ibot = igll + 1;
        // if (itop < 0) itop = 0;
        if (ibot > numGll - 1) ibot = numGll - 1;
        double gll_top = mRlGLL(itop);
        double gll_mid = mRlGLL(imid);
        double gll_bot = mRlGLL(ibot);
        double gll_mid_value = 2. / (gll_top - gll_bot);
        for (size_t row = 0; row < mRl.rows(); row++) {
            // integrate
            for (size_t ilayer = colSurf; ilayer < colMoho; ilayer++) {
                double phy_top = mRl(row, ilayer);
                double phy_bot = mRl(row, ilayer + 1);
                // top to mid
                double top = std::min(phy_top, gll_top);
                double bot = std::max(phy_bot, gll_mid);
                if (top > bot) {
                    double gllt = gll_mid_value / (gll_top - gll_mid) * (gll_top - top);
                    double gllb = gll_mid_value / (gll_top - gll_mid) * (gll_top - bot);
                    double area = .5 * (gllt + gllb) * (top - bot);
                    mVpGLL(row, igll) += mVp(row, ilayer) * area;
                    mVsGLL(row, igll) += mVs(row, ilayer) * area;
                    mRhGLL(row, igll) += mRh(row, ilayer) * area; 
                }
                // mid to bot
                top = std::min(phy_top, gll_mid);
                bot = std::max(phy_bot, gll_bot);
                if (top > bot) {
                    double gllt = gll_mid_value / (gll_mid - gll_bot) * (top - gll_bot);
                    double gllb = gll_mid_value / (gll_mid - gll_bot) * (bot - gll_bot);
                    double area = .5 * (gllt + gllb) * (top - bot);
                    mVpGLL(row, igll) += mVp(row, ilayer) * area;
                    mVsGLL(row, igll) += mVs(row, ilayer) * area;
                    mRhGLL(row, igll) += mRh(row, ilayer) * area; 
                } 
            }
        }
    }
    
    // double xmax = -1e100;
    // size_t rmax = -1;
    // double xmin = 1e100;
    // size_t rmin = -1;
    // for (size_t row = 0; row < mRl.rows(); row++) {
    //     double xdiff = mRhGLL(row, 0);
    //     if (xdiff > xmax) {
    //         xmax = xdiff;
    //         rmax = row;
    //     }
    //     if (xdiff < xmin) {
    //         xmin = xdiff;
    //         rmin = row;
    //     }
    // }
    // std::cout << rmin << std::endl;
    // std::cout << rmax << std::endl;
    // 
    // exit(0);
    // std::fstream fs;
    // for (size_t row = 0; row < mRl.rows(); row++) {
    //     std::stringstream ss;
    //     ss << "x/vp/" << row << ".txt";
    //     fs.open(ss.str(), std::fstream::out);
    //     fs << mRl.row(row).block(0, colSurf, 1, colMoho - colSurf + 1) << std::endl;
    //     fs << mVp.row(row).block(0, colSurf, 1, colMoho - colSurf + 1) << std::endl;
    //     fs << mRlGLL.transpose() << std::endl;
    //     fs << mVpGLL.row(row) << std::endl << std::endl;
    //     fs.close();
    // }
    // for (size_t row = 0; row < mRl.rows(); row++) {
    //     std::stringstream ss;
    //     ss << "x/vs/" << row << ".txt";
    //     fs.open(ss.str(), std::fstream::out);
    //     fs << mRl.row(row).block(0, colSurf, 1, colMoho - colSurf + 1) << std::endl;
    //     fs << mVs.row(row).block(0, colSurf, 1, colMoho - colSurf + 1) << std::endl;
    //     fs << mRlGLL.transpose() << std::endl;
    //     fs << mVsGLL.row(row) << std::endl << std::endl;
    //     fs.close();
    // }
    // for (size_t row = 0; row < mRl.rows(); row++) {
    //     std::stringstream ss;
    //     ss << "x/rho/" << row << ".txt";
    //     fs.open(ss.str(), std::fstream::out);
    //     fs << mRl.row(row).block(0, colSurf, 1, colMoho - colSurf + 1) << std::endl;
    //     fs << mRh.row(row).block(0, colSurf, 1, colMoho - colSurf + 1) << std::endl;
    //     fs << mRlGLL.transpose() << std::endl;
    //     fs << mRhGLL.row(row) << std::endl << std::endl;
    //     fs.close();
    // }
}

bool Volumetric3D_crust1::get3dProperties(double r, double theta, double phi, double rElemCenter,
    double &vpv, double &vph, double &vsv, double &vsh, double &rho) const {
    // not in crust 
    if (rElemCenter > mRSurf || rElemCenter < mRMoho) {
        vpv = vph = vsv = vsh = rho = 0.;
        return false;
    }    
    
    // horizontal interpolation
    std::vector<size_t> ilat, ilon;
    std::vector<double> wlat, wlon;
    if (mGeographic) theta = pi / 2. - XMath::theta2Lat(theta, mRSurf - r) * degree;
    interpThetaPhi(theta, phi, mNPointInterp, ilat, ilon, wlat, wlon);
    
    // element inside but point slightly outside
    if (r >= mRSurf * 0.999999) r = mRSurf * 0.999999;
    if (r <= mRMoho * 1.000001) r = mRMoho * 1.000001;
    
    vpv = 0.;
    vsv = 0.;
    rho = 0.;
    for (size_t i = 0; i < mNPointInterp; i++) {
        for (size_t j = 0; j < mNPointInterp; j++) {
            double weight = wlat[i] * wlon[j];
            size_t rowdata = ilat[i] * sNLon + ilon[j];
            bool found = false;
            for (size_t iLayer = 1; iLayer < mRlGLL.rows(); iLayer++) {
                if (r > mRlGLL(iLayer)) {
                    double vp_top = mVpGLL(rowdata, iLayer - 1);
                    double vs_top = mVsGLL(rowdata, iLayer - 1);
                    double rh_top = mRhGLL(rowdata, iLayer - 1);
                    double vp_bot = mVpGLL(rowdata, iLayer);
                    double vs_bot = mVsGLL(rowdata, iLayer);
                    double rh_bot = mRhGLL(rowdata, iLayer);
                    double top = mRlGLL(iLayer - 1);
                    double bot = mRlGLL(iLayer);
                    double vp = (vp_top - vp_bot) / (top - bot) * (r - bot) + vp_bot;
                    double vs = (vs_top - vs_bot) / (top - bot) * (r - bot) + vs_bot;
                    double rh = (rh_top - rh_bot) / (top - bot) * (r - bot) + rh_bot;
                    vpv += vp * weight;
                    vsv += vs * weight;
                    rho += rh * weight;
                    found = true;
                    break;
                }
            } 
            if (!found) {
                throw std::runtime_error("Volumetric3D_crust1::get3dProperties || Please report this bug.");
            }
        }
    }
    
    // when sediment is considered, Vs at some locations 
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
    ss << "  Num. Element Layers   =   " << mElementBoundaries.size() - 1 << std::endl;
    ss << "  Num. Interp. Points   =   " << mNPointInterp << std::endl;
    ss << "  Use Geographic        =   " << (mGeographic ? "YES" : "NO") << std::endl;
    ss << "======================= 3D Volumetric ======================\n" << std::endl;
    return ss.str();
}

void Volumetric3D_crust1::interpThetaPhi(double theta, double phi, size_t np, 
    std::vector<size_t> &ilat, std::vector<size_t> &ilon, 
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
    size_t ilon_all[sNLon * 2];
    for (size_t i = 0; i <= sNLat; i++) lat_all[i] = i * 1.0; 
    for (size_t i = 0; i < sNLon * 2; i++) {
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
    for (size_t ip = 0; ip < np; ip++) {
        double distmin = 1e100;
        size_t closest_ilat = -1;
        for (size_t i = 0; i < sNLat + 1; i++) {
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
    for (size_t ip = 0; ip < np; ip++) {
        double distmin = 1e100;
        size_t closest_ilon = -1;
        for (size_t i = 0; i < sNLon * 2; i++) {
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

void Volumetric3D_crust1::setupExodusModel(const ExodusModel *exModel) {
    // find all z-coordinates on the axis
    for (size_t i = 0; i < exModel->getNumNodes(); i++) {
        double s = exModel->getNodalS(i);
        double z = exModel->getNodalZ(i);
        if (s <= tinyDouble && z > 0) mElementBoundaries.push_back(z);
    }
    std::sort(mElementBoundaries.begin(), mElementBoundaries.end(), std::greater<double>());
}

void Volumetric3D_crust1::initialize(const std::vector<std::string> &params) {
    try {
        size_t ipar = 0;
        const std::string source = "Volumetric3D_crust1::initialize";
        XMath::castValue(mIncludeSediment, params.at(ipar++), source);
        XMath::castValue(mNPointInterp, params.at(ipar++), source);
        XMath::castValue(mGeographic, params.at(ipar++), source);
        XMath::castValue(mRMoho, params.at(ipar++), source); mRMoho *= 1e3;
        XMath::castValue(mIncludeIce, params.at(ipar++), source);
    } catch (std::out_of_range) {
        // nothing
    }
    
    // delete element boundaries out of crust
    for (size_t i = 0; i < mElementBoundaries.size(); i++) {
        if (mElementBoundaries[i] < mRMoho - 1.) {
            mElementBoundaries.erase(mElementBoundaries.begin() + i, mElementBoundaries.end());
            break;
        }
    }
    
    // check boundaries
    if (mElementBoundaries.size() < 2) 
        throw std::runtime_error("Volumetric3D_crust1::setupExodusModel || No element layer found in crust.");
    if (std::abs(mRSurf - mElementBoundaries[0]) > 1.)
        throw std::runtime_error("Volumetric3D_crust1::setupExodusModel || Conflict in surface radius.");    
    if (std::abs(mRMoho - mElementBoundaries[mElementBoundaries.size() - 1]) > 1.)
        throw std::runtime_error("Volumetric3D_crust1::setupExodusModel || Conflict in Moho radius.");
        
    // initialize
    initialize();
}

