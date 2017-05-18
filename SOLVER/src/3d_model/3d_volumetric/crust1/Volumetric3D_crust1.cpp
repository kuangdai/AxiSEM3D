// Volumetric3D_crust1.h
// created by Kuangdai on 16-May-2016 
// crustal model CRUST 1.0 
// http://igppweb.ucsd.edu/~gabi/crust1.html

#include "Volumetric3D_crust1.h"
#include <sstream>
#include <fstream>
#include "XMPI.h"
#include "XMath.h"
#include "Geodesy.h"
#include "Parameters.h"
#include "SpectralConstants.h"
#include <algorithm>
#include "ExodusModel.h"

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
        if (!fsbnd || !fsvp_ || !fsvs_ || !fsrho) {
            throw std::runtime_error("Volumetric3D_crust1::initialize || "
                "Error opening crust1.0 data files at directory: ||" + path);
        }
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
    
    // layers
    if (mIncludeIce) {
        mIncludeSediment = true;
    }
    int colSurf = columnSurf();
    int colMoho = 8;

    // linear mapping to sphere
    const RDColX &rmoho = RDColX::Constant(nrow + sNLon, mRMoho);
    const RDColX &rdiff = (RDColX::Constant(nrow + sNLon, mRSurf - mRMoho).array() 
        / (mRl.col(colSurf) - mRl.col(colMoho)).array()).matrix();
    RDMatXX copyRl = mRl;    
    for (int i = 0; i < sNLayer; i++) {
        mRl.col(i).array() = rdiff.array() * (copyRl.col(i) - copyRl.col(colMoho)).array() + rmoho.array();
    }
    
    // above: physical layers
    ///////////////////////////////////////////////////////
    // blow: elemental layers
    
    // form gll boundaries
    int nEleCrust = mElementBoundaries.size() - 1;
    const RDColP &eta = SpectralConstants::getP_GLL();
    int numGll = 1 + nPol * nEleCrust;
    mRlGLL = RDColX(numGll);
    mRlGLL(0) = mRSurf; // surface
    int index = 1;
    for (int iele = 0; iele < nEleCrust; iele++) {
        double eBoundTop = mElementBoundaries[iele];
        double eBoundBot = mElementBoundaries[iele + 1];
        double eHeight = eBoundTop - eBoundBot;
        for (int jpol = 1; jpol <= nPol; jpol++) {
            double z = eBoundTop - (eta(jpol) - eta(0)) / (eta(nPol) - eta(0)) * eHeight;
            mRlGLL(index++) = z;
        }
    }
    
    // zero properties
    mVpGLL = mVsGLL = mRhGLL = RDMatXX::Zero(nrow + sNLon, numGll);
    
    // form values at GLL boundaries
    for (int igll = 0; igll < numGll; igll++) {
        int itop = igll >= 1 ? igll - 1 : 0;
        int imid = igll;
        int ibot = igll + 1;
        if (ibot > numGll - 1) {
            ibot = numGll - 1;
        }
        double gll_top = mRlGLL(itop);
        double gll_mid = mRlGLL(imid);
        double gll_bot = mRlGLL(ibot);
        double gll_mid_value = 2. / (gll_top - gll_bot);
        for (int row = 0; row < mRl.rows(); row++) {
            // integrate
            for (int ilayer = colSurf; ilayer < colMoho; ilayer++) {
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
}

void Volumetric3D_crust1::initialize(const std::vector<std::string> &params) {
    try {
        int ipar = 0;
        const std::string source = "Volumetric3D_crust1::initialize";
        Parameters::castValue(mIncludeSediment, params.at(ipar++), source);
        Parameters::castValue(mIncludeIce, params.at(ipar++), source);
        Parameters::castValue(mGeographic, params.at(ipar++), source);
        Parameters::castValue(mRMoho, params.at(ipar++), source); mRMoho *= 1e3;
    } catch (std::out_of_range) {
        // nothing
    }
    
    // delete element boundaries out of crust
    for (int i = 0; i < mElementBoundaries.size(); i++) {
        if (mElementBoundaries[i] < mRMoho - 1.) {
            mElementBoundaries.erase(mElementBoundaries.begin() + i, mElementBoundaries.end());
            break;
        }
    }
    
    // check boundaries
    if (mElementBoundaries.size() < 2) {
        throw std::runtime_error("Volumetric3D_crust1::setupExodusModel || No element layer found in crust.");
    } 
    if (std::abs(mRSurf - mElementBoundaries[0]) > 1.) {
        throw std::runtime_error("Volumetric3D_crust1::setupExodusModel || Conflict in surface radius.");    
    }
    if (std::abs(mRMoho - mElementBoundaries[mElementBoundaries.size() - 1]) > 1.) {
        throw std::runtime_error("Volumetric3D_crust1::setupExodusModel || Conflict in Moho radius.");
    }
        
    // initialize
    initialize();
}

bool Volumetric3D_crust1::get3dProperties(double r, double theta, double phi, double rElemCenter,
    std::vector<MaterialProperty> &properties, 
    std::vector<MaterialRefType> &refTypes,
    std::vector<double> &values) const {
    
    // header
    properties.clear();
    properties.push_back(Volumetric3D::MaterialProperty::VP);
    properties.push_back(Volumetric3D::MaterialProperty::VS);
    properties.push_back(Volumetric3D::MaterialProperty::RHO);
    refTypes = std::vector<MaterialRefType>(3, Volumetric3D::MaterialRefType::Absolute);
    values = std::vector<double>(3, 0.);
    
    // not in crust 
    if (rElemCenter > mRSurf || rElemCenter < mRMoho) {
        return false;
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
    
    // interpolation on sphere
    int llat[2], llon[2];
    double wlat[2], wlon[2];
    XMath::interpLinear(lat, mGridLat, llat[0], wlat[0]);
    XMath::interpLinear(lon, mGridLon, llon[0], wlon[0]);    
    llat[1] = llat[0] + 1;
    llon[1] = llon[0] + 1;
    wlat[1] = 1. - wlat[0];
    wlon[1] = 1. - wlon[1];
    if (llon[1] == sNLon) {
        llon[1] = 0;
    }
    
    // element inside but point slightly outside
    if (r >= mRSurf * 0.999999) r = mRSurf * 0.999999;
    if (r <= mRMoho * 1.000001) r = mRMoho * 1.000001;
    
    double v_p = 0.;
    double v_s = 0.;
    double rho = 0.;
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            double weight = wlat[i] * wlon[j];
            int rowdata = llat[i] * sNLon + llon[j];
            bool found = false;
            for (int iLayer = 1; iLayer < mRlGLL.rows(); iLayer++) {
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
                    v_p += vp * weight;
                    v_s += vs * weight;
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
    
    // when sediment is considered, Vs at some locations may be too small
    v_s = std::max(vsv, 500.);
    v_p = std::max(vpv, sqrt(2) * vsv);
    
    values[0] = v_p;
    values[1] = v_s;
    values[2] = rho;
    return true;
}

std::string Volumetric3D_crust1::verbose() const {
    std::stringstream ss;
    ss << "\n======================= 3D Volumetric ======================" << std::endl;
    ss << "  Model Name            =   Crust 1.0" << std::endl;
    ss << "  Scope                 =   Crust" << std::endl;
    ss << "  Radii (km)            =   [" << mRMoho / 1e3 << ", " << mRSurf / 1e3 << "]" << std::endl;
    ss << "  Reference Type        =   Absolute" << std::endl;
    ss << "  Affected Propertis    =   VP VS RHO" << std::endl;
    ss << "  Include Sediment      =   " << (mIncludeSediment ? "YES" : "NO") << std::endl;
    ss << "  Include Ice           =   " << (mIncludeIce ? "YES" : "NO") << std::endl;
    ss << "  Use Geographic        =   " << (mGeographic ? "YES" : "NO") << std::endl;
    ss << "  Num. Element Layers   =   " << mElementBoundaries.size() - 1 << std::endl;
    ss << "======================= 3D Volumetric ======================\n" << std::endl;
    return ss.str();
}

void Volumetric3D_crust1::setupExodusModel(const ExodusModel *exModel) {
    // find all z-coordinates on the axis
    for (int i = 0; i < exModel->getNumNodes(); i++) {
        double s = exModel->getNodalS(i);
        double z = exModel->getNodalZ(i);
        if (s <= tinyDouble && z > 0) mElementBoundaries.push_back(z);
    }
    std::sort(mElementBoundaries.begin(), mElementBoundaries.end(), std::greater<double>());
}
