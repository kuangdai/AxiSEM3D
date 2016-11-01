// Volumetric3D_s40rts.cpp
// created by Kuangdai on 16-May-2016 
// mantle model s40rts
// http://www.earth.lsa.umich.edu/~jritsema/research.html

#include "Volumetric3D_s40rts.h"
#include <sstream>

extern "C" {
    void __s40rts_MOD_initialize_s40rts(double *rcmb, double *rmoho, double *rearth);
    void __s40rts_MOD_finalize_s40rts();
    bool __s40rts_MOD_perturb_s40rts(double *r, double *theta, double *phi, double *r_center,
        double *vp, double *vs);
};


void Volumetric3D_s40rts::initialize() {
    __s40rts_MOD_initialize_s40rts(&mRCMB, &mRMoho, &mRSurf);
}

void Volumetric3D_s40rts::initialize(const std::vector<double> &params) {
    try {
        int ipar = 0;
        mScaleRho = params.at(ipar++);
        mRCMB = params.at(ipar++) * 1e3;
        mRMoho = params.at(ipar++) * 1e3;
        mRSurf = params.at(ipar++) * 1e3;
    } catch (std::out_of_range) {
        // nothing
    }
    initialize();
}

void Volumetric3D_s40rts::finalize() {
    __s40rts_MOD_finalize_s40rts();
}

bool Volumetric3D_s40rts::get3dProperties(double r, double theta, double phi, double rElemCenter,
    double &dvpv, double &dvph, double &dvsv, double &dvsh, double &drho) const {
    double dvp, dvs;
    bool result = __s40rts_MOD_perturb_s40rts(&r, &theta, &phi, &rElemCenter, &dvp, &dvs);
    dvpv = dvph = dvp;
    dvsv = dvsh = dvs;
    drho = mScaleRho * dvs;
    return result;
}

std::string Volumetric3D_s40rts::verbose() const {
    std::stringstream ss;
    ss << "\n======================= 3D Volumetric ======================" << std::endl;
    ss << "  Model Name           =   s40rts" << std::endl;
    ss << "  Scope                =   Mantle" << std::endl;
    ss << "  Radii (km)           =   [" << mRCMB / 1e3 << ", " << mRMoho / 1e3 << "]" << std::endl;
    ss << "  Reference Type       =   Reference1D" << std::endl;
    ss << "  Max. Fourier Order   =   40" << std::endl;
    ss << "  Anisotropic          =   NO" << std::endl;
    ss << "  3D Density           =   " << (mScaleRho != 0. ? "YES" : "NO") << std::endl;
    ss << "======================= 3D Volumetric ======================\n" << std::endl;
    return ss.str();
}

