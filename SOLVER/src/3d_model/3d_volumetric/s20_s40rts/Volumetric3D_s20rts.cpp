// Volumetric3D_s20rts.cpp
// created by Kuangdai on 16-May-2016 
// mantle model s20rts
// http://www.earth.lsa.umich.edu/~jritsema/research.html

#include "Volumetric3D_s20rts.h"
#include <sstream>

extern "C" {
    void __s20rts_MOD_initialize_s20rts(double *rcmb, double *rmoho, double *rearth);
    void __s20rts_MOD_finalize_s20rts();
    bool __s20rts_MOD_perturb_s20rts(double *r, double *theta, double *phi, double *r_center,
        double *vp, double *vs);
};


void Volumetric3D_s20rts::initialize() {
    __s20rts_MOD_initialize_s20rts(&mRCMB, &mRMoho, &mRSurf);
}

void Volumetric3D_s20rts::initialize(const std::vector<double> &params) {
    if (params.size() >= 1) mScaleRho = params[0];
    if (params.size() >= 4) {
        mRCMB = params[1];
        mRMoho = params[2];
        mRSurf = params[3];
    }
    initialize();
}

void Volumetric3D_s20rts::finalize() {
    __s20rts_MOD_finalize_s20rts();
}

bool Volumetric3D_s20rts::get3dProperties(double r, double theta, double phi, double rElemCenter,
    double &dvpv, double &dvph, double &dvsv, double &dvsh, double &drho) const {
    double dvp, dvs;
    bool result = __s20rts_MOD_perturb_s20rts(&r, &theta, &phi, &rElemCenter, &dvp, &dvs);
    dvpv = dvph = dvp;
    dvsv = dvsh = dvs;
    drho = mScaleRho * dvs;
    return result;
}

std::string Volumetric3D_s20rts::verbose() const {
    std::stringstream ss;
    ss << "\n======================= 3D Volumetric ======================" << std::endl;
    ss << "  Model Name           =   s20rts" << std::endl;
    ss << "  Scope                =   Mantle" << std::endl;
    ss << "  Radii (km)           =   [" << mRCMB / 1e3 << ", " << mRMoho / 1e3 << "]" << std::endl;
    ss << "  Reference Type       =   Reference1D" << std::endl;
    ss << "  Max. Fourier Order   =   20" << std::endl;
    ss << "  Anisotropic          =   NO" << std::endl;
    ss << "  3D Density           =   " << (mScaleRho != 0. ? "YES" : "NO") << std::endl;
    ss << "======================= 3D Volumetric ======================\n" << std::endl;
    return ss.str();
}

