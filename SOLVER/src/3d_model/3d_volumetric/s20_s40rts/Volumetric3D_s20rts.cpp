// Volumetric3D_s20rts.cpp
// created by Kuangdai on 16-May-2016 
// mantle model s20rts
// http://www.earth.lsa.umich.edu/~jritsema/research.html

#include "Volumetric3D_s20rts.h"
#include <sstream>
#include <fstream>
#include "XMPI.h"
#include "XMath.h"

extern "C" {
    void __s20rts_MOD_initialize_s20rts(double *rcmb, double *rmoho, double *rearth, 
        double meta_data_p12[], double meta_data_s20[]);
    void __s20rts_MOD_finalize_s20rts();
    bool __s20rts_MOD_perturb_s20rts(double *r, double *theta, double *phi, double *r_center,
        double *vp, double *vs);
};


void Volumetric3D_s20rts::initialize() {
    const int np12 = 3549;
    const int ns20 = 9261;
    double meta_data_p12[np12], meta_data_s20[ns20];
    if (XMPI::root()) {
        std::string path = projectDirectory + "/src/3d_model/3d_volumetric/s20_s40rts/data";
        std::fstream fs;
        fs.open(path + "/P12.dat", std::fstream::in);
        for (int i = 0; i < np12; i++) fs >> meta_data_p12[i];
        fs.close();
        fs.open(path + "/S20RTS.dat", std::fstream::in);
        for (int i = 0; i < ns20; i++) fs >> meta_data_s20[i];
        fs.close();
    }
    XMPI::bcast(meta_data_p12, np12);
    XMPI::bcast(meta_data_s20, ns20);
    __s20rts_MOD_initialize_s20rts(&mRCMB, &mRMoho, &mRSurf, meta_data_p12, meta_data_s20);
}

void Volumetric3D_s20rts::initialize(const std::vector<std::string> &params) {
    try {
        int ipar = 0;
        const std::string source = "Volumetric3D_s20rts::initialize";
        XMath::castValue(mScaleRho, params.at(ipar++), source);
        XMath::castValue(mRCMB, params.at(ipar++), source); mRCMB *= 1e3;
        XMath::castValue(mRMoho, params.at(ipar++), source); mRMoho *= 1e3;
        XMath::castValue(mRSurf, params.at(ipar++), source); mRSurf *= 1e3;
    } catch (std::out_of_range) {
        // nothing
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

