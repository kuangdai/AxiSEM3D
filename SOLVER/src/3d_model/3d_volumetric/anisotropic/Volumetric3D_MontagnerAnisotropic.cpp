// Volumetric3D_MontagnerAnisotropic.h
// created by Kuangdai on 11-Dec-2019 
// full anisotropic model in the upper mantle by Jean-Paul Montagner

#include "Volumetric3D_MontagnerAnisotropic.h"
#include <fstream>
#include "XMPI.h"

void Volumetric3D_MontagnerAnisotropic::initialize() {
    if (XMPI::root()) {
        // open beta_values.txt
        std::string fname = projectDirectory + 
        "/src/3d_model/3d_volumetric/anisotropic/beta_values.txt";
        std::fstream fs(fname, std::fstream::in);
        if (!fs) {
            throw std::runtime_error("Volumetric3D_MontagnerAnisotropic::initialize || "
            "Error opening model data file: || " + fname);
        }
        
        // read comments
        std::string line;
        getline(fs, line);
        getline(fs, line);
        getline(fs, line);
        
        // read beta, note that index starts from 1
        for (int ipar = 1; ipar <= mN_par; ipar++) {
            for (int idep = 1; idep <= mN_dep; idep++) {
                for (int ilat = 1; ilat <= mN_lat; ilat++) {
                    for (int ilon = 1; ilon <= mN_lon; ilon++) {
                        fs >> mBeta[ipar][idep][ilat][ilon];
                    }    
                }    
            }
        }
    }
    // broadcast
    XMPI::bcast(&(mBeta[0][0][0][0]), 
    (mN_par + 1) * (mN_dep + 1) * (mN_lat + 1) * (mN_lon + 1));
    
    
    // // testing code
    // std::vector<MaterialProperty> properties;
    // std::vector<MaterialRefType> refTypes;
    // std::vector<double> values;
    // get3dProperties(6000e3, 48*degree, 155*degree, 6000e3, properties, refTypes, values);
    // for (int i=0; i<22; i++) {
    //     XMPI::cout << values[i] << XMPI::endl;
    // }
    // exit(0);
}

bool Volumetric3D_MontagnerAnisotropic::get3dProperties(
    double r, double theta, double phi, double rElemCenter,
    std::vector<MaterialProperty> &properties, 
    std::vector<MaterialRefType> &refTypes,
    std::vector<double> &values) const {
        
    // output header
    properties.push_back(Volumetric3D::MaterialProperty::RHO);
    properties.push_back(Volumetric3D::MaterialProperty::C11);
    properties.push_back(Volumetric3D::MaterialProperty::C12);
    properties.push_back(Volumetric3D::MaterialProperty::C13);
    properties.push_back(Volumetric3D::MaterialProperty::C14);
    properties.push_back(Volumetric3D::MaterialProperty::C15);
    properties.push_back(Volumetric3D::MaterialProperty::C16);
    properties.push_back(Volumetric3D::MaterialProperty::C22);
    properties.push_back(Volumetric3D::MaterialProperty::C23);
    properties.push_back(Volumetric3D::MaterialProperty::C24);
    properties.push_back(Volumetric3D::MaterialProperty::C25);
    properties.push_back(Volumetric3D::MaterialProperty::C26);
    properties.push_back(Volumetric3D::MaterialProperty::C33);
    properties.push_back(Volumetric3D::MaterialProperty::C34);
    properties.push_back(Volumetric3D::MaterialProperty::C35);
    properties.push_back(Volumetric3D::MaterialProperty::C36);
    properties.push_back(Volumetric3D::MaterialProperty::C44);
    properties.push_back(Volumetric3D::MaterialProperty::C45);
    properties.push_back(Volumetric3D::MaterialProperty::C46);
    properties.push_back(Volumetric3D::MaterialProperty::C55);
    properties.push_back(Volumetric3D::MaterialProperty::C56);
    properties.push_back(Volumetric3D::MaterialProperty::C66);
    refTypes = std::vector<MaterialRefType>(22, MaterialRefType::Absolute);
    
    
    //////////////////////////////////////////
    // the rest is traslated from specfem3d_globe
    //////////////////////////////////////////
        
        
    int ndepth = mN_dep;
    double pxy0 = 5.;
    double x0 = 0.;
    double y0 = 0.;
    
    int nx0 = mN_lat;
    int ny0 = mN_lon;
    int nz0 = mN_dep;
    
    double tet = theta / degree;
    double ph = phi / degree;
    
    // avoid edge effects
    if (tet == 0.0) tet = 0.000001;
    if (tet == 180.) tet = 0.999999*tet;
    if (ph == 0.0) ph = 0.000001;
    if (ph == 360.) ph = 0.999999*ph;
    
    // depth
    double depth = 6371. - r / 1e3;
    if (depth <= mDepths[nz0] || depth >= mDepths[1]) {
        return false;
    } 
    
    int icolat = int(int(tet + pxy0)/pxy0);
    int ilon = int(int(ph + pxy0)/pxy0);
    
    int icz0 = 0;
    for (int idep = 1; idep <= ndepth; idep++) {
        if (mDepths[idep] > depth) icz0 = icz0 + 1;
    }
    
    int ict0 = icolat;
    int ict1 = ict0 + 1;
    int icp0 = ilon;
    int icp1 = icp0 + 1;
    int icz1 = icz0 + 1;
    
    // check that parameters make sense
    if (ict0 < 1 || ict0 > nx0) return false;
    if (ict1 < 1 || ict1 > nx0) return false;
    if (icp0 < 1 || icp0 > ny0) return false;
    if (icp1 < 1 || icp1 > ny0) return false;
    if (icz0 < 1 || icz0 > nz0) return false;
    if (icz1 < 1 || icz1 > nz0) return false;
    
    double anispara[mN_par + 1][3][5];
    for (int ipar = 1; ipar <= mN_par; ipar++) {
        anispara[ipar][1][1] = mBeta[ipar][icz0][ict0][icp0];
        anispara[ipar][2][1] = mBeta[ipar][icz1][ict0][icp0];
        anispara[ipar][1][2] = mBeta[ipar][icz0][ict0][icp1];
        anispara[ipar][2][2] = mBeta[ipar][icz1][ict0][icp1];
        anispara[ipar][1][3] = mBeta[ipar][icz0][ict1][icp0];
        anispara[ipar][2][3] = mBeta[ipar][icz1][ict1][icp0];
        anispara[ipar][1][4] = mBeta[ipar][icz0][ict1][icp1];
        anispara[ipar][2][4] = mBeta[ipar][icz1][ict1][icp1];
    }
    
    double tei = pxy0*ict0 + x0 - pxy0;
    double fi = pxy0*icp0 + y0 - pxy0;
    
    double d1 = sqrt(((tei - tet)*(tei - tet)) + ((fi - ph)*(fi - ph))*(sin((tet + tei)*degree/2.)*sin((tet + tei)*degree/2.)));
    double d2 = sqrt(((tei - tet + pxy0)*(tei - tet + pxy0)) + ((fi - ph)*(fi - ph))*(sin((tet + tei + pxy0)*degree/2.)*sin((tet + tei + pxy0)*degree/2.)));
    double d3 = sqrt(((tei - tet)*(tei - tet)) + ((fi - ph + pxy0)*(fi - ph + pxy0))*(sin((tet + tei)*degree/2.)*sin((tet + tei)*degree/2.)));
    double d4 = sqrt(((tei - tet + pxy0)*(tei - tet + pxy0)) + ((fi - ph + pxy0)*(fi - ph + pxy0))*(sin((tet + tei + pxy0)*degree/2.)*sin((tet + tei + pxy0)*degree/2.)));
    
    double sd = d2*d3*d4 + d1*d2*d4 + d1*d3*d4 + d1*d2*d3;
    double thickness = mDepths[icz0] - mDepths[icz1];
    double dprof1 = mDepths[icz0] - depth;
    double dprof2 = depth - mDepths[icz1];
    double eps = 0.01;
    
    double elpar[mN_par+1];
    for (int ipar = 1; ipar <= mN_par; ipar++) {
        double pc1, pc2, pc3, pc4;
        if (thickness < eps) {
            pc1 = anispara[ipar][1][1];
            pc2 = anispara[ipar][1][2];
            pc3 = anispara[ipar][1][3];
            pc4 = anispara[ipar][1][4];
        } else {
            double dpr1 = dprof1/thickness;
            double dpr2 = dprof2/thickness;
            pc1 = anispara[ipar][1][1]*dpr2+anispara[ipar][2][1]*dpr1;
            pc2 = anispara[ipar][1][2]*dpr2+anispara[ipar][2][2]*dpr1;
            pc3 = anispara[ipar][1][3]*dpr2+anispara[ipar][2][3]*dpr1;
            pc4 = anispara[ipar][1][4]*dpr2+anispara[ipar][2][4]*dpr1;
        }
        double param = pc1*d2*d3*d4 + pc2*d1*d3*d4 + pc3*d1*d2*d4 + pc4*d1*d2*d3;
        param = param/sd;
        elpar[ipar] = param;
    }
    
    double rho = elpar[1];
    double A = elpar[2];
    double C = elpar[3];
    double F = elpar[4];
    double AL = elpar[5];
    double AN = elpar[6];
    double BC = elpar[7];
    double BS = elpar[8];
    double GC = elpar[9];
    double GS = elpar[10];
    double HC = elpar[11];
    double HS = elpar[12];
    double EC = elpar[13];
    double ES = elpar[14];
    double C1p = 0.0;
    double S1p = 0.0;
    double C1sv = 0.0;
    double S1sv = 0.0;
    double C1sh = 0.0;
    double S1sh = 0.0;
    double C3 = 0.0;
    double S3 = 0.0;
    
    double d11 = A + EC + BC;
    double d12 = A - 2.*AN - EC;
    double d13 = F + HC;
    double d14 = S3 + 2.*S1sh + 2.*S1p;
    double d15 = 2.*C1p + C3;
    double d16 = -BS/2. - ES;
    double d22 = A + EC - BC;
    double d23 = F - HC;
    double d24 = 2.*S1p - S3;
    double d25 = 2.*C1p - 2.*C1sh - C3;
    double d26 = -BS/2. + ES;
    double d33 = C;
    double d34 = 2.*(S1p - S1sv);
    double d35 = 2.*(C1p - C1sv);
    double d36 = -HS;
    double d44 = AL - GC;
    double d45 = -GS;
    double d46 = C1sh - C3;
    double d55 = AL + GC;
    double d56 = S3 - S1sh;
    double d66 = AN - EC;
    
    // copy to output
    values = {rho,
        d11, d12, d13, d14, d15, d16,
        d22, d23, d24, d25, d26,
        d33, d34, d35, d36,
        d44, d45, d46,
        d55, d56,
        d66
    };
    
    // SI units
    values[0] *= 1e3;
    for (int i = 1; i < 22; i++) {
        values[i] *= 1e9;
    }
    return true;
}

std::string Volumetric3D_MontagnerAnisotropic::verbose() const {
    std::stringstream ss;
    ss << "\n======================= 3D Volumetric ======================" << std::endl;
    ss << "  Model Name           =   montagner_anisotropic" << std::endl;
    ss << "  Scope                =   Upper Mantle" << std::endl;
    ss << "  Radii (km)           =   [" << 6371-mDepths[1] << ", " << 6371-mDepths[mN_dep] << "]" << std::endl;
    ss << "  Reference Type       =   Absolute" << std::endl;
    ss << "  Affected Propertis   =   Cijkl and RHO" << std::endl;
    ss << "======================= 3D Volumetric ======================\n" << std::endl;
    return ss.str();
}
    
