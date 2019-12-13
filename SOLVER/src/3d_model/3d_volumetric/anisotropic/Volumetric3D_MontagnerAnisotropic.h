// Volumetric3D_MontagnerAnisotropic.h
// created by Kuangdai on 11-Dec-2019 
// full anisotropic model in the upper mantle by Jean-Paul Montagner

#pragma once
#include "Volumetric3D.h"

class Volumetric3D_MontagnerAnisotropic: public Volumetric3D {
public:
    
    void initialize();
    
    bool get3dProperties(double r, double theta, double phi, double rElemCenter,
        std::vector<MaterialProperty> &properties, 
        std::vector<MaterialRefType> &refTypes,
        std::vector<double> &values) const;
    
    std::string verbose() const;
    
private:
    // dimensions
    static const int mN_par = 14;
    static const int mN_dep = 34;
    static const int mN_lat = 37;
    static const int mN_lon = 73;
    
    // anchoring depths
    const double mDepths[mN_dep + 1] = {
        -1., 670., 647., 623., 600., 578., 556., 533., 511., 489., 467., 444.,
       422., 400., 400., 370., 340., 310., 280., 250., 220., 220., 200.,
       180., 160., 140., 120., 100.,  80.,  80.,  69.,  58.,  47.,  35., 24.};
       
    // beta array, read from file  
    double mBeta[mN_par+1][mN_dep+1][mN_lat+1][mN_lon+1];
};
