// Attenuation3D.cpp
// created by Kuangdai on 29-Apr-2016 
// base class of 3D attenuation 

#include "Attenuation3D.h"

Attenuation3D::Attenuation3D(int nsls, const RColX &alpha, 
    const RColX &beta, const RColX &gamma):
Attenuation(nsls, alpha, beta, gamma) {
    // nothing
}

