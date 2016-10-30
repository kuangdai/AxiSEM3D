// Attenuation1D.cpp
// created by Kuangdai on 29-Apr-2016 
// base class of 1D attenuation 

#include "Attenuation1D.h"

Attenuation1D::Attenuation1D(int nsls, const RColX &alpha, 
    const RColX &beta, const RColX &gamma):
Attenuation(nsls, alpha, beta, gamma) {
    // nothing
}

