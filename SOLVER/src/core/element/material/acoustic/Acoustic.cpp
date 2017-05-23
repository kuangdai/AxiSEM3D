// Acoustic.cpp
// created by Kuangdai on 23-Apr-2016 
// base class of acoustic constitutive relation

#include "Acoustic1D.h"
#include "Acoustic3D.h"

Acoustic *Acoustic::createAcoustic(const RDMatXN &KFluid, bool elem1D) {
    if (elem1D) {
        RMatPP kstruct;
        for (int ipol = 0; ipol < nPntEdge; ipol++) {
            kstruct.block(ipol, 0, 1, nPntEdge) = 
            KFluid.block(0, nPntEdge * ipol, 1, nPntEdge).cast<Real>();
        }
        return new Acoustic1D(kstruct);
    } else {
        return new Acoustic3D(KFluid.cast<Real>());
    }
}

