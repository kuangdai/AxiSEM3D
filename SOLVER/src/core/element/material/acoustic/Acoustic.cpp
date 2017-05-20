// Acoustic.cpp
// created by Kuangdai on 23-Apr-2016 
// base class of acoustic constitutive relation

#include "Acoustic1D.h"
#include "Acoustic3D.h"

Acoustic *Acoustic::createAcoustic(const RDMatXN &KFluid, bool elem1D) {
    if (elem1D) {
        return new Acoustic1D(KFluid);
    } else {
        return new Acoustic3D(KFluid);
    }
}

