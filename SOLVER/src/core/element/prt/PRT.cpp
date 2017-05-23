// PRT.cpp
// created by Kuangdai on 19-May-2017 
// particle relabelling transformation

#include "PRT.h"
#include "PRT_1D.h"
#include "PRT_3D.h"

PRT *PRT::createPRT(const RDMatXN4 &X, bool elem1D) {
    if (X.array().abs().maxCoeff() < tinyDouble) {
        return 0;
    }
    
    if (elem1D) {
        std::array<RMatPP, 4> xstruct;
        for (int idim = 0; idim < 4; idim++) {
            for (int ipol = 0; ipol < nPntEdge; ipol++) {
                xstruct[idim].block(ipol, 0, 1, nPntEdge) = 
                X.block(0, nPE * idim + nPntEdge * ipol, 1, nPntEdge).cast<Real>();
            }
        }
        return new PRT_1D(xstruct);
    } else {
        return new PRT_3D(X.cast<Real>());
    }
}

