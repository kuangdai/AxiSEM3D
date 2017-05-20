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
        return new PRT_1D(X);
    } else {
        return new PRT_3D(X);
    }
}

