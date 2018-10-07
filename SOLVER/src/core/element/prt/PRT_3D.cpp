// PRT_3D.cpp
// created by Kuangdai on 19-May-2017 
// 3D particle relabelling transformation

#include "PRT_3D.h"
#include "SolidElement.h"
#include "FluidElement.h"

#include "SolverFFTW_N3.h"
#include "SolverFFTW_N6.h"
#include "SolverFFTW_N9.h"

PRT_3D::PRT_3D(const RMatXN4 &X) {
    int Nr = X.rows();
    mXFlat0 = X.block(0, 0 * nPE, Nr, nPE);
    mXFlat1 = X.block(0, 1 * nPE, Nr, nPE);
    mXFlat2 = X.block(0, 2 * nPE, Nr, nPE);
    mXFlat3 = X.block(0, 3 * nPE, Nr, nPE);
}

void PRT_3D::sphericalToUndulated(FluidResponse &response) const {
    int Nr = response.mNr;
    RMatXN3 &strain = SolverFFTW_N3::getC2R_RMat(Nr);
    strain.block(0, 0 * nPE, Nr, nPE) = mXFlat0.schur(strain.block(0, 0 * nPE, Nr, nPE)) 
                                      + mXFlat1.schur(strain.block(0, 2 * nPE, Nr, nPE));
    strain.block(0, 1 * nPE, Nr, nPE) = mXFlat0.schur(strain.block(0, 1 * nPE, Nr, nPE)) 
                                      + mXFlat2.schur(strain.block(0, 2 * nPE, Nr, nPE));                                          
    strain.block(0, 2 * nPE, Nr, nPE) = mXFlat3.schur(strain.block(0, 2 * nPE, Nr, nPE));
}

void PRT_3D::undulatedToSpherical(FluidResponse &response) const {
    int Nr = response.mNr;
    RMatXN3 &stress = SolverFFTW_N3::getR2C_RMat(Nr);
    stress.block(0, 2 * nPE, Nr, nPE) = mXFlat1.schur(stress.block(0, 0 * nPE, Nr, nPE))
                                      + mXFlat2.schur(stress.block(0, 1 * nPE, Nr, nPE))
                                      + mXFlat3.schur(stress.block(0, 2 * nPE, Nr, nPE));
    stress.block(0, 0 * nPE, Nr, nPE) = mXFlat0.schur(stress.block(0, 0 * nPE, Nr, nPE)); 
    stress.block(0, 1 * nPE, Nr, nPE) = mXFlat0.schur(stress.block(0, 1 * nPE, Nr, nPE)); 
}

void PRT_3D::sphericalToUndulated(SolidResponse &response) const {
    int Nr = response.mNr;
    const RMatXN9 &sph = SolverFFTW_N9::getC2R_RMat(Nr);
    RMatXN6 &und = SolverFFTW_N6::getC2R_RMat(Nr);
    und.block(0, nPE * 0, Nr, nPE) = mXFlat0.schur(sph.block(0, nPE * 0, Nr, nPE))
                                   + mXFlat1.schur(sph.block(0, nPE * 2, Nr, nPE));
    und.block(0, nPE * 1, Nr, nPE) = mXFlat0.schur(sph.block(0, nPE * 4, Nr, nPE))
                                   + mXFlat2.schur(sph.block(0, nPE * 5, Nr, nPE));
    und.block(0, nPE * 2, Nr, nPE) = mXFlat3.schur(sph.block(0, nPE * 8, Nr, nPE));
    und.block(0, nPE * 3, Nr, nPE) = mXFlat0.schur(sph.block(0, nPE * 7, Nr, nPE))
                                   + mXFlat2.schur(sph.block(0, nPE * 8, Nr, nPE))
                                   + mXFlat3.schur(sph.block(0, nPE * 5, Nr, nPE));
    und.block(0, nPE * 4, Nr, nPE) = mXFlat0.schur(sph.block(0, nPE * 6, Nr, nPE))
                                   + mXFlat1.schur(sph.block(0, nPE * 8, Nr, nPE))
                                   + mXFlat3.schur(sph.block(0, nPE * 2, Nr, nPE));
    und.block(0, nPE * 5, Nr, nPE) = mXFlat0.schur(sph.block(0, nPE * 3, Nr, nPE) 
                                                 + sph.block(0, nPE * 1, Nr, nPE))
                                   + mXFlat1.schur(sph.block(0, nPE * 5, Nr, nPE))
                                   + mXFlat2.schur(sph.block(0, nPE * 2, Nr, nPE));
}

void PRT_3D::undulatedToSpherical(SolidResponse &response) const{
    int Nr = response.mNr;
    const RMatXN6 &und = SolverFFTW_N6::getR2C_RMat(Nr);
    RMatXN9 &sph = SolverFFTW_N9::getR2C_RMat(Nr);
    sph.block(0, nPE * 0, Nr, nPE) = mXFlat0.schur(und.block(0, nPE * 0, Nr, nPE));
    sph.block(0, nPE * 1, Nr, nPE) = mXFlat0.schur(und.block(0, nPE * 5, Nr, nPE));
    sph.block(0, nPE * 2, Nr, nPE) = mXFlat1.schur(und.block(0, nPE * 0, Nr, nPE)) 
                                   + mXFlat3.schur(und.block(0, nPE * 4, Nr, nPE))
                                   + mXFlat2.schur(und.block(0, nPE * 5, Nr, nPE));
    sph.block(0, nPE * 3, Nr, nPE) = sph.block(0, nPE * 1, Nr, nPE);                                  
    sph.block(0, nPE * 4, Nr, nPE) = mXFlat0.schur(und.block(0, nPE * 1, Nr, nPE));
    sph.block(0, nPE * 5, Nr, nPE) = mXFlat2.schur(und.block(0, nPE * 1, Nr, nPE)) 
                                   + mXFlat3.schur(und.block(0, nPE * 3, Nr, nPE))
                                   + mXFlat1.schur(und.block(0, nPE * 5, Nr, nPE));
    sph.block(0, nPE * 6, Nr, nPE) = mXFlat0.schur(und.block(0, nPE * 4, Nr, nPE));
    sph.block(0, nPE * 7, Nr, nPE) = mXFlat0.schur(und.block(0, nPE * 3, Nr, nPE));
    sph.block(0, nPE * 8, Nr, nPE) = mXFlat3.schur(und.block(0, nPE * 2, Nr, nPE)) 
                                   + mXFlat2.schur(und.block(0, nPE * 3, Nr, nPE))
                                   + mXFlat1.schur(und.block(0, nPE * 4, Nr, nPE));
}

void PRT_3D::sphericalToUndulated9(SolidResponse &response) const {
    int Nr = response.mNr;
    RMatXN9 &sph = SolverFFTW_N9::getC2R_RMat(Nr);
    RMatXN9 &und = SolverFFTW_N9::getR2C_RMat(Nr);
    und.block(0, nPE * 0, Nr, nPE) = mXFlat0.schur(sph.block(0, nPE * 0, Nr, nPE))
                                   + mXFlat1.schur(sph.block(0, nPE * 2, Nr, nPE));
    und.block(0, nPE * 1, Nr, nPE) = mXFlat0.schur(sph.block(0, nPE * 3, Nr, nPE))
                                   + mXFlat1.schur(sph.block(0, nPE * 5, Nr, nPE));
    und.block(0, nPE * 2, Nr, nPE) = mXFlat0.schur(sph.block(0, nPE * 6, Nr, nPE))
                                   + mXFlat1.schur(sph.block(0, nPE * 8, Nr, nPE));
    und.block(0, nPE * 3, Nr, nPE) = mXFlat0.schur(sph.block(0, nPE * 1, Nr, nPE))
                                   + mXFlat2.schur(sph.block(0, nPE * 2, Nr, nPE));
    und.block(0, nPE * 4, Nr, nPE) = mXFlat0.schur(sph.block(0, nPE * 4, Nr, nPE))
                                   + mXFlat2.schur(sph.block(0, nPE * 5, Nr, nPE));                                                                 
    und.block(0, nPE * 5, Nr, nPE) = mXFlat0.schur(sph.block(0, nPE * 7, Nr, nPE))
                                   + mXFlat2.schur(sph.block(0, nPE * 8, Nr, nPE));
    und.block(0, nPE * 6, Nr, nPE) = mXFlat3.schur(sph.block(0, nPE * 2, Nr, nPE));
    und.block(0, nPE * 7, Nr, nPE) = mXFlat3.schur(sph.block(0, nPE * 5, Nr, nPE));
    und.block(0, nPE * 8, Nr, nPE) = mXFlat3.schur(sph.block(0, nPE * 8, Nr, nPE));
    // replace
    sph.block(0, 0, Nr, nPE * 9) = und.block(0, 0, Nr, nPE * 9);
}

void PRT_3D::checkCompatibility(int Nr) const{
    int myNr = mXFlat0.rows();
    if (Nr != myNr) {
        throw std::runtime_error("PRT_3D::checkCompatibility || Incompatible size.");
    }
}
