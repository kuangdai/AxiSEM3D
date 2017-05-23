// TransverselyIsotropic3D.cpp
// created by Kuangdai on 22-Apr-2016 
// transversely isotropic 3D material

#include "TransverselyIsotropic3D.h"
#include "Attenuation3D.h"
#include "SolidElement.h"
#include "SolverFFTW_N6.h"

void TransverselyIsotropic3D::strainToStress(SolidResponse &response) const {
    int Nr = response.mNr;
    // to avoid dynamic allocation, use stressTIsoR.block(0, 3&4 * nPE, Nr, nPE) as temp memory
    const RMatXN6 &strainTIsoR = SolverFFTW_N6::getC2R_RMat(Nr);
    RMatXN6 &stressTIsoR = SolverFFTW_N6::getR2C_RMat(Nr);
    stressTIsoR.block(0, 3 * nPE, Nr, nPE) = strainTIsoR.block(0, 0 * nPE, Nr, nPE) + strainTIsoR.block(0, 1 * nPE, Nr, nPE);
    stressTIsoR.block(0, 4 * nPE, Nr, nPE) = mA.schur(stressTIsoR.block(0, 3 * nPE, Nr, nPE)) + mF.schur(strainTIsoR.block(0, 2 * nPE, Nr, nPE));
    
    stressTIsoR.block(0, 0 * nPE, Nr, nPE) = stressTIsoR.block(0, 4 * nPE, Nr, nPE) - mN2.schur(strainTIsoR.block(0, 1 * nPE, Nr, nPE));
    stressTIsoR.block(0, 1 * nPE, Nr, nPE) = stressTIsoR.block(0, 4 * nPE, Nr, nPE) - mN2.schur(strainTIsoR.block(0, 0 * nPE, Nr, nPE));
    stressTIsoR.block(0, 2 * nPE, Nr, nPE) = mC.schur(strainTIsoR.block(0, 2 * nPE, Nr, nPE)) + mF.schur(stressTIsoR.block(0, 3 * nPE, Nr, nPE));
    stressTIsoR.block(0, 3 * nPE, Nr, nPE) = mL.schur(strainTIsoR.block(0, 3 * nPE, Nr, nPE));
    stressTIsoR.block(0, 4 * nPE, Nr, nPE) = mL.schur(strainTIsoR.block(0, 4 * nPE, Nr, nPE));
    stressTIsoR.block(0, 5 * nPE, Nr, nPE) = mN.schur(strainTIsoR.block(0, 5 * nPE, Nr, nPE));
    if (mAttenuation) {
        mAttenuation->applyToStress(stressTIsoR);
        mAttenuation->updateMemoryVariables(strainTIsoR);
    }
}

void TransverselyIsotropic3D::checkCompatibility(int Nr) const {
    Elastic3D::checkCompatibility(Nr);
    int myNr = mA.rows();
    if (Nr != myNr) {
        throw std::runtime_error("TransverselyIsotropic3D::checkCompatibility || Incompatible size.");
    }
}

