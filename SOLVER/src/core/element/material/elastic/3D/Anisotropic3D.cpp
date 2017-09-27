// Anisotropic3D.cpp
// created by Kuangdai on 25-Sep-2017 
// full anisotropic 3D material 

#include "Anisotropic3D.h"
#include "Attenuation3D.h"
#include "SolidElement.h"
#include "SolverFFTW_N6.h"

void Anisotropic3D::strainToStress(SolidResponse &response) const {
    int Nr = response.mNr;
    const RMatXN6 &strainTIsoR = SolverFFTW_N6::getC2R_RMat(Nr);
    RMatXN6 &stressTIsoR = SolverFFTW_N6::getR2C_RMat(Nr);
    stressTIsoR.block(0, 0 * nPE, Nr, nPE) = mC11.schur(strainTIsoR.block(0, 0 * nPE, Nr, nPE))
                                           + mC12.schur(strainTIsoR.block(0, 1 * nPE, Nr, nPE))
                                           + mC13.schur(strainTIsoR.block(0, 2 * nPE, Nr, nPE))
                                           + mC14.schur(strainTIsoR.block(0, 3 * nPE, Nr, nPE))
                                           + mC15.schur(strainTIsoR.block(0, 4 * nPE, Nr, nPE))
                                           + mC16.schur(strainTIsoR.block(0, 5 * nPE, Nr, nPE));
    stressTIsoR.block(0, 1 * nPE, Nr, nPE) = mC12.schur(strainTIsoR.block(0, 0 * nPE, Nr, nPE))
                                           + mC22.schur(strainTIsoR.block(0, 1 * nPE, Nr, nPE))
                                           + mC23.schur(strainTIsoR.block(0, 2 * nPE, Nr, nPE))
                                           + mC24.schur(strainTIsoR.block(0, 3 * nPE, Nr, nPE))
                                           + mC25.schur(strainTIsoR.block(0, 4 * nPE, Nr, nPE))
                                           + mC26.schur(strainTIsoR.block(0, 5 * nPE, Nr, nPE));    
    stressTIsoR.block(0, 2 * nPE, Nr, nPE) = mC13.schur(strainTIsoR.block(0, 0 * nPE, Nr, nPE))
                                           + mC23.schur(strainTIsoR.block(0, 1 * nPE, Nr, nPE))
                                           + mC33.schur(strainTIsoR.block(0, 2 * nPE, Nr, nPE))
                                           + mC34.schur(strainTIsoR.block(0, 3 * nPE, Nr, nPE))
                                           + mC35.schur(strainTIsoR.block(0, 4 * nPE, Nr, nPE))
                                           + mC36.schur(strainTIsoR.block(0, 5 * nPE, Nr, nPE));
    stressTIsoR.block(0, 3 * nPE, Nr, nPE) = mC14.schur(strainTIsoR.block(0, 0 * nPE, Nr, nPE))
                                           + mC24.schur(strainTIsoR.block(0, 1 * nPE, Nr, nPE))
                                           + mC34.schur(strainTIsoR.block(0, 2 * nPE, Nr, nPE))
                                           + mC44.schur(strainTIsoR.block(0, 3 * nPE, Nr, nPE))
                                           + mC45.schur(strainTIsoR.block(0, 4 * nPE, Nr, nPE))
                                           + mC46.schur(strainTIsoR.block(0, 5 * nPE, Nr, nPE));
    stressTIsoR.block(0, 4 * nPE, Nr, nPE) = mC15.schur(strainTIsoR.block(0, 0 * nPE, Nr, nPE))
                                           + mC25.schur(strainTIsoR.block(0, 1 * nPE, Nr, nPE))
                                           + mC35.schur(strainTIsoR.block(0, 2 * nPE, Nr, nPE))
                                           + mC45.schur(strainTIsoR.block(0, 3 * nPE, Nr, nPE))
                                           + mC55.schur(strainTIsoR.block(0, 4 * nPE, Nr, nPE))
                                           + mC56.schur(strainTIsoR.block(0, 5 * nPE, Nr, nPE));  
    stressTIsoR.block(0, 5 * nPE, Nr, nPE) = mC16.schur(strainTIsoR.block(0, 0 * nPE, Nr, nPE))
                                           + mC26.schur(strainTIsoR.block(0, 1 * nPE, Nr, nPE))
                                           + mC36.schur(strainTIsoR.block(0, 2 * nPE, Nr, nPE))
                                           + mC46.schur(strainTIsoR.block(0, 3 * nPE, Nr, nPE))
                                           + mC56.schur(strainTIsoR.block(0, 4 * nPE, Nr, nPE))
                                           + mC66.schur(strainTIsoR.block(0, 5 * nPE, Nr, nPE));                                                                                                                     
    if (mAttenuation) {
        mAttenuation->applyToStress(stressTIsoR);
        mAttenuation->updateMemoryVariables(strainTIsoR);
    }
}

void Anisotropic3D::checkCompatibility(int Nr) const {
    Elastic3D::checkCompatibility(Nr);
    int myNr = mC11.rows();
    if (Nr != myNr) {
        throw std::runtime_error("Anisotropic3D::checkCompatibility || Incompatible size.");
    }
}

