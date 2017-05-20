// PRT_1D.cpp
// created by Kuangdai on 19-May-2017 
// 1D particle relabelling transformation 

#include "PRT_1D.h"
#include "SolidElement.h"
#include "FluidElement.h"

#include "SolverFFTW_N3.h"
#include "SolverFFTW_N6.h"
#include "SolverFFTW_N9.h"

PRT_1D::PRT_1D(const RDMatXN4 &X) {
    for (int idim = 0; idim < 4; idim++) {
        for (int ipol = 0; ipol < nPntEdge; ipol++) {
            mXStruct[idim].block(ipol, 0, 1, nPntEdge) 
                = X.block(0, nPE * idim + nPntEdge * ipol, 1, nPntEdge).cast<Real>();
        }
    }
}

void PRT_1D::sphericalToUndulated(FluidElementResponse &response) const {
    for (int alpha = 0; alpha < response.mNu; alpha++) {
        ar3_CMatPP &strain = response.mStrain[alpha];
        strain[0] = mXStruct[0].schur(strain[0]) 
                  + mXStruct[1].schur(strain[2]);
        strain[1] = mXStruct[0].schur(strain[1]) 
                  + mXStruct[2].schur(strain[2]);
        strain[2] = mXStruct[3].schur(strain[2]);
    }
}

void PRT_1D::undulatedToSpherical(FluidElementResponse &response) const {
    for (int alpha = 0; alpha < response.mNu; alpha++) {
        ar3_CMatPP &stress = response.mStress[alpha];
        stress[2] = mXStruct[1].schur(stress[0]) 
                  + mXStruct[2].schur(stress[1]) 
                  + mXStruct[3].schur(stress[2]);
        stress[0] = mXStruct[0].schur(stress[0]);
        stress[1] = mXStruct[0].schur(stress[1]);
    }
}

void PRT_1D::sphericalToUndulated(SolidElementResponse &response) const {
    for (int alpha = 0; alpha < response.mNu; alpha++) {
        const ar9_CMatPP &sph = response.mStrainC9[alpha];
        ar6_CMatPP &und = response.mStrainC6[alpha];
        und[0] = mXStruct[0].schur(sph[0])
               + mXStruct[1].schur(sph[2]);
        und[1] = mXStruct[0].schur(sph[4])
               + mXStruct[2].schur(sph[5]);
        und[2] = mXStruct[3].schur(sph[8]);
        und[3] = mXStruct[0].schur(sph[7])
               + mXStruct[2].schur(sph[8])
               + mXStruct[3].schur(sph[5]);
        und[4] = mXStruct[0].schur(sph[6])
               + mXStruct[1].schur(sph[8])
               + mXStruct[3].schur(sph[2]);
        und[5] = mXStruct[0].schur(sph[3] 
                                 + sph[1])
               + mXStruct[1].schur(sph[5])
               + mXStruct[2].schur(sph[2]);
    }
}

void PRT_1D::undulatedToSpherical(SolidElementResponse &response) const{
    for (int alpha = 0; alpha < response.mNu; alpha++) {
        const ar6_CMatPP &und = response.mStressC6[alpha];
        ar9_CMatPP &sph = response.mStressC9[alpha];
        sph[0] = mXStruct[0].schur(und[0]);
        sph[1] = mXStruct[0].schur(und[5]);
        sph[2] = mXStruct[1].schur(und[0]) 
               + mXStruct[3].schur(und[4])
               + mXStruct[2].schur(und[5]);
        sph[3] = sph[1];                                  
        sph[4] = mXStruct[0].schur(und[1]);
        sph[5] = mXStruct[2].schur(und[1]) 
               + mXStruct[3].schur(und[3])
               + mXStruct[1].schur(und[5]);
        sph[6] = mXStruct[0].schur(und[4]);
        sph[7] = mXStruct[0].schur(und[3]);
        sph[8] = mXStruct[3].schur(und[2]) 
               + mXStruct[2].schur(und[3])
               + mXStruct[1].schur(und[4]);
    }
}
