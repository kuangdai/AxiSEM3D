// FluidElement.cpp
// created by Kuangdai on 29-Mar-2016 
// base class of fluid elements

#include "FluidElement.h"
#include "Point.h"
#include "Gradient.h"
#include "PRT.h"
#include "Acoustic.h"
#include "CrdTransTIsoFluid.h"
#include "FieldFFT.h"

#include "MultilevelTimer.h"

FluidElement::FluidElement(Gradient *grad, PRT *prt,
    const std::array<Point *, nPntElem> &points, 
    Acoustic *acous): 
Element(grad, prt, points), mAcoustic(acous), mCrdTransTIso(0) {
    mAcoustic->checkCompatibility(mMaxNr);
    // TISO
    mInTIso = mHasPRT;
    if (mInTIso) {
        mCrdTransTIso = new CrdTransTIsoFluid(formThetaMat());
    } 
    // 3D
    bool acous1D = mAcoustic->is1D();
    if (mHasPRT) {
        if (acous1D != mPRT->is1D()) {
            throw std::runtime_error("FluidElement::FluidElement || "
                "Particle Relabelling and Elasticity are generated in different spaces.");  
        }
    } 
    mElem3D = !acous1D;
}

FluidElement::~FluidElement() {
    delete mAcoustic;
    if (mInTIso) {
        delete mCrdTransTIso;
    }
}
    
void FluidElement::computeStiff() const {
    // setup static
    sResponse.setNr(mMaxNr);
    
    // get displ from points
    int ipnt = 0;
    for (int ipol = 0; ipol <= nPol; ipol++) {
        for (int jpol = 0; jpol <= nPol; jpol++) {
            mPoints[ipnt++]->scatterDisplToElement(sResponse.mDispl, ipol, jpol, mMaxNu);
        }
    }
        
    // compute stiff
    displToStiff();
    
    // set stiff to points
    ipnt = 0;
    for (int ipol = 0; ipol <= nPol; ipol++) {
        for (int jpol = 0; jpol <= nPol; jpol++) {
            mPoints[ipnt++]->gatherStiffFromElement(sResponse.mStiff, ipol, jpol);
        }
    }
}

double FluidElement::measure(int count) const {
    // random disp
    int ipnt = 0;
    for (int ipol = 0; ipol <= nPol; ipol++) {
        for (int jpol = 0; jpol <= nPol; jpol++) {
            mPoints[ipnt++]->randomDispl((Real)1e-6);
        }
    }
    
    // measure stiffness
    MyBoostTimer timer;
    timer.start();
    for (int i = 0; i < count; i++) {
        computeStiff();
    }
    double elapsed_time = timer.elapsed();
    
    // reset point
    ipnt = 0;
    for (int ipol = 0; ipol <= nPol; ipol++) {
        for (int jpol = 0; jpol <= nPol; jpol++) {
            mPoints[ipnt++]->resetZero();
        }
    }
    return elapsed_time / count;
}

void FluidElement::test() const {
    // zero disp
    sResponse.setNr(mMaxNr);
    for (int alpha = 0; alpha <= mMaxNu; alpha++) {
        sResponse.mDispl[alpha].setZero();
    }
    // stiffness matrix
    int totalDim = (mMaxNu + 1) * nPntElem;
    RMatXX K = RMatXX::Zero(totalDim, totalDim);
    bool axial = this->axial();
    
    for (int alpha = 0; alpha <= mMaxNu; alpha++) {
        if (mMaxNr % 2 == 0 && alpha == mMaxNu) {continue;}
        for (int ipol = 0; ipol <= nPol; ipol++) {
            for (int jpol = 0; jpol <= nPol; jpol++) {
                // delta function
                if (axial && ipol == 0) {
                    if (alpha > 0) {continue;}
                }
                sResponse.mDispl[alpha](ipol, jpol) = one;
                if (alpha == 0) {sResponse.mDispl[alpha](ipol, jpol) = two;}
                
                // compute stiff 
                displToStiff();
                
                // positive-definite
                Real sr = sResponse.mStiff[alpha](ipol, jpol).real();
                if (sr <= zero) {
                    // add code here to debug
                    throw std::runtime_error("FluidElement::test || "
                        "Stiffness matrix is not positive definite.");  
                }
                    
                // store stiffness
                int row = alpha * nPntElem + ipol * nPntEdge + jpol;
                for (int alpha1 = 0; alpha1 <= mMaxNu; alpha1++) {
                    if (mMaxNr % 2 == 0 && alpha1 == mMaxNu) {continue;}
                    for (int ipol1 = 0; ipol1 <= nPol; ipol1++) {
                        for (int jpol1 = 0; jpol1 <= nPol; jpol1++) {
                            if (axial && ipol1 == 0) {
                                if (alpha1 > 0) {continue;}
                            }
                            int col = alpha1 * nPntElem + ipol1 * nPntEdge + jpol1;
                            K(row, col) = sResponse.mStiff[alpha1](ipol1, jpol1).real();
                        }
                    }
                }
                
                // restore zero 
                sResponse.mDispl[alpha](ipol, jpol) = zero;
            }
        }
    }
    
    // test self-adjointness 
    Real maxK = K.array().abs().maxCoeff();
    Real tole = maxK * tinyReal;
    for (int i = 0; i < totalDim; i++) {
        for (int j = i + 1; j < totalDim; j++) {
            Real diff = std::abs(K(i, j) - K(j, i));
            if (diff > tole) {
                // add code here to debug
                throw std::runtime_error("FluidElement::test || "
                    "Stiffness matrix is not self-adjoint."); 
            }
        }
    }
}

void FluidElement::computeGroundMotion(Real phi, const RMatPP &weights, RRow3 &u_spz) const {
    throw std::runtime_error("FluidElement::computeGroundMotion || "
        "Not implemented."); 
}

void FluidElement::computeStrain(Real phi, const RMatPP &weights, RRow6 &strain) const {
    throw std::runtime_error("FluidElement::computeStrain || "
        "Not implemented."); 
}

void FluidElement::forceTIso() {
    if (mCrdTransTIso == 0) {
        mCrdTransTIso = new CrdTransTIsoFluid(formThetaMat());
        mInTIso = true;
    } 
}

void FluidElement::feedDispOnSide(int side, CMatXX_RM &buffer, int row) const {
    throw std::runtime_error("FluidElement::getDispOnSide || "
        "Not implemented."); 
}

std::string FluidElement::verbose() const {
    if (mHasPRT) {
        return "FluidElement$" + mPRT->verbose() + "$" + mAcoustic->verbose();
    } else {
        return "FluidElement$" + mAcoustic->verbose();
    }
}

void FluidElement::displToStiff() const {
    mGradient->computeGrad(sResponse.mDispl, sResponse.mStrain, sResponse.mNu, sResponse.mNyquist);
    if (mInTIso) {
        mCrdTransTIso->transformSPZ_RTZ(sResponse.mStrain, sResponse.mNu);
    }
    if (mElem3D) {
        FieldFFT::transformF2P(sResponse.mStrain, sResponse.mNr);
    }
    if (mHasPRT) {
        mPRT->sphericalToUndulated(sResponse);
    }    
    mAcoustic->strainToStress(sResponse);
    if (mHasPRT) {
        mPRT->undulatedToSpherical(sResponse);
    }
    if (mElem3D) {
        FieldFFT::transformP2F(sResponse.mStress, sResponse.mNr);
    }
    if (mInTIso) {
        mCrdTransTIso->transformRTZ_SPZ(sResponse.mStress, sResponse.mNu);
    }
    mGradient->computeQuad(sResponse.mStiff, sResponse.mStress, sResponse.mNu, sResponse.mNyquist);
}

//-------------------------- static --------------------------//
FluidResponse FluidElement::sResponse;
void FluidElement::initWorkspace(int maxMaxNu) {
    sResponse.mDispl = vec_CMatPP(maxMaxNu + 1, CMatPP::Zero());
    sResponse.mStiff = vec_CMatPP(maxMaxNu + 1, CMatPP::Zero());
    sResponse.mStrain = vec_ar3_CMatPP(maxMaxNu + 1, zero_ar3_CMatPP);
    sResponse.mStress = vec_ar3_CMatPP(maxMaxNu + 1, zero_ar3_CMatPP);
}

