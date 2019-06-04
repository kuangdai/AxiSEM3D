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
    // setup static
    sResponse.setNr(mMaxNr);
    
    // get displ from points
    int ipnt = 0;
    for (int ipol = 0; ipol <= nPol; ipol++) {
        for (int jpol = 0; jpol <= nPol; jpol++) {
            mPoints[ipnt++]->scatterDisplToElement(sResponse.mDispl, ipol, jpol, mMaxNu);
        }
    }
    
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
    
    // compute ground motion pointwise
    u_spz.setZero();
    for (int ipol = 0; ipol <= nPol; ipol++) {
        for (int jpol = 0; jpol <= nPol; jpol++) {
            // if (std::abs(weights(ipol, jpol)) < tinyDouble) continue;
            Real up0 = sResponse.mStress[0][0](ipol, jpol).real();
            Real up1 = sResponse.mStress[0][1](ipol, jpol).real();
            Real up2 = sResponse.mStress[0][2](ipol, jpol).real();
            for (int alpha = 1; alpha <= mMaxNu - (int)(mMaxNr % 2 == 0); alpha++) {
                Complex expval = two * exp((Real)alpha * phi * ii);
                up0 += (expval * sResponse.mStress[alpha][0](ipol, jpol)).real();
                up1 += (expval * sResponse.mStress[alpha][1](ipol, jpol)).real();
                up2 += (expval * sResponse.mStress[alpha][2](ipol, jpol)).real();
            }
            u_spz(0) += weights(ipol, jpol) * up0;
            u_spz(1) += weights(ipol, jpol) * up1;
            u_spz(2) += weights(ipol, jpol) * up2;
        }
    }
}

#include "SolidElement.h"
void FluidElement::computeStrain(Real phi, const RMatPP &weights, RRow6 &strain) const {
    // setup static
    sResponse.setNr(mMaxNr);
    
    // get displ from points
    int ipnt = 0;
    for (int ipol = 0; ipol <= nPol; ipol++) {
        for (int jpol = 0; jpol <= nPol; jpol++) {
            mPoints[ipnt++]->scatterDisplToElement(sResponse.mDispl, ipol, jpol, mMaxNu);
        }
    }
    
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
    
    //////////////////////
    if (mHasPRT) {
        throw std::runtime_error("FluidElement::computeStrain || "
            "Not implemented."); 
    } else {
        mGradient->computeGrad6(sResponse.mStress, SolidElement::sResponse.mStrain6, 
            sResponse.mNu, sResponse.mNyquist);
    }
    
    //////////////
    strain.setZero();
    for (int ipol = 0; ipol <= nPol; ipol++) {
        for (int jpol = 0; jpol <= nPol; jpol++) {
            if (std::abs(weights(ipol, jpol)) < tinyDouble) continue;
            Real s0 = SolidElement::sResponse.mStrain6[0][0](ipol, jpol).real();
            Real s1 = SolidElement::sResponse.mStrain6[0][1](ipol, jpol).real();
            Real s2 = SolidElement::sResponse.mStrain6[0][2](ipol, jpol).real();
            Real s3 = SolidElement::sResponse.mStrain6[0][3](ipol, jpol).real();
            Real s4 = SolidElement::sResponse.mStrain6[0][4](ipol, jpol).real();
            Real s5 = SolidElement::sResponse.mStrain6[0][5](ipol, jpol).real();
            for (int alpha = 1; alpha <= mMaxNu - (int)(mMaxNr % 2 == 0); alpha++) {
                Complex expval = two * exp((Real)alpha * phi * ii);
                s0 += (expval * SolidElement::sResponse.mStrain6[alpha][0](ipol, jpol)).real();
                s1 += (expval * SolidElement::sResponse.mStrain6[alpha][1](ipol, jpol)).real();
                s2 += (expval * SolidElement::sResponse.mStrain6[alpha][2](ipol, jpol)).real();
                s3 += (expval * SolidElement::sResponse.mStrain6[alpha][3](ipol, jpol)).real();
                s4 += (expval * SolidElement::sResponse.mStrain6[alpha][4](ipol, jpol)).real();
                s5 += (expval * SolidElement::sResponse.mStrain6[alpha][5](ipol, jpol)).real();
            }
            strain(0) += weights(ipol, jpol) * s0;
            strain(1) += weights(ipol, jpol) * s1;
            strain(2) += weights(ipol, jpol) * s2;
            strain(3) += weights(ipol, jpol) * s3;
            strain(4) += weights(ipol, jpol) * s4;
            strain(5) += weights(ipol, jpol) * s5;
        }
    }
}

void FluidElement::computeCurl(Real phi, const RMatPP &weights, RRow3 &curl) const {
    // no curl in fluid
    curl.setZero();
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

