// SolidElement.cpp
// created by Kuangdai on 29-Mar-2016 
// base class of solid element

#include "SolidElement.h"
#include "Point.h"
#include "Gradient.h"
#include "PRT.h"
#include "Elastic.h"
#include "CrdTransTIsoSolid.h"
#include "FieldFFT.h"

#include "MultilevelTimer.h"

SolidElement::SolidElement(Gradient *grad, PRT *prt,
    const std::array<Point *, nPntElem> &points, 
    Elastic *elas): 
Element(grad, prt, points), mElastic(elas), mCrdTransTIso(0) {
    mElastic->checkCompatibility(mMaxNr);
    // TISO
    mInTIso = mHasPRT;
    if (mInTIso) {
        mCrdTransTIso = new CrdTransTIsoSolid(formThetaMat());
    } 
    // 3D
    bool elas1D = mElastic->is1D();
    bool prt1D = false;
    if (mHasPRT) {
        prt1D = mPRT->is1D();
    } else {
        prt1D = true;
    }
    mElem3D = !(prt1D && elas1D);
}

SolidElement::~SolidElement() {
    delete mElastic;
    if (mInTIso) {
        delete mCrdTransTIso;
    }
}

void SolidElement::computeStiff() const {
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

double SolidElement::measure(int count) const {
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

void SolidElement::test() const {
    // zero disp
    for (int alpha = 0; alpha <= mMaxNu; alpha++) {
        sResponse.mDispl[alpha][0].setZero();
        sResponse.mDispl[alpha][1].setZero();
        sResponse.mDispl[alpha][2].setZero();    
    }
    // stiffness matrix
    int totalDim = (mMaxNu + 1) * 3 * nPntElem;
    RMatXX K = RMatXX::Zero(totalDim, totalDim);
    bool axial = this->axial();
    
    for (int alpha = 0; alpha <= mMaxNu; alpha++) {
        if (mMaxNr % 2 == 0 && alpha == mMaxNu) {continue;}
        for (int idim = 0; idim <= 2; idim++) {
            for (int ipol = 0; ipol <= nPol; ipol++) {
                for (int jpol = 0; jpol <= nPol; jpol++) {
                    // delta function
                    if (axial && ipol == 0) {
                        if (alpha == 0 && idim != 2) {continue;}
                        if (alpha == 1 && idim == 2) {continue;}
                        if (alpha >= 2) {continue;}
                    }
                    sResponse.mDispl[alpha][idim](ipol, jpol) = one;
                    if (alpha == 0) {sResponse.mDispl[alpha][idim](ipol, jpol) = two;}
                    
                    // compute stiff 
                    displToStiff();
                    
                    // positive-definite
                    Real sr = sResponse.mStiff[alpha][idim](ipol, jpol).real();
                    if (sr <= zero) {
                        // add code here to debug
                        throw std::runtime_error("SolidElement::test || "
                            "Stiffness matrix is not positive definite.");  
                    } 
                        
                    // store stiffness 
                    int row = alpha * nPntElem * 3 + idim * nPntElem + ipol * nPntEdge + jpol;
                    for (int alpha1 = 0; alpha1 <= mMaxNu; alpha1++) {
                        if (mMaxNr % 2 == 0 && alpha1 == mMaxNu) {continue;}
                        for (int idim1 = 0; idim1 <= 2; idim1++) {
                            for (int ipol1 = 0; ipol1 <= nPol; ipol1++) {
                                for (int jpol1 = 0; jpol1 <= nPol; jpol1++) {
                                    if (axial && ipol1 == 0) {
                                        if (alpha1 == 0 && idim1 != 2) {continue;}
                                        if (alpha1 == 1 && idim1 == 2) {continue;}
                                        if (alpha1 >= 2) {continue;}
                                    }
                                    int col = alpha1 * nPntElem * 3 + idim1 * nPntElem + ipol1 * nPntEdge + jpol1;
                                    K(row, col) = sResponse.mStiff[alpha1][idim1](ipol1, jpol1).real();
                                }
                            }
                        }
                    }
                    
                    // reset disp to zero 
                    sResponse.mDispl[alpha][idim](ipol, jpol) = zero;
                    // reset memory variables to zero
                    mElastic->resetZero();
                }
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
                throw std::runtime_error("SolidElement::test || "
                    "Stiffness matrix is not self-adjoint.");  
                // int alpha = i / (nPntElem * 3);
                // int idim = (i - alpha * nPntElem * 3) / nPntElem;
                // int ipol = (i - alpha * nPntElem * 3 - idim * nPntElem) / nPntEdge;
                // int jpol = i - alpha * nPntElem * 3 - idim * nPntElem - ipol * nPntEdge;
                // 
                // int alpha1 = j / (nPntElem * 3);
                // int idim1 = (j - alpha1 * nPntElem * 3) / nPntElem;
                // int ipol1 = (j - alpha1 * nPntElem * 3 - idim1 * nPntElem) / nPntEdge;
                // int jpol1 = j - alpha1 * nPntElem * 3 - idim1 * nPntElem - ipol1 * nPntEdge;
                // 
                // std::cout << getDomainTag() << std::endl;
                // std::cout << alpha << " " << idim << " " << ipol << " " << jpol << std::endl;
                // std::cout << alpha1 << " " << idim1 << " " << ipol1 << " " << jpol1 << std::endl << std::endl;
            }
        }
    }
}

void SolidElement::computeGroundMotion(Real phi, const RMatPP &weights, RRow3 &u_spz) const {
    // get displ from points
    int ipnt = 0;
    for (int ipol = 0; ipol <= nPol; ipol++) {
        for (int jpol = 0; jpol <= nPol; jpol++) {
            mPoints[ipnt++]->scatterDisplToElement(sResponse.mDispl, ipol, jpol, mMaxNu);
        }
    }
    // compute ground motion pointwise
    u_spz.setZero();
    for (int ipol = 0; ipol <= nPol; ipol++) {
        for (int jpol = 0; jpol <= nPol; jpol++) {
            if (std::abs(weights(ipol, jpol)) < tinyDouble) continue;
            Real up0 = sResponse.mDispl[0][0](ipol, jpol).real();
            Real up1 = sResponse.mDispl[0][1](ipol, jpol).real();
            Real up2 = sResponse.mDispl[0][2](ipol, jpol).real();
            for (int alpha = 1; alpha <= mMaxNu - (int)(mMaxNr % 2 == 0); alpha++) {
                Complex expval = two * exp((Real)alpha * phi * ii);
                up0 += (expval * sResponse.mDispl[alpha][0](ipol, jpol)).real();
                up1 += (expval * sResponse.mDispl[alpha][1](ipol, jpol)).real();
                up2 += (expval * sResponse.mDispl[alpha][2](ipol, jpol)).real();
            }
            u_spz(0) += weights(ipol, jpol) * up0;
            u_spz(1) += weights(ipol, jpol) * up1;
            u_spz(2) += weights(ipol, jpol) * up2;
        }
    }
}

std::string SolidElement::verbose() const {
    if (mHasPRT) {
        return "SolidElement$" + mPRT->verbose() + "$" + mElastic->verbose();
    } else {
        return "SolidElement$" + mElastic->verbose();
    }
}

void SolidElement::displToStiff() const {
    if (mHasPRT) {
        mGradient->computeGrad9(sResponse.mDispl, sResponse.mStrain9, sResponse.mNu, sResponse.mNyquist);
        if (mInTIso) {
            mCrdTransTIso->transformSPZ_RTZ(sResponse.mStrain9, sResponse.mNu);
        }    
        if (mElem3D) {
            FieldFFT::transformF2P(sResponse.mStrain9, sResponse.mNr);
        }
        mPRT->sphericalToUndulated(sResponse);
    } else {
        mGradient->computeGrad6(sResponse.mDispl, sResponse.mStrain6, sResponse.mNu, sResponse.mNyquist);
        if (mInTIso) {
            mCrdTransTIso->transformSPZ_RTZ(sResponse.mStrain6, sResponse.mNu);
        }    
        if (mElem3D) {
            FieldFFT::transformF2P(sResponse.mStrain6, sResponse.mNr);
        }
    }
    mElastic->strainToStress(sResponse);
    if (mHasPRT) {
        mPRT->undulatedToSpherical(sResponse);
        if (mElem3D) {
            FieldFFT::transformP2F(sResponse.mStress9, sResponse.mNr);
        }
        if (mInTIso) {
            mCrdTransTIso->transformRTZ_SPZ(sResponse.mStress9, sResponse.mNu);
        }
        mGradient->computeQuad9(sResponse.mStiff, sResponse.mStress9, sResponse.mNu, sResponse.mNyquist);
    } else {
        if (mElem3D) {
            FieldFFT::transformP2F(sResponse.mStress6, sResponse.mNr);
        }
        if (mInTIso) {
            mCrdTransTIso->transformRTZ_SPZ(sResponse.mStress6, sResponse.mNu);
        }
        mGradient->computeQuad6(sResponse.mStiff, sResponse.mStress6, sResponse.mNu, sResponse.mNyquist);
    }
    
}

//-------------------------- static --------------------------//
SolidResponse SolidElement::sResponse;
void SolidElement::initWorkspace(int maxMaxNu) {
    sResponse.mDispl = vec_ar3_CMatPP(maxMaxNu + 1, zero_ar3_CMatPP);
    sResponse.mStiff = vec_ar3_CMatPP(maxMaxNu + 1, zero_ar3_CMatPP);
    
    sResponse.mStrain6 = vec_ar6_CMatPP(maxMaxNu + 1, zero_ar6_CMatPP);
    sResponse.mStress6 = vec_ar6_CMatPP(maxMaxNu + 1, zero_ar6_CMatPP);
    
    sResponse.mStrain9 = vec_ar9_CMatPP(maxMaxNu + 1, zero_ar9_CMatPP);
    sResponse.mStress9 = vec_ar9_CMatPP(maxMaxNu + 1, zero_ar9_CMatPP);
}

