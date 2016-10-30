// FluidElement.cpp
// created by Kuangdai on 29-Mar-2016 
// base class of fluid elements

#include "FluidElement.h"
#include "Acoustic.h"
#include "Point.h"
#include "Gradient.h"
#include <boost/timer/timer.hpp>

FluidElement::FluidElement(Gradient *grad, const std::array<Point *, nPntElem> &points, 
    Acoustic *acous): 
Element(grad, points), mAcoustic(acous) {
    mAcoustic->checkCompatibility(mMaxNr);
}

FluidElement::~FluidElement() {
    delete mAcoustic;
}
    
void FluidElement::computeStiff() const {
    // get displ from points
    int ipnt = 0;
    for (int ipol = 0; ipol <= nPol; ipol++)
        for (int jpol = 0; jpol <= nPol; jpol++)
            mPoints[ipnt++]->scatterDisplToElement(sDispl, ipol, jpol, mMaxNu);
        
    // compute stiff
    displToStiff(sDispl, sStiff);
    
    // set stiff to points
    ipnt = 0;
    for (int ipol = 0; ipol <= nPol; ipol++)
        for (int jpol = 0; jpol <= nPol; jpol++)
            mPoints[ipnt++]->gatherStiffFromElement(sStiff, ipol, jpol);
}

double FluidElement::measure(int count, bool user) const {
    // random disp
    int ipnt = 0;
    for (int ipol = 0; ipol <= nPol; ipol++)
        for (int jpol = 0; jpol <= nPol; jpol++)
            mPoints[ipnt++]->randomDispl((Real)1e-6);
    
    // measure stiffness
    boost::timer::cpu_timer timer;
    for (int i = 0; i < count; i++) computeStiff();
    double elapsed_time = user ? timer.elapsed().user * 1.0 : timer.elapsed().wall * 1.0;
    
    // reset point
    ipnt = 0;
    for (int ipol = 0; ipol <= nPol; ipol++)
        for (int jpol = 0; jpol <= nPol; jpol++)
            mPoints[ipnt++]->resetZero();
    return elapsed_time / count;
}

void FluidElement::test() const {
    // zero disp
    for (int alpha = 0; alpha <= mMaxNu; alpha++)
        sDispl[alpha].setZero();
    // stiffness matrix
    int totalDim = (mMaxNu + 1) * nPntElem;
    RMatXX K = RMatXX::Zero(totalDim, totalDim);
    // axial 
    bool axial = getPoint(0)->axial();
    
    for (int alpha = 0; alpha <= mMaxNu; alpha++) {
        if (mMaxNr % 2 == 0 && alpha == mMaxNu) continue;
        for (int ipol = 0; ipol <= nPol; ipol++) {
            for (int jpol = 0; jpol <= nPol; jpol++) {
                // delta function
                if (axial && ipol == 0) {
                    if (alpha > 0) continue;
                }
                sDispl[alpha](ipol, jpol) = one;
                if (alpha == 0) sDispl[alpha](ipol, jpol) = two;
                
                // compute stiff 
                displToStiff(sDispl, sStiff);
                
                // positive-definite
                Real sr = sStiff[alpha](ipol, jpol).real();
                if (sr <= zero) {
                    // add code here to debug
                    throw std::runtime_error("FluidElement::test || "
                        "Stiffness matrix is not positive definite.");  
                }
                    
                // store stiffness
                int row = alpha * nPntElem + ipol * nPntEdge + jpol;
                for (int alpha1 = 0; alpha1 <= mMaxNu; alpha1++) {
                    if (mMaxNr % 2 == 0 && alpha1 == mMaxNu) continue;
                    for (int ipol1 = 0; ipol1 <= nPol; ipol1++) {
                        for (int jpol1 = 0; jpol1 <= nPol; jpol1++) {
                            if (axial && ipol1 == 0) {
                                if (alpha1 > 0) continue;
                            }
                            int col = alpha1 * nPntElem + ipol1 * nPntEdge + jpol1;
                            K(row, col) = sStiff[alpha1](ipol1, jpol1).real();
                        }
                    }
                }
                
                // restore zero 
                sDispl[alpha](ipol, jpol) = zero;
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
    // get chi
    int ipnt = 0;
    for (int ipol = 0; ipol <= nPol; ipol++)
        for (int jpol = 0; jpol <= nPol; jpol++)
            mPoints[ipnt++]->scatterDisplToElement(sDispl, ipol, jpol, mMaxNu);
    // u = nabla(chi) / rho       
    mGradient->gradScalar(sDispl, sStrain, mMaxNu, mMaxNr % 2 == 0);
    mAcoustic->strainToStress(sStrain, sStress, mMaxNu);
    // compute ground motion pointwise
    u_spz.setZero();
    for (int ipol = 0; ipol <= nPol; ipol++) {
        for (int jpol = 0; jpol <= nPol; jpol++) {
            if (std::abs(weights(ipol, jpol)) < tinyDouble) continue;
            Real up0 = sStress[0][0](ipol, jpol).real();
            Real up1 = sStress[0][1](ipol, jpol).real();
            Real up2 = sStress[0][2](ipol, jpol).real();
            for (int alpha = 1; alpha <= mMaxNu - (int)(mMaxNr % 2 == 0); alpha++) {
                Complex expval = two * exp((Real)alpha * phi * ii);
                up0 += (expval * sStress[alpha][0](ipol, jpol)).real();
                up1 += (expval * sStress[alpha][1](ipol, jpol)).real();
                up2 += (expval * sStress[alpha][2](ipol, jpol)).real();
            }
            u_spz(0) += weights(ipol, jpol) * up0;
            u_spz(1) += weights(ipol, jpol) * up1;
            u_spz(2) += weights(ipol, jpol) * up2;
        }
    }
}

std::string FluidElement::verbose() const {
    return "FluidElement$" + mAcoustic->verbose();
}

void FluidElement::displToStiff(const vec_CMatPP &displ, vec_CMatPP &stiff) const {
    mGradient->gradScalar(displ, sStrain, mMaxNu, mMaxNr % 2 == 0);
    mAcoustic->strainToStress(sStrain, sStress, mMaxNu);
    mGradient->quadScalar(sStress, stiff, mMaxNu, mMaxNr % 2 == 0);
}

//-------------------------- static --------------------------//
vec_CMatPP FluidElement::sDispl;
vec_CMatPP FluidElement::sStiff;
vec_ar3_CMatPP FluidElement::sStrain;
vec_ar3_CMatPP FluidElement::sStress;
void FluidElement::initWorkspace(int maxMaxNu) {
    sDispl = vec_CMatPP(maxMaxNu + 1, CMatPP::Zero());
    sStiff = vec_CMatPP(maxMaxNu + 1, CMatPP::Zero());
    sStrain = vec_ar3_CMatPP(maxMaxNu + 1, zero_ar3_CMatPP);
    sStress = vec_ar3_CMatPP(maxMaxNu + 1, zero_ar3_CMatPP);
}

