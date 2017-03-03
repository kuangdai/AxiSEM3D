// FluidPoint.cpp
// created by Kuangdai on 4-Apr-2016 
// fluid gll points 

#include "FluidPoint.h"
#include "Mass.h"
#include "XTimer.h"

FluidPoint::FluidPoint(int nr, bool axial, const RDCol2 &crds, Mass *mass):
Point(nr, axial, crds), mMass(mass) {
    mDispl = CColX::Zero(mNu + 1, 1);
    mVeloc = CColX::Zero(mNu + 1, 1);
    mAccel = CColX::Zero(mNu + 1, 1);
    mStiff = CColX::Zero(mNu + 1, 1);
    mMass->checkCompatibility(nr);
}

FluidPoint::~FluidPoint() {
    delete mMass;
}

void FluidPoint::updateNewmark(Real dt) {
    // mask stiff 
    maskField(mStiff);
    // compute accel inplace
    mMass->computeAccel(mStiff);
    // mask accel (masking must be called twice if mass is 3D)
    maskField(mStiff);
    // update dt
    Real half_dt = half * dt;
    Real half_dt_dt = half_dt * dt;
    mVeloc += half_dt * (mAccel + mStiff);
    mAccel = mStiff;
    mDispl += dt * mVeloc + half_dt_dt * mAccel;  
    // zero stiffness for next time step
    mStiff.setZero();
}

void FluidPoint::resetZero() {
    mStiff.setZero();
    mDispl.setZero();
    mVeloc.setZero();
    mAccel.setZero();
}

void FluidPoint::randomDispl(Real factor, int seed) {
    std::srand(seed); 
    mDispl.setRandom(); 
    mDispl *= factor;
    maskField(mDispl);
}

void FluidPoint::randomStiff(Real factor, int seed) {
    std::srand(seed); 
    mStiff.setRandom(); 
    mStiff *= factor;
    maskField(mStiff);
}

std::string FluidPoint::verbose() const {
    return "FluidPoint$" + mMass->verbose();
}

double FluidPoint::measure(int count) {
    Real dt = .1;
    randomStiff((Real)1e6); 
    MyBoostTimer timer;
    timer.start();
    for (int i = 0; i < count; i++) {
        maskField(mStiff);
        mMass->computeAccel(mStiff);
        maskField(mStiff);
        Real half_dt = half * dt;
        Real half_dt_dt = half_dt * dt;
        mVeloc += half_dt * (mAccel + mStiff);
        mAccel = mStiff;
        mDispl += dt * mVeloc + half_dt_dt * mAccel; 
        mVeloc.setZero();
    }
    double elapsed_time = timer.elapsed();
    resetZero();
    return elapsed_time / count;
}

void FluidPoint::test() {
    // mass matrix
    int totalDim = mNu + 1;
    RMatXX M = RMatXX::Zero(totalDim, totalDim);
    
    resetZero();
    for (int alpha = 0; alpha <= mNu; alpha++) {
        if (mNr % 2 == 0 && alpha == mNu) continue;
        // delta function
        if (axial() && alpha > 0) continue;
        mStiff(alpha) = one;
        if (alpha == 0) mStiff(alpha) = two;
                
        // compute stiff 
        updateNewmark(1.);
                
        // positive-definite
        Real sr = mAccel(alpha).real();
        if (sr <= zero) {
            // add code here to debug
            throw std::runtime_error("FluidPoint::test || "
                "Mass matrix is not positive definite.");  
        }
                    
        // store mass
        int row = alpha;
        for (int alpha1 = 0; alpha1 <= mNu; alpha1++) {
            if (mNr % 2 == 0 && alpha1 == mNu) continue;
            if (axial() && alpha > 0) continue;
            int col = alpha1;
            M(row, col) = mAccel(alpha1).real();
        }
                
        // restore zero 
        mStiff(alpha) = czero;
    }
    resetZero();

    // test self-adjointness 
    Real maxM = M.array().abs().maxCoeff();
    Real tole = maxM * tinyReal;
    for (int i = 0; i < totalDim; i++) {
        for (int j = i + 1; j < totalDim; j++) {
            Real diff = std::abs(M(i, j) - M(j, i));
            if (diff > tole) {
                // add code here to debug
                throw std::runtime_error("FluidPoint::test || "
                    "Mass matrix is not self-adjoint."); 
            }
        }
    }
}

void FluidPoint::feedBuffer(CColX &buffer, int &row) {
    int rows = mStiff.rows();
    buffer.block(row, 0, rows, 1) = mStiff;
    row += rows;
}

void FluidPoint::extractBuffer(CColX &buffer, int &row) {
    int rows = mStiff.rows();
    mStiff += buffer.block(row, 0, rows, 1);
    row += rows;
}

void FluidPoint::scatterDisplToElement(vec_CMatPP &displ, int ipol, int jpol, int maxNu) const {
    // lower orders
    int nyquist = (int)(mNr % 2 == 0);
    for (int alpha = 0; alpha <= mNu - nyquist; alpha++)
        displ[alpha](ipol, jpol) = mDispl(alpha);
    // mask Nyquist
    if (nyquist) displ[mNu](ipol, jpol) = czero;
    // mask higher orders
    for (int alpha = mNu + 1; alpha <= maxNu; alpha++)
        displ[alpha](ipol, jpol) = czero;
}

void FluidPoint::gatherStiffFromElement(const vec_CMatPP &stiff, int ipol, int jpol) {
    // lower orders
    int nyquist = (int)(mNr % 2 == 0);
    for (int alpha = 0; alpha <= mNu - nyquist; alpha++)
        mStiff(alpha) -= stiff[alpha](ipol, jpol);
    // mask Nyquist
    if (nyquist) mStiff(mNu) = czero;
}

void FluidPoint::maskField(CColX &field) {
    field.row(0).imag().setZero();
    // axial boundary condition
    if (mAxial) field.bottomRows(mNu).setZero();
    // mask Nyquist
    if (mNr % 2 == 0) field(mNu) = czero;
}

#include "SolverFFTW_1.h"
void FluidPoint::learnWisdom(double cutoff) {
    // compute real displ
    SolverFFTW_1::getC2R_CMat(mNr) = mDispl;
    SolverFFTW_1::computeC2R(mNr);
    RColX &displ = SolverFFTW_1::getC2R_RMat(mNr);
    
    // check max displ
    Real maxDispl = displ.array().abs().maxCoeff();
    if (maxDispl <= mMaxDisplWisdom) return;
    mMaxDisplWisdom = maxDispl;
    
    // keep original
    RColX &displOrig = SolverFFTW_1::getR2C_RMat(mNr);
    displOrig = displ;
    Real normTol = cutoff * displOrig.norm();
    
    // try smaller orders
    for (int newNu = 0; newNu < mNu; newNu++) {
        SolverFFTW_1::getC2R_CMat(mNr) = mDispl; // C2R destroys input
        SolverFFTW_1::getC2R_CMat(mNr).bottomRows(mNu - newNu).setZero();
        SolverFFTW_1::computeC2R(mNr);
        RColX &dispNew = SolverFFTW_1::getC2R_RMat(mNr);
        if ((dispNew - displOrig).norm() <= normTol) {
            mNuWisdom = newNu;
            return;
        }
    }
    mNuWisdom = mNu;
}


