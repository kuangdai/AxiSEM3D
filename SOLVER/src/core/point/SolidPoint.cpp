// SolidPoint.h
// created by Kuangdai on 3-Apr-2016 
// solid gll points 

#include "SolidPoint.h"
#include "Mass.h"
#include "MultilevelTimer.h"

SolidPoint::SolidPoint(int nr, bool axial, const RDCol2 &crds, Mass *mass):
Point(nr, axial, crds), mMass(mass) {
    mDispl = CMatX3::Zero(mNu + 1, 3);
    mVeloc = CMatX3::Zero(mNu + 1, 3);
    mAccel = CMatX3::Zero(mNu + 1, 3);
    mStiff = CMatX3::Zero(mNu + 1, 3);
    mMass->checkCompatibility(nr);
    mNuWisdom.fill(mNu);
}

SolidPoint::~SolidPoint() {
    delete mMass;
}

void SolidPoint::updateNewmark(double dt) {
    // mask stiff 
    maskField(mStiff);
    // compute accel inplace
    mMass->computeAccel(mStiff);
    // mask accel (masking must be called twice if mass is 3D)
    maskField(mStiff);
    // update dt
    double half_dt = half * dt;
    double half_dt_dt = half_dt * dt;
    mVeloc += (Real)half_dt * (mAccel + mStiff);
    mAccel = mStiff;
    mDispl += (Real)dt * mVeloc + (Real)half_dt_dt * mAccel;  
    // zero stiffness for next time step
    mStiff.setZero();
}

void SolidPoint::resetZero() {
    mStiff.setZero();
    mDispl.setZero();
    mVeloc.setZero();
    mAccel.setZero();
}

void SolidPoint::randomDispl(Real factor, int seed, int max_order) {
    if (seed >= 0) {
        std::srand(seed);
    }
    if (max_order < 0 || max_order > mNu) {
        mDispl.setRandom(); 
    } else {
        mDispl.topRows(max_order + 1).setRandom();
    }
    mDispl *= factor;
    maskField(mDispl);
}

void SolidPoint::randomStiff(Real factor, int seed, int max_order) {
    if (seed >= 0) {
        std::srand(seed);
    }
    if (max_order < 0 || max_order > mNu) {
        mStiff.setRandom(); 
    } else {
        mStiff.topRows(max_order + 1).setRandom();
    }
    mStiff *= factor;
    maskField(mStiff);
}

std::string SolidPoint::verbose() const {
    return "SolidPoint$" + mMass->verbose();
}

double SolidPoint::measure(int count) {
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

void SolidPoint::test() {
    // mass matrix
    int totalDim = (mNu + 1) * 3;
    RMatXX M = RMatXX::Zero(totalDim, totalDim);
    
    resetZero();
    for (int alpha = 0; alpha <= mNu; alpha++) {
        if (mNr % 2 == 0 && alpha == mNu) {continue;}
        for (int idim = 0; idim <= 2; idim++) {
            // delta function
            if (axial()) {
                if (alpha == 0 && idim != 2) {continue;}
                if (alpha == 1 && idim == 2) {continue;}
                if (alpha >= 2) {continue;}
            }
            mStiff(alpha, idim) = one;
            if (alpha == 0) {mStiff(alpha, idim) = two;}
                    
            // compute stiff 
            updateNewmark(1.);
                    
            // positive-definite
            Real sr = mAccel(alpha, idim).real();
            if (sr <= zero) {
                // add code here to debug
                throw std::runtime_error("SolidPoint::test || "
                    "Mass matrix is not positive definite.");  
            }
                        
            // store mass
            int row = alpha * 3 + idim;
            for (int alpha1 = 0; alpha1 <= mNu; alpha1++) {
                if (mNr % 2 == 0 && alpha1 == mNu) {continue;}
                for (int idim1 = 0; idim1 <= 2; idim1++) {
                    if (axial()) {
                        if (alpha1 == 0 && idim1 != 2) {continue;}
                        if (alpha1 == 1 && idim1 == 2) {continue;}
                        if (alpha1 >= 2) {continue;}
                    }
                    int col = alpha1 * 3 + idim1;
                    M(row, col) = mAccel(alpha1, idim1).real();
                }
            }
                    
            // restore zero 
            mStiff(alpha, idim) = czero;
        }
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
                throw std::runtime_error("SolidPoint::test || "
                    "Mass matrix is not self-adjoint."); 
            }
        }
    }
}

void SolidPoint::feedBuffer(CColX &buffer, int &row) {
    int size = mStiff.size();
    buffer.block(row, 0, size, 1) = Eigen::Map<CColX>(mStiff.data(), size);
    row += size;
}

void SolidPoint::extractBuffer(CColX &buffer, int &row) {
    int size = mStiff.size();
    mStiff += Eigen::Map<CMatX3>(buffer.block(row, 0, size, 1).data(), mStiff.rows(), 3);
    row += size;
}

void SolidPoint::scatterDisplToElement(vec_ar3_CMatPP &displ, int ipol, int jpol, int maxNu) const {
    // lower orders
    int nyquist = (int)(mNr % 2 == 0); 
    for (int alpha = 0; alpha <= mNu - nyquist; alpha++) {
        displ[alpha][0](ipol, jpol) = mDispl(alpha, 0);
        displ[alpha][1](ipol, jpol) = mDispl(alpha, 1);
        displ[alpha][2](ipol, jpol) = mDispl(alpha, 2);
    }
    // mask Nyquist 
    if (nyquist) {
        displ[mNu][0](ipol, jpol) = czero;
        displ[mNu][1](ipol, jpol) = czero;
        displ[mNu][2](ipol, jpol) = czero;
    }
    // mask higher orders
    for (int alpha = mNu + 1; alpha <= maxNu; alpha++) {
        displ[alpha][0](ipol, jpol) = czero;
        displ[alpha][1](ipol, jpol) = czero;
        displ[alpha][2](ipol, jpol) = czero;
    }
}

void SolidPoint::gatherStiffFromElement(const vec_ar3_CMatPP &stiff, int ipol, int jpol) {
    // lower orders
    int nyquist = (int)(mNr % 2 == 0); 
    for (int alpha = 0; alpha <= mNu - nyquist; alpha++) {
        mStiff(alpha, 0) -= stiff[alpha][0](ipol, jpol);
        mStiff(alpha, 1) -= stiff[alpha][1](ipol, jpol);
        mStiff(alpha, 2) -= stiff[alpha][2](ipol, jpol);
    }
    // mask Nyquist 
    if (nyquist) {
        mStiff.row(mNu).setZero();
    }
}

void SolidPoint::addToStiff(const CMatX3 &source) {
    // make sure the length of "source" does not exceed mNu + 1 
    mStiff.topRows(source.rows()) += source;
}

void SolidPoint::maskField(CMatX3 &field) {
    field.row(0).imag().setZero();
    // axial boundary condition
    if (mAxial) {
        // alpha = 0
        field(0, 0) = czero;
        field(0, 1) = czero;
        if (mNu >= 1) {
            // alpha = 1
            Complex s0 = field(1, 0);
            Complex s1 = field(1, 1);
            field(1, 0) = half * (s0 - ii * s1);
            field(1, 1) = half * (s1 + ii * s0);
            field(1, 2) = czero;
            // alpha > 1
            field.bottomRows(mNu - 1).setZero();
        }
    }
    // mask Nyquist 
    if (mNr % 2 == 0) {
        field.row(mNu).setZero();
    }
}

void SolidPoint::learnWisdom(Real cutoff) {
    for (int idim = 0; idim < 3; idim++) {
        // L2 norm
        Real L2norm = mDispl.col(idim).squaredNorm();
        // Hilbert norm
        Real h2norm = L2norm - .5 * mDispl.col(idim).row(0).squaredNorm();
        if (h2norm <= mMaxDisplWisdom(idim)) {
            continue;
        }
        mMaxDisplWisdom(idim) = h2norm;
        
        // try smaller orders
        Real tol = h2norm * cutoff * cutoff;
        bool found = false;
        for (int newNu = 0; newNu < mNu; newNu++) {
            Real diff = L2norm - mDispl.col(idim).topRows(newNu + 1).squaredNorm();
            if (diff <= tol) {
                mNuWisdom(idim) = newNu;
                found = true;
                break;
            }
        }
        if (!found) {
            mNuWisdom(idim) = mNu;
        }
    }
}

int SolidPoint::getNuWisdom() const {
    // int maxloc = 0;
    // mMaxDisplWisdom.maxCoeff(&maxloc);
    return mNuWisdom.maxCoeff();
}

