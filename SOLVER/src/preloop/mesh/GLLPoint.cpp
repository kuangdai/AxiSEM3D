// GLLPoint.cpp
// created by Kuangdai on 6-May-2016 
// general gll point 

#include "GLLPoint.h"
#include "Mass1D.h"
#include "Mass3D.h"
#include "SolidPoint.h"
#include "FluidPoint.h"
#include "SolidFluidPoint.h"
#include "SFCoupling1D.h"
#include "SFCoupling3D.h"
#include "Domain.h"
#include "XMath.h"
#include "Geodesy.h"

#include "MassOcean1D.h"
#include "MassOcean3D.h"
#include "OceanLoad3D.h"

GLLPoint::GLLPoint(): mNr(0) {
    // nothing
}

void GLLPoint::setup(int nr, bool axial, bool surface, const RDCol2 &crds, double distTol) {
    if (mNr > 0) {
        // already set
        if (mNr != nr || mIsAxial != axial || mOnSurface != surface || (crds - mCoords).norm() > distTol) {
            throw std::runtime_error("GLLPoint::setup || Conflict in GLL point setup.");
        }
        mReferenceCount++;    
    } else {
        // not set yet
        mNr = nr;
        mIsAxial = axial;
        mOnSurface = surface;
        mCoords = crds;
        mMassSolid = RDColX::Zero(mNr);
        mMassFluid = RDColX::Zero(mNr);
        mSFNormal = RDMatX3::Zero(mNr, 3);
        mSFNormal_assmble = RDMatX3::Zero(mNr, 3);
        mOceanDepth = RDColX::Zero(mNr);
        mSurfNormal = RDMatX3::Zero(mNr, 3);
        mReferenceCount = 1;
    }
}

int GLLPoint::release(Domain &domain) const {
    // solid fluid
    bool isSolid = mMassSolid.norm() > tinyDouble;
    bool isFluid = mMassFluid.norm() > tinyDouble; 
    
    // solid point ptr
    SolidPoint *solid;
    if (isSolid) {
        Mass *mass;
        if (mOceanDepth.array().abs().maxCoeff() > tinyDouble) {
            if (XMath::equalRows(mMassSolid) && XMath::equalRows(mOceanDepth) && XMath::equalRows(mSurfNormal)) {
                double theta = Geodesy::theta(mCoords);
                double mass_solid = mMassSolid(0);
                double mass_ocean = OceanLoad3D::mWaterDensity * mOceanDepth(0) * mSurfNormal.row(0).norm();
                mass = new MassOcean1D(mass_solid, mass_ocean, theta);
            } else {
                RDColX mass_ocean(mNr);
                RDMatX3 unit_normal(mNr, 3);
                for (int i = 0; i < mNr; i++) {
                    double area = mSurfNormal.row(i).norm();
                    mass_ocean(i) = OceanLoad3D::mWaterDensity * mOceanDepth(i) * area;
                    unit_normal.row(i) = mSurfNormal.row(i) / area;
                }
                mass = new MassOcean3D(mMassSolid, mass_ocean, unit_normal);
            }
        } else {
            if (XMath::equalRows(mMassSolid)) {
                mass = new Mass1D((Real)(1. / mMassSolid(0)));
            } else {
                const RDColX &invMass = mMassSolid.array().pow(-1.).matrix(); 
                mass = new Mass3D(invMass.cast<Real>());
            }
        }
        solid = new SolidPoint(mNr, mIsAxial, mCoords, mass);
    }
    
    // fluid point ptr
    FluidPoint *fluid;
    if (isFluid) {
        Mass *mass;
        if (XMath::equalRows(mMassFluid)) {
            mass = new Mass1D((Real)(1. / mMassFluid(0)));
        } else {
            const RDColX &invMass = mMassFluid.array().pow(-1.).matrix(); 
            mass = new Mass3D(invMass.cast<Real>());
        }
        fluid = new FluidPoint(mNr, mIsAxial, mCoords, mass);
    }
    
    // released as different point classes
    if (isSolid && isFluid) {
        SFCoupling *couple;
        if (XMath::equalRows(mSFNormal_assmble) && XMath::equalRows(mMassFluid)) {
            double ns = mSFNormal_assmble(0, 0);
            double np = mSFNormal_assmble(0, 1);
            double nz = mSFNormal_assmble(0, 2);
            // in 1D cases, the normal vector should be in-slice 
            if (std::abs(np) > tinyDouble * sqrt(ns * ns + nz * nz)) {
                throw std::runtime_error("GLLPoint::release || Invalid 1D solid-fluid boundary normal.");
            }
            couple = new SFCoupling1D((Real)mSFNormal(0, 0), (Real)mSFNormal(0, 2), 
                (Real)(mSFNormal_assmble(0, 0) / mMassFluid(0)), 
                (Real)(mSFNormal_assmble(0, 2) / mMassFluid(0)));
        } else {
            RDMatX3 normal_invmf = mSFNormal_assmble;
            normal_invmf.col(0).array() /= mMassFluid.array();
            normal_invmf.col(1).array() /= mMassFluid.array();
            normal_invmf.col(2).array() /= mMassFluid.array();
            couple = new SFCoupling3D(mSFNormal.cast<Real>(), normal_invmf.cast<Real>());
        }
        SolidFluidPoint *sfpoint = new SolidFluidPoint(solid, fluid, couple);
        domain.addSFPoint(sfpoint);
        return domain.addPoint(sfpoint);
    } else if (isSolid) {
        return domain.addPoint(solid);         
    } else if (isFluid) {
        return domain.addPoint(fluid);         
    } else {
        throw std::runtime_error("GLLPoint::release || Point is in neither solid nor fluid domain.");
    } 
}

void GLLPoint::feedBuffer(RDMatXX &buffer, int col) {
    int nr = mSFNormal_assmble.rows();
    buffer.block(0, col, nr, 1) = mMassSolid;
    buffer.block(nr, col, nr, 1) = mMassFluid;
    buffer.block(nr * 2, col, nr * 3, 1) = Eigen::Map<RDColX>(mSFNormal_assmble.data(), nr * 3);
    buffer.block(nr * 5, col, nr * 3, 1) = Eigen::Map<RDColX>(mSurfNormal.data(), nr * 3);
    buffer(nr * 8, col) = (double)mReferenceCount;
}

void GLLPoint::extractBuffer(RDMatXX &buffer, int col) {
    int nr = mSFNormal_assmble.rows();
    mMassSolid += buffer.block(0, col, nr, 1);
    mMassFluid += buffer.block(nr, col, nr, 1);
    mSFNormal_assmble += Eigen::Map<RDMatX3>(buffer.block(nr * 2, col, nr * 3, 1).data(), nr, 3);
    mSurfNormal += Eigen::Map<RDMatX3>(buffer.block(nr * 5, col, nr * 3, 1).data(), nr, 3);
    mReferenceCount += round(buffer(nr * 8, col));
}



