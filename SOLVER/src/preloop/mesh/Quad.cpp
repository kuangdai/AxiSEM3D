// Quad.cpp
// created by Kuangdai on 5-May-2016 
// general Quad element 
// created from Exodus model, elements of AxiSEM3D mesh

#include "Quad.h"
#include "ExodusModel.h"

#include "SphericalMapping.h"
#include "LinearMapping.h"
#include "SemiSphericalMapping.h"
#include "SpectralConstants.h"

#include "Material.h"
#include "Relabelling.h"

#include "Elastic.h"
#include "Acoustic.h"

#include "Gradient.h"
#include "PreloopGradient.h"
#include "SolidElement.h"
#include "FluidElement.h"
#include "Domain.h"

#include "GLLPoint.h" 
#include "NrField.h"

#include "XMath.h"
#include "Geodesy.h"

#include "OceanLoad3D.h"

#include "PreloopFFTW.h"
#include <cfloat>

Quad::Quad(const ExodusModel &exModel, int quadTag, const NrField &nrf): 
mQuadTag(quadTag) {
    
    // connectivity and coords
    const IMatX4 &connect = exModel.getConnectivity();
    for (int i = 0; i < 4; i++) {
        int nodeTag = connect(mQuadTag, i);
        mNodalCoords(0, i) = exModel.getNodalS(nodeTag);
        mNodalCoords(1, i) = exModel.getNodalZ(nodeTag);
        mGlobalNodeTags(i) = nodeTag;
    }
    
    // geometric mapping
    double distTol = exModel.getDistTolerance();
    double etype = exModel.getElementalVariables("element_type", mQuadTag);
    if (etype < .5) {
        // etype = 0.0, spherical
        mMapping = new SphericalMapping();
        double r0 = mNodalCoords.col(0).norm();
        double r1 = mNodalCoords.col(1).norm();
        double r2 = mNodalCoords.col(2).norm();
        double r3 = mNodalCoords.col(3).norm();
        if (std::abs(r0 - r1) < distTol && std::abs(r2 - r3) < distTol) {
            mCurvedOuter = r2 > r0 ? 2 : 0;
        } else if (std::abs(r1 - r2) < distTol && std::abs(r3 - r0) < distTol) {
            mCurvedOuter = r3 > r1 ? 3 : 1;       
        } else {
            throw std::runtime_error("Quad::Quad || Invalid spherical element shape.");
        }
    } else if (etype < 1.5) {
        // etype = 1.0, linear
        mMapping = new LinearMapping();
        mCurvedOuter = -1;
    } else {
        // etype = 2.0, semi-spherical
        mMapping = new SemiSphericalMapping();
        double r0 = mNodalCoords.col(0).norm();
        double r1 = mNodalCoords.col(1).norm();
        double r2 = mNodalCoords.col(2).norm();
        double r3 = mNodalCoords.col(3).norm();
        if (std::abs(r0 - r1) < distTol && r0 > r2) {
            mCurvedOuter = 0;
        } else if (std::abs(r1 - r2) < distTol && r1 > r3) {
            mCurvedOuter = 1;
        } else if (std::abs(r2 - r3) < distTol && r2 > r0) {
            mCurvedOuter = 2;
        } else if (std::abs(r3 - r0) < distTol && r3 > r1) {
            mCurvedOuter = 3;
        } else {
            throw std::runtime_error("Quad::Quad || Invalid semi-spherical element shape.");
        }
    }
    
    // solid fluid
    mIsFluid = exModel.getElementalVariables("fluid", mQuadTag) > .5;
    
    // axial boundary
    mIsAxial = false;
    mAxialSide = exModel.getSideAxis(mQuadTag);
    if (mAxialSide >= 0) {
        mIsAxial = true;
        if (mMapping->getType() != Mapping::MappingTypes::Linear && (mAxialSide + mCurvedOuter) % 2 == 0) {
            throw std::runtime_error("Quad::Quad || Conflict in axial setting.");
        }
        if (mAxialSide != 3) {
            throw std::runtime_error("Quad::Quad || Axial side must be 3.");
        }
    }

    // solid-fluid boundary
    ///////////////////////////////////////////////
    ////// this piece of code is correct and more general,
    ////// which should be re-enabled after fixing the issue in mesher
    mOnSFBoundary = false;
    mSFSide = exModel.getSideSolidFluid(mQuadTag);
    if (mSFSide >= 0) {
        mOnSFBoundary = true;
        if (!exModel.isCartesian()) {
            if (mMapping->getType() == Mapping::MappingTypes::Linear) {
                throw std::runtime_error("Quad::Quad || Conflict in solid-fluid boundary.");
            }
            if (mMapping->getType() == Mapping::MappingTypes::Spherical && (mSFSide + mCurvedOuter) % 2 != 0) {
                throw std::runtime_error("Quad::Quad || Conflict in solid-fluid boundary.");
            }
            if (mMapping->getType() == Mapping::MappingTypes::SemiSpherical && mSFSide != mCurvedOuter) {
                throw std::runtime_error("Quad::Quad || Conflict in solid-fluid boundary.");
            }
        }
    }
    ///////////////////////////////////////////////
    
    ///////////////////////////////////////////////
    ////// this piece of code is temporary
    // change radius of CMB and ICB if necessary
    // double rCMB = exModel.getR_CMB();
    // double rICB = exModel.getR_ICB();
    // 
    // // init
    // mOnSFBoundary = false;
    // mSFSide = -1;
    // 
    // // vertex radius
    // double r0 = mNodalCoords.col(0).norm();
    // double r1 = mNodalCoords.col(1).norm();
    // double r2 = mNodalCoords.col(2).norm();
    // double r3 = mNodalCoords.col(3).norm();
    // 
    // // CMB
    // if (std::abs(r0 - rCMB) < distTol && std::abs(r1 - rCMB) < distTol) {
    //     mOnSFBoundary = true;
    //     mSFSide = 0;
    // }
    // if (std::abs(r1 - rCMB) < distTol && std::abs(r2 - rCMB) < distTol) {
    //     mOnSFBoundary = true;
    //     mSFSide = 1;
    // }
    // if (std::abs(r2 - rCMB) < distTol && std::abs(r3 - rCMB) < distTol) {
    //     mOnSFBoundary = true;
    //     mSFSide = 2;
    // }
    // if (std::abs(r3 - rCMB) < distTol && std::abs(r0 - rCMB) < distTol) {
    //     mOnSFBoundary = true;
    //     mSFSide = 3;
    // }
    // 
    // // ICB
    // if (std::abs(r0 - rICB) < distTol && std::abs(r1 - rICB) < distTol) {
    //     mOnSFBoundary = true;
    //     mSFSide = 0;
    // }
    // if (std::abs(r1 - rICB) < distTol && std::abs(r2 - rICB) < distTol) {
    //     mOnSFBoundary = true;
    //     mSFSide = 1;
    // }
    // if (std::abs(r2 - rICB) < distTol && std::abs(r3 - rICB) < distTol) {
    //     mOnSFBoundary = true;
    //     mSFSide = 2;
    // }
    // if (std::abs(r3 - rICB) < distTol && std::abs(r0 - rICB) < distTol) {
    //     mOnSFBoundary = true;
    //     mSFSide = 3;
    // }
    
    
    ///////////////////////////////////////////////


    // surface boundary
    mOnSurface = false;
    mSurfaceSide = exModel.getSideSurface(mQuadTag);
    if (mSurfaceSide >= 0) {
        mOnSurface = true;
        if (!exModel.isCartesian()) {
            if (mMapping->getType() == Mapping::MappingTypes::Linear) {
                throw std::runtime_error("Quad::Quad || Conflict in surface setting.");
            }
            if (mMapping->getType() == Mapping::MappingTypes::Spherical && mSurfaceSide != mCurvedOuter) {
                throw std::runtime_error("Quad::Quad || Conflict in surface setting.");
            }
            if (mMapping->getType() == Mapping::MappingTypes::SemiSpherical && mSurfaceSide != mCurvedOuter) {
                throw std::runtime_error("Quad::Quad || Conflict in surface setting.");
            }
        }
        // if (mIsFluid) {
        //     throw std::runtime_error("Quad::Quad || Fluid element on surface. Not implemented.");
        // }
        // if (mOnSFBoundary) {
        //     throw std::runtime_error("Quad::Quad || Element on both surface and solid-fluid boundary.");
        // }
    }
    
    // nr field 
    for (int j = 0; j < 4; j++) {
        mNearAxisNodes(j) = exModel.getVicinalAxis(quadTag)[j];
        mNodalAveGLLSpacing(j) = exModel.getAveGLLSpacing(mGlobalNodeTags(j));
    } 
    formNrField(nrf, distTol);
    
    // integral factor
    formIntegralFactor();    
    
    // 1D material
    mMaterial = new Material(this, exModel);
    
    // relabelling
    // mRelabelling = new Relabelling(this);
    mRelabelling = 0;
    
    // dt 
    mDeltaTRef = exModel.getElementalVariables("dt", mQuadTag);
    
    // init ocean depth
    for (int ipol = 0; ipol <= nPol; ipol++) {
        for (int jpol = 0; jpol <= nPol; jpol++) {
            int ipnt = ipol * nPntEdge + jpol;
            mOceanDepth[ipnt] = RDColX::Zero(mPointNr(ipol, jpol));    
        }
    }
}

Quad::~Quad() {
    delete mMapping;
    delete mMaterial;
    if (mRelabelling) {
        delete mRelabelling;
    }
}

void Quad::addVolumetric3D(const std::vector<Volumetric3D *> &m3D, 
    double srcLat, double srcLon, double srcDep, double phi2D) {
    mMaterial->addVolumetric3D(m3D, srcLat, srcLon, srcDep, phi2D);
}

void Quad::addGeometric3D(const std::vector<Geometric3D *> &g3D, 
    double srcLat, double srcLon, double srcDep, double phi2D) {
    if (g3D.size() > 0) {
        mRelabelling = new Relabelling(this);
        mRelabelling->addUndulation(g3D, srcLat, srcLon, srcDep, phi2D);
    }    
}

void Quad::setOceanLoad3D(const OceanLoad3D &o3D, 
    double srcLat, double srcLon, double srcDep, double phi2D) {
    if (!mOnSurface) {
        return;
    }
    for (int ipol = 0; ipol <= nPol; ipol++) {
        for (int jpol = 0; jpol <= nPol; jpol++) {
            bool surface = mOnSurface && (
                (mSurfaceSide == 0 && jpol == 0) || 
                (mSurfaceSide == 1 && ipol == nPol) || 
                (mSurfaceSide == 2 && jpol == nPol) ||
                (mSurfaceSide == 3 && ipol == 0));
            if (surface) {
                int nr_read = getPointNr(ipol, jpol);
                int ipnt = ipol * nPntEdge + jpol;
                const RDCol2 &xieta = SpectralConstants::getXiEta(ipol, jpol, mIsAxial);
                const RDMatX3 &rtpS = computeGeocentricGlobal(srcLat, srcLon, srcDep, xieta, nr_read, phi2D);
                for (int alpha = 0; alpha < nr_read; alpha++) {
                    double t = rtpS(alpha, 1);
                    double p = rtpS(alpha, 2);
                    mOceanDepth[ipnt](alpha) = o3D.getOceanDepth(t, p); 
                }
            }    
        }
    }
}

bool Quad::hasRelabelling() const {
    if (!mRelabelling) {
        return false;
    }
    return !(mRelabelling->isZero());
}

double Quad::getDeltaT() const {
    // courant number
    double courant = getCourant();
    
    // vmax on slices
    RDColX vmax = mMaterial->getVMax();
    
    // hmin on slices
    RDColX hmin = getHminSlices();
    
    // dt on slices
    RDColX dt_slices = courant * hmin.schur(vmax.array().pow(-1.).matrix());
    double dt_min_org = dt_slices.minCoeff();
    
    // we have checked this in relabelling
    return dt_min_org;
    
    // int denseFactor = 11;
    // // check minimum on densed sampling slices
    // double dt_min_all = XMath::trigonResampling(denseFactor * dt_slices.rows(), dt_slices).minCoeff();
    // if (dt_min_all < .0) {
    //     throw std::runtime_error("Quad::getDeltaT || Strong Gibbs phenomenon detected. ||"
    //         "This exception means your 3-D models (usually of Geometric type, e.g., Moho undulation) ||"
    //         "contain very strong and sharp gradients such as a spike. The Gibbs phenomenon is too ||"
    //         "strong to be handled by AxiSEM3D and will definitely cause instability. ||"
    //         "Please try to smooth the problematic model externally.");
    // }
    // 
    // // return average of the two
    // return (dt_min_org + dt_min_all) / 2.;
}

void Quad::setupGLLPoints(std::vector<GLLPoint *> &gllPoints, const IMatPP &myPointTags, double distTol) {
    // compute mass on points
    const arPP_RDColX &mass = mMaterial->computeElementalMass();

    for (int ipol = 0; ipol <= nPol; ipol++) {
        for (int jpol = 0; jpol <= nPol; jpol++) {
            int ipnt = ipol * nPntEdge + jpol;
            // get point tag
            int pointTag = myPointTags(ipol, jpol);
            
            // setup properties
            bool axial = mIsAxial && ipol == 0;
            bool surface = mOnSurface && (
                (mSurfaceSide == 0 && jpol == 0) || 
                (mSurfaceSide == 1 && ipol == nPol) || 
                (mSurfaceSide == 2 && jpol == nPol) ||
                (mSurfaceSide == 3 && ipol == 0));
            const RDCol2 &xieta = SpectralConstants::getXiEta(ipol, jpol, mIsAxial);
            const RDCol2 &crds = mapping(xieta);
            gllPoints[pointTag]->setup(mPointNr(ipol, jpol), axial, surface, crds, distTol);
            
            // add mass 
            if (mIsFluid) {
                gllPoints[pointTag]->addMassFluid(mass[ipnt]);
            } else {
                gllPoints[pointTag]->addMassSolid(mass[ipnt]);
            }
            
            //////// boundary terms ////////
            bool sfbry = mOnSFBoundary && (
                (mSFSide == 0 && jpol == 0) || 
                (mSFSide == 1 && ipol == nPol) || 
                (mSFSide == 2 && jpol == nPol) ||
                (mSFSide == 3 && ipol == 0));
            
            // normal
            if (sfbry) {
                RDMatX3 normal = computeNormal(mSFSide, ipol, jpol);
                // inverse in solid domain
                if (!mIsFluid) normal *= -1.; 
                // both solid and fluid elements contribute, so divided by 2
                gllPoints[pointTag]->addSFNormal(normal * .5);
            }
            
            if (surface) {
                // surface area
                const RDMatX3 &normal = computeNormal(mSurfaceSide, ipol, jpol);
                gllPoints[pointTag]->addSurfNormal(normal);
                gllPoints[pointTag]->setOceanDepth(mOceanDepth[ipnt]);    
            }
        }
    }
}

int Quad::release(Domain &domain, const IMatPP &myPointTags, const AttBuilder *attBuild) const {
    if (mIsFluid) {
        return releaseFluid(domain, myPointTags);
    } else {
        return releaseSolid(domain, myPointTags, attBuild);
    } 
}

int Quad::releaseSolid(Domain &domain, const IMatPP &myPointTags, const AttBuilder *attBuild) const {
    bool elem1D = mMaterial->isSolidPar1D(attBuild != 0);
    if (hasRelabelling()) {
        elem1D = elem1D && mRelabelling->isPar1D();
    }
    std::array<Point *, nPntElem> points;
    for (int ipol = 0; ipol <= nPol; ipol++) {
        for (int jpol = 0; jpol <= nPol; jpol++) {
            points[ipol * nPntEdge + jpol] = domain.getPoint(myPointTags(ipol, jpol));
        }
    }
    Gradient *grad = createGraident();
    PRT *prt = mRelabelling ? mRelabelling->createPRT(elem1D) : 0;
    Elastic *elas = mMaterial->createElastic(elem1D, attBuild);
    Element *elem = new SolidElement(grad, prt, points, elas);
    return domain.addElement(elem);
}

int Quad::releaseFluid(Domain &domain, const IMatPP &myPointTags) const {
    bool elem1D = mMaterial->isFluidPar1D();
    if (hasRelabelling()) {
        elem1D = elem1D && mRelabelling->isPar1D();
    }
    std::array<Point *, nPntElem> points;
    for (int ipol = 0; ipol <= nPol; ipol++) {
        for (int jpol = 0; jpol <= nPol; jpol++) {
            points[ipol * nPntEdge + jpol] = domain.getPoint(myPointTags(ipol, jpol));
        }
    }
    Gradient *grad = createGraident();
    PRT *prt = mRelabelling ? mRelabelling->createPRT(elem1D) : 0;
    Acoustic *acous = mMaterial->createAcoustic(elem1D);
    Element *elem = new FluidElement(grad, prt, points, acous);
    return domain.addElement(elem);
}

RDCol2 Quad::mapping(const RDCol2 &xieta) const {
    return mMapping->mapping(mNodalCoords, xieta, mCurvedOuter);
}

RDMat22 Quad::jacobian(const RDCol2 &xieta) const {
    return mMapping->jacobian(mNodalCoords, xieta, mCurvedOuter);
}

double Quad::detJacobian(const RDCol2 &xieta) const {
    return mMapping->detJacobian(mNodalCoords, xieta, mCurvedOuter);
}

bool Quad::invMapping(const RDCol2 &sz, RDCol2 &xieta) const {
    return mMapping->invMapping(mNodalCoords, sz, mCurvedOuter, xieta);
}

RDRow4 Quad::computeWeightsCG4() const {
    RDRow4 weights_cg4;
    RDMatPP ifact;
    XMath::structuredUseFirstRow(mIntegralFactor, ifact);
    weights_cg4(0) = (ifact(0, 0) + ifact(0, 1)
                   + ifact(1, 0) + ifact(1, 1)
            + 0.5 * (ifact(0, 2) + ifact(1, 2)
                   + ifact(2, 0) + ifact(2, 1))
            + 0.25 * ifact(2, 2)) / ifact(1, 1);
            
    weights_cg4(1) = (ifact(0, 3) + ifact(0, 4)
                   + ifact(1, 3) + ifact(1, 4)
            + 0.5 * (ifact(0, 2) + ifact(1, 2)
                   + ifact(2, 3) + ifact(2, 4))
            + 0.25 * ifact(2, 2)) / ifact(1, 3);
            
    weights_cg4(2) = (ifact(3, 0) + ifact(3, 1)
                   + ifact(4, 0) + ifact(4, 1)
            + 0.5 * (ifact(2, 0) + ifact(2, 1)
                   + ifact(3, 2) + ifact(4, 2))
            + 0.25 * ifact(2, 2)) / ifact(3, 1);
            
    weights_cg4(3) = (ifact(3, 3) + ifact(3, 4)
                   + ifact(4, 3) + ifact(4, 4)
            + 0.5 * (ifact(2, 3) + ifact(2, 4)
                   + ifact(3, 2) + ifact(4, 2))
            + 0.25 * ifact(2, 2)) / ifact(3, 3);            
    return weights_cg4;
}

RDMatX3 Quad::computeGeocentricGlobal(double srcLat, double srcLon, double srcDep,
    const RDCol2 &xieta, int npnt, double phi2D) const {
    RDMatX3 rtpG_Nr(npnt, 3);
    RDCol3 rtpS;
    Geodesy::rtheta(mapping(xieta), rtpS(0), rtpS(1));
    // debug relabelling
    // rtpS(1) = Geodesy::theta(mapping(RDCol2::Zero()));
    double dphi = 2. * pi / npnt;
    for (int i = 0; i < npnt; i++) {
        rtpS(2) = phi2D < -DBL_MAX * .9 ? dphi * i : phi2D;
        rtpG_Nr.row(i) = Geodesy::rotateSrc2Glob(rtpS, srcLat, srcLon, srcDep).transpose();
    }
    return rtpG_Nr;
}

double Quad::computeCenterRadius() const {
    // xi = eta = 0
    return mapping(RDCol2::Zero()).norm();
}

void Quad::computeGradientScalar(const vec_CDMatPP &u, vec_ar3_CDMatPP &u_i) const {
    // create Gradient object
    RDMatPP dsdxii, dsdeta, dzdxii, dzdeta, inv_s;
    for (int ipol = 0; ipol <= nPol; ipol++) {
        for (int jpol = 0; jpol <= nPol; jpol++) {
            const RDCol2 &xieta = SpectralConstants::getXiEta(ipol, jpol, mIsAxial);
            const RDMat22 &J = jacobian(xieta);
            double detJ = J.determinant();
            dsdxii(ipol, jpol) = J(0, 0) / detJ;
            dsdeta(ipol, jpol) = -J(0, 1) / detJ;
            dzdxii(ipol, jpol) = -J(1, 0) / detJ;
            dzdeta(ipol, jpol) = J(1, 1) / detJ;
            double s = mapping(xieta)(0);
            inv_s(ipol, jpol) = (mIsAxial && ipol == 0) ? 0. : 1. / s;
        }
    }
    PreloopGradient grad(dsdxii, dsdeta, dzdxii, dzdeta, inv_s, mIsAxial);
    // compute gradient
    grad.gradScalar(u, u_i, mNr / 2, mNr % 2 == 0);
}

Gradient *Quad::createGraident() const {
    RDMatPP dsdxii, dsdeta, dzdxii, dzdeta, inv_s;
    for (int ipol = 0; ipol <= nPol; ipol++) {
        for (int jpol = 0; jpol <= nPol; jpol++) {
            const RDCol2 &xieta = SpectralConstants::getXiEta(ipol, jpol, mIsAxial);
            const RDMat22 &J = jacobian(xieta);
            double detJ = J.determinant();
            dsdxii(ipol, jpol) = J(0, 0) / detJ;
            dsdeta(ipol, jpol) = -J(0, 1) / detJ;
            dzdxii(ipol, jpol) = -J(1, 0) / detJ;
            dzdeta(ipol, jpol) = J(1, 1) / detJ;
            double s = mapping(xieta)(0);
            inv_s(ipol, jpol) = (mIsAxial && ipol == 0) ? 0. : 1. / s;
        }
    }
    return new Gradient(dsdxii, dsdeta, dzdxii, dzdeta, inv_s, mIsAxial);  
}

void Quad::formIntegralFactor() {
    for (int ipol = 0; ipol <= nPol; ipol++) {
        for (int jpol = 0; jpol <= nPol; jpol++) {
            const RDCol2 &xieta = SpectralConstants::getXiEta(ipol, jpol, mIsAxial);
            const RDCol2 &weights = SpectralConstants::getWeights(ipol, jpol, mIsAxial);
            double wxi_weta = weights(0) * weights(1);
            double detJ = detJacobian(xieta);
            double s = mapping(xieta)(0);
            if (mIsAxial) {
                if (ipol == 0) {
                    const RDMat22 &J = jacobian(xieta);
                    mIntegralFactor(ipol * nPntEdge + jpol) = wxi_weta * J(0, 0) * detJ; 
                } else {
                    mIntegralFactor(ipol * nPntEdge + jpol) = wxi_weta * s / (1. + xieta(0)) * detJ; 
                }
            } else {
                mIntegralFactor(ipol * nPntEdge + jpol) = wxi_weta * s * detJ;  
            }
        }
    }
}

void Quad::formNrField(const NrField &nrf, double distTol) {
    for (int ipol = 0; ipol <= nPol; ipol++) {
        for (int jpol = 0; jpol <= nPol; jpol++) {
            // interpolate
            const RDCol2 &xieta = SpectralConstants::getXiEta(ipol, jpol, mIsAxial);
            RDCol2 crds = mapping(xieta);
            // erase numerical error
            double offset = std::min((double)tinySingle, distTol / 1e3);
            crds(0) = round(crds(0) / offset) * offset;
            crds(1) = round(crds(1) / offset) * offset;
            mPointNr(ipol, jpol) = nrf.getNrAtPoint(crds);
            
            // upper limit
            double spacing = Mapping::interpolate(mNodalAveGLLSpacing, xieta);
            double circ = 2. * pi * crds(0);
            int upper = (int)(circ / spacing);
            if (upper < 3) upper = 3; // axis
            mPointNr(ipol, jpol) = std::min(mPointNr(ipol, jpol), upper);
            
            ////////// deal with even numbers ////////
            bool forceOdd = false;
            if (mPointNr(ipol, jpol) % 2 == 0) {
                if (mIsAxial) {
                    // axis
                    mPointNr(ipol, jpol)++;
                    forceOdd = true;
                } else if (mNearAxisNodes.maxCoeff() >= 0) {
                    // near axis 
                    for (int i = 0; i < 4; i++) {
                        int n0 = mNearAxisNodes(i);
                        int n1 = mNearAxisNodes(Mapping::period0123(i + 1));
                        if (n0 >= 0 && n1 >= 0) {
                            // on side
                            bool onSide = 
                            (n0 == 0 && jpol == 0) || 
                            (n0 == 1 && ipol == nPol) || 
                            (n0 == 2 && jpol == nPol) ||
                            (n0 == 3 && ipol == 0);
                            if (onSide) {
                                mPointNr(ipol, jpol)++;
                                forceOdd = true;
                                break;
                            }
                        } else if (n0 >= 0 && n1 < 0) {
                            // on point
                            bool onPoint = 
                            (n0 == 0 && ipol == 0 && jpol == 0) || 
                            (n0 == 1 && ipol == nPol && jpol == 0) || 
                            (n0 == 2 && ipol == nPol && jpol == nPol) ||
                            (n0 == 3 && ipol == 0 && jpol == nPol);
                            if (onPoint) {
                                mPointNr(ipol, jpol)++;
                                forceOdd = true;
                                break;
                            }
                        }
                    }
                }
            }
            
            // luck number
            if (nrf.useLuckyNumber()) {
                mPointNr(ipol, jpol) = PreloopFFTW::nextLuckyNumber(mPointNr(ipol, jpol), forceOdd);
            }
        }
    }
    mNr = mPointNr.maxCoeff();
}

double Quad::getCourant() const {
    // 1D reference
    double vmaxRef = mMaterial->getVMaxRef();
    double hminRef = DBL_MAX;
    for (int i = 0; i < 4; i++) {
        double sideLen = (mNodalCoords.col(i) - mNodalCoords.col(Mapping::period0123(i + 1))).norm();
        hminRef = std::min(hminRef, sideLen);
    }
    // dt = courant * hmin / vmax
    return mDeltaTRef * vmaxRef / hminRef;
}

RDColX Quad::getHminSlices() const {
    // point coordinates
    std::vector<RDCol2> coordsRef;
    for (int ipol = 0; ipol <= nPol; ipol++) {
        for (int jpol = 0; jpol <= nPol; jpol++) {
            const RDCol2 &xieta = SpectralConstants::getXiEta(ipol, jpol, mIsAxial);
            coordsRef.push_back(mapping(xieta));
        }
    }
    // hmin of reference mesh
    int nslices = getNr();
    RDColX hmin = RDColX::Constant(nslices, XMath::findClosestDist(coordsRef));
    // hmin of deformed mesh
    if (hasRelabelling()) {
        const RDMatXN &deltaR = mRelabelling->getDeltaR();
        for (int islice = 0; islice < nslices; islice++) {
            std::vector<RDCol2> coordsSlice = coordsRef;
            for (int ipol = 0; ipol <= nPol; ipol++) {
                for (int jpol = 0; jpol <= nPol; jpol++) {
                    int ipnt = ipol * nPntEdge + jpol;
                    double theta = Geodesy::theta(coordsRef[ipnt]);
                    coordsSlice[ipnt](0) += deltaR(islice, ipnt) * sin(theta);
                    coordsSlice[ipnt](1) += deltaR(islice, ipnt) * cos(theta);
                }
            }
            hmin(islice) = XMath::findClosestDist(coordsSlice);
        }
    }
    return hmin;
}

RDMatX3 Quad::computeNormal(int side, int ipol, int jpol) const {
    // check element-side type
    if (mMapping->getType() == Mapping::MappingTypes::Spherical && (side + mCurvedOuter) % 2 != 0) {
        throw std::runtime_error("Quad::computeNormal || Computing normal on non-spherical side.");
    }
    if (mMapping->getType() == Mapping::MappingTypes::SemiSpherical && side != mCurvedOuter) {
        throw std::runtime_error("Quad::computeNormal || Computing normal on non-spherical side.");
    }
    
    // compute common boundary terms
    double r0, r1, theta0, theta1;
    Geodesy::rtheta(mNodalCoords.col(side), r0, theta0);
    Geodesy::rtheta(mNodalCoords.col(Mapping::period0123(side + 1)), r1, theta1);
    double rsf = .5 * (r0 + r1);
    double half_r_dtheta = .5 * rsf * std::abs(theta1 - theta0);
    double half_r2_dtheta = .5 * rsf * rsf * std::abs(theta1 - theta0);  
    
    // compute normal of undulated solid-fluid boundary
    const RDCol2 &xieta = SpectralConstants::getXiEta(ipol, jpol, mIsAxial);
    const RDCol2 &crds = mapping(xieta);
    double sint = crds(0) / crds.norm();
    double cost = crds(1) / crds.norm();
    RDMat33 Q = RDMat33::Zero();
    Q(0, 0) = cost;
    Q(0, 2) = sint;
    Q(1, 1) = 1.;
    Q(2, 0) = -sint;
    Q(2, 2) = cost;
    RDMatX3 nRTZ(mPointNr(ipol, jpol), 3);
    if (hasRelabelling()) {
        nRTZ = mRelabelling->getSFNormalRTZ(ipol, jpol);
    } else {
        nRTZ.col(0).fill(0.);
        nRTZ.col(1).fill(0.);
        nRTZ.col(2).fill(1.);
    }
    RDMatX3 normal = nRTZ * Q.transpose();
    // multiply area
    const RDCol2 &weights = SpectralConstants::getWeights(ipol, jpol, mIsAxial);
    double wsf = (side == 0 || side == 2) ? weights(0) : weights(1);
    if (mIsAxial) {
        if (ipol == 0) {
            const RDMat22 &J = jacobian(xieta);
            normal *= wsf * J(0, 0) * half_r_dtheta;
        } else {
            normal *= wsf / (1. + xieta(0)) * sint * half_r2_dtheta;
        }
    } else {
        normal *= wsf * sint * half_r2_dtheta;
    }
    if (side != mCurvedOuter) {
        normal *= -1.;
    }
    return normal;
}

RDRowN Quad::getUndulationOnSlice(double phi) const {
    if (hasRelabelling()) {
        return XMath::computeFourierAtPhi(mRelabelling->getDeltaR(), phi);
    } else {
        return RDRowN::Zero();
    }
}

RDRowN Quad::getMaterialOnSlice(const std::string &parName, int refType, double phi) const {
    return XMath::computeFourierAtPhi(mMaterial->getProperty(parName, refType), phi);
}

void Quad::getSpatialRange(double &s_max, double &s_min, double &z_max, double &z_min) const {
    s_max = mNodalCoords.row(0).maxCoeff();
    s_min = mNodalCoords.row(0).minCoeff();
    z_max = mNodalCoords.row(1).maxCoeff();
    z_min = mNodalCoords.row(1).minCoeff();
}

bool Quad::nearMe(double s, double z) const {
    if (s > mNodalCoords.row(0).maxCoeff() + tinySingle || s < mNodalCoords.row(0).minCoeff() - tinySingle) {
        return false;
    }
    if (z > mNodalCoords.row(1).maxCoeff() + tinySingle || z < mNodalCoords.row(1).minCoeff() - tinySingle) {
        return false;
    }
    return true;
}

