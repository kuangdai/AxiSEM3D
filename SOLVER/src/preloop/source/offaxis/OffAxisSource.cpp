// OffAxisSource.h
// created by Kuangdai on 11-Nov-2017 
// base class of off-axis source

#include "OffAxisSource.h"
#include "Geodesy.h"
#include "Quad.h"
#include "Domain.h"
#include "Element.h"
// TODO: implement this in core
// #include "OffAxisSourceTerm.h"
#include "Mesh.h"
#include "XMath.h"
#include "SpectralConstants.h"
#include "XMPI.h"
#include "MultilevelTimer.h"

OffAxisSource::OffAxisSource(double depth, double lat, double lon,
    double srcLat, double srcLon, double srcDep):
mDepth(depth), mLatitude(lat), mLongitude(lon) {
    // handle singularity at poles
    if (std::abs(mLatitude - 90.) < tinyDouble) {
        mLatitude = 90.;
        mLongitude = 0.;
    }
    if (std::abs(mLatitude + 90.) < tinyDouble) {
        mLatitude = -90.;
        mLongitude = 0.;
    }
    // compute theta and phi in source-centered coordinate system
    RDCol3 rtpG, rtpS;
    rtpG(0) = 1.;
    rtpG(1) = Geodesy::lat2Theta_d(mLatitude, mDepth);
    rtpG(2) = Geodesy::lon2Phi(mLongitude);
    rtpS = Geodesy::rotateGlob2Src(rtpG, srcLat, srcLon, srcDep);
    mThetaSrc = rtpS(1);
    mPhiSrc = rtpS(2);
}

void OffAxisSource::release(Domain &domain, const Mesh &mesh) const {
    MultilevelTimer::begin("Locate Off-axis Source", 2);
    // locate local
    int myrank = XMPI::nproc();
    int locTag;
    RDColP interpFactXii, interpFactEta;
    if (locate(mesh, locTag, interpFactXii, interpFactEta)) {
        myrank = XMPI::rank();
    }

    // min recRank
    int myrank_min = XMPI::min(myrank);
    if (myrank_min == XMPI::nproc()) {
        throw std::runtime_error("OffAxisSource::release || Error locating off-axis source.");
    }
    MultilevelTimer::end("Locate Off-axis Source", 2);

    MultilevelTimer::begin("Compute Off-axis Source", 2);
    // release to me
    if (myrank_min == XMPI::rank()) {
        // compute OffAxisSource term
        vec_arPP_CMatX3 fouriers;
        const Quad *myQuad = mesh.getQuad(locTag);
        computeSourceFourier(*myQuad, interpFactXii, interpFactEta, mPhiSrc, fouriers);
        // add to domain
        Element *myElem = domain.getElement(myQuad->getElementTag());
        // TODO: implement this in core
        //       the existing class SourceTerm is not general enough
        //       because the source time funciton varies for each source
        // domain.addOffAxisSourceTerm(new OffAxisSourceTerm(myElem, fouriers));
    }
    MultilevelTimer::end("Compute Off-axis Source", 2);
}

bool OffAxisSource::locate(const Mesh &mesh, int &locTag, 
    RDColP &interpFactXii, RDColP &interpFactEta) const {
    MultilevelTimer::begin("R OffAxisSource", 3);
    RDCol2 srcCrds = RDCol2::Zero();
    double r = mesh.computeRadiusRef(mDepth, mLatitude, mLongitude);
    srcCrds(0) = r * sin(mThetaSrc);
    srcCrds(1) = r * cos(mThetaSrc);
    MultilevelTimer::end("R OffAxisSource", 3);

    // check range of subdomain
    if (srcCrds(0) > mesh.sMax() + tinySingle || srcCrds(0) < mesh.sMin() - tinySingle) {
        return false;
    }
    if (srcCrds(1) > mesh.zMax() + tinySingle || srcCrds(1) < mesh.zMin() - tinySingle) {
        return false;
    }
    // find host element
    RDCol2 srcXiEta;
    for (int iloc = 0; iloc < mesh.getNumQuads(); iloc++) {
        const Quad *quad = mesh.getQuad(iloc);
        if (quad->isFluid() || !quad->nearMe(srcCrds(0), srcCrds(1))) {
            continue;
        }
        if (quad->invMapping(srcCrds, srcXiEta)) {
            if (std::abs(srcXiEta(0)) <= 1.000001 && std::abs(srcXiEta(1)) <= 1.000001) {
                locTag = iloc;
                XMath::interpLagrange(srcXiEta(0), nPntEdge,
                    (quad->isAxial() ? SpectralConstants::getP_GLJ().data()
                                     : SpectralConstants::getP_GLL().data()), 
                    interpFactXii.data());
                XMath::interpLagrange(srcXiEta(1), nPntEdge,
                    SpectralConstants::getP_GLL().data(), interpFactEta.data());
                return true;
            }
        }
    }
    return false;
}

#include "Parameters.h"
// #include "OffAxisPointForce.h"
#include <fstream>
#include <boost/algorithm/string.hpp>

void OffAxisSource::buildInparam(std::vector<OffAxisSource> *&offsrc, 
    const Parameters &par, int verbose) {
    // TODO: 
    // read location of off-axis sources
    // read source-time function
    // build off-axis source objects
}
