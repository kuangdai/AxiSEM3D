// Receiver.cpp
// created by Kuangdai on 13-May-2016 
// receiver   

#include "Receiver.h"
#include "XMath.h"
#include "Geodesy.h"
#include "Domain.h"
#include "Mesh.h"
#include "Quad.h"
#include "Element.h"
#include "SpectralConstants.h"
#include "PointwiseRecorder.h"

#include <sstream>
#include <iomanip>
#include "XMPI.h"

Receiver::Receiver(const std::string &name, const std::string &network, 
    double theta_lat, double phi_lon, bool geographic, 
    double depth, bool dumpStrain, double srcLat, double srcLon, double srcDep):
mName(name), mNetwork(network), mDepth(depth), mDumpStrain(dumpStrain) {
    RDCol3 rtpG, rtpS;
    if (geographic) {
        rtpG(0) = 1.;
        rtpG(1) = Geodesy::lat2Theta_d(theta_lat, mDepth);
        rtpG(2) = Geodesy::lon2Phi(phi_lon);
        rtpS = Geodesy::rotateGlob2Src(rtpG, srcLat, srcLon, srcDep);
    } else {
        rtpS(0) = 1.;
        rtpS(1) = theta_lat * degree;
        rtpS(2) = phi_lon * degree;
        rtpG = Geodesy::rotateSrc2Glob(rtpS, srcLat, srcLon, srcDep);
    }
    mTheta = rtpS(1);
    mPhi = rtpS(2);
    mLat = Geodesy::theta2Lat_d(rtpG(1), mDepth);
    mLon = Geodesy::phi2Lon(rtpG(2));
    mBackAzimuth = Geodesy::backAzimuth(srcLat, srcLon, srcDep, mLat, mLon, mDepth);
    // // test
    // XMPI::cout << name << " " << network << " "; 
    // XMPI::cout << mLat << " " << mLon << " " << " 0.0 " << mDepth << XMPI::endl;
}

void Receiver::release(PointwiseRecorder &recorderPW, const Domain &domain, 
    int elemTag, const RDMatPP &interpFact) {
    Element *myElem = domain.getElement(elemTag);
    if (mDumpStrain) {
        myElem->forceTIso();
    }
    recorderPW.addReceiver(mName, mNetwork, mPhi, interpFact, myElem, mTheta, mBackAzimuth,
        mLat, mLon, mDepth, mDumpStrain);
}

bool Receiver::locate(const Mesh &mesh, int &elemTag, RDMatPP &interpFact) const {
    RDCol2 recCrds, srcXiEta;
    double r = mesh.computeRadiusRef(mDepth, mLat, mLon);
    recCrds(0) = r * sin(mTheta);
    recCrds(1) = r * cos(mTheta);
    if (recCrds(0) > mesh.sMax() + tinySingle || recCrds(0) < mesh.sMin() - tinySingle) {
        return false;
    }
    if (recCrds(1) > mesh.zMax() + tinySingle || recCrds(1) < mesh.zMin() - tinySingle) {
        return false;
    }
    for (int iloc = 0; iloc < mesh.getNumQuads(); iloc++) {
        const Quad *quad = mesh.getQuad(iloc);
        if (!quad->nearMe(recCrds(0), recCrds(1))) {
            continue;
        }
        if (quad->invMapping(recCrds, srcXiEta)) {
            if (std::abs(srcXiEta(0)) <= 1.000001 && std::abs(srcXiEta(1)) <= 1.000001) {
                elemTag = quad->getElementTag();
                RDColP interpXi, interpEta;
                XMath::interpLagrange(srcXiEta(0), nPntEdge, 
                    quad->isAxial() ? SpectralConstants::getP_GLJ().data(): 
                    SpectralConstants::getP_GLL().data(), interpXi.data());
                XMath::interpLagrange(srcXiEta(1), nPntEdge, 
                    SpectralConstants::getP_GLL().data(), interpEta.data());
                interpFact = interpXi * interpEta.transpose();
                return true;
            }
        }
    }
    return false;
}

std::string Receiver::verbose(bool geographic, int wname, int wnet) const {
    std::stringstream ss;
    ss << std::setw(wname) << mName << "   ";
    ss << std::setw(wnet) << mNetwork << "   ";
    if (geographic) {
        ss << std::setw(8) << mLat << "   ";
        ss << std::setw(8) << mLon << "   ";
    } else {
        ss << std::setw(8) << mTheta / degree << "   ";
        ss << std::setw(8) << mPhi / degree << "   ";
    }
    ss << mDepth;
    return ss.str();
}

