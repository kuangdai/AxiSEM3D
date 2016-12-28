// Receiver.cpp
// created by Kuangdai on 13-May-2016 
// receiver   

#include "Receiver.h"
#include "XMath.h"
#include "SeismometerRTZ.h"
#include "SeismometerENZ.h"
#include "RecorderAscii.h"
#include "RecorderBinary.h"
#include "Station.h"
#include "Domain.h"
#include "Mesh.h"
#include "Quad.h"
#include "Element.h"
#include "SpectralConstants.h"

#include <sstream>
#include <iomanip>
#include "XMPI.h"

Receiver::Receiver(const std::string &name, const std::string &network, 
    double theta_lat, double phi_lon, bool geographic, 
    double depth, double srcLat, double srcLon, double srcDep):
mName(name), mNetwork(network), mDepth(depth) {
    RDCol3 rtpG, rtpS;
    if (geographic) {
        rtpG(0) = 1.;
        rtpG(1) = XMath::lat2Theta(theta_lat, mDepth);
        rtpG(2) = XMath::lon2Phi(phi_lon);
        rtpS = XMath::rotateGlob2Src(rtpG, srcLat, srcLon, srcDep);
    } else {
        rtpS(0) = 1.;
        rtpS(1) = theta_lat * degree;
        rtpS(2) = phi_lon * degree;
        rtpG = XMath::rotateSrc2Glob(rtpS, srcLat, srcLon, srcDep);
    }
    mTheta = rtpS(1);
    mPhi = rtpS(2);
    mLat = XMath::theta2Lat(rtpG(1), mDepth);
    mLon = XMath::phi2Lon(rtpG(2));
    mBackAzimuth = XMath::backAzimuth(srcLat, srcLon, srcDep, mLat, mLon, mDepth);
    // // test
    // XMPI::cout << "6371_" << round(mTheta / degree) << "_" << round(mPhi / degree) << " TEST "; 
    // XMPI::cout << mLat << " " << mLon << " " << " 0.0 " << mDepth << XMPI::endl;
}

void Receiver::release(Domain &domain, const Mesh &mesh, 
    int recordInterval, int component,
    const std::string &path, bool binary, bool append, int bufferSize,
    int elemTag, const RDMatPP &interpFact) {
                    
    Element *myElem = domain.getElement(elemTag);
                    
    // seismometer
    Seismometer *seis;
    if (component == 0) 
        seis = new SeismometerRTZ(mPhi, interpFact.cast<Real>(), myElem, mTheta);
    else if (component == 1)
        seis = new SeismometerENZ(mPhi, interpFact.cast<Real>(), myElem, mTheta, mBackAzimuth * degree);
    else 
        seis = new Seismometer(mPhi, interpFact.cast<Real>(), myElem);
                
    // recorder
    Recorder *rec;
    std::string fname = path + "/" + mNetwork + "_" + mName;
    if (binary) 
        rec = new RecorderBinary(bufferSize, fname, append);
    else 
        rec = new RecorderAscii(bufferSize, fname, append);
    
    // station    
    domain.addStation(new Station(recordInterval, seis, rec));
}

bool Receiver::locate(const Mesh &mesh, int &elemTag, RDMatPP &interpFact) const {
    RDCol2 recCrds, srcXiEta;
    double r = mesh.computeRadiusRef(mDepth, mLat, mLon);
    recCrds(0) = r * sin(mTheta);
    recCrds(1) = r * cos(mTheta);
    if (recCrds(0) > mesh.sMax() + tinySingle || recCrds(0) < mesh.sMin() - tinySingle) return false;
    if (recCrds(1) > mesh.zMax() + tinySingle || recCrds(1) < mesh.zMin() - tinySingle) return false;
    for (int iloc = 0; iloc < mesh.getNumQuads(); iloc++) {
        const Quad *quad = mesh.getQuad(iloc);
        if (!quad->nearMe(recCrds(0), recCrds(1))) continue;
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

