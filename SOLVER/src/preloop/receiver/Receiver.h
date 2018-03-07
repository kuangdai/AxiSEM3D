// Receiver.h
// created by Kuangdai on 13-May-2016 
// receiver   

#pragma once

#include <string>
#include "eigenp.h"

class Domain;
class Mesh;
class PointwiseRecorder;

class Receiver {
public:
    Receiver(const std::string &name, const std::string &network, 
        double theta_lat, double phi_lon, bool geographic, 
        double depth, bool dumpStrain, double srcLat, double srcLon, double srcDep);
    
    void release(PointwiseRecorder &recorderPW, const Domain &domain, 
        int elemTag, const RDMatPP &interpFact);     
    
    bool locate(const Mesh &mesh, int &elemTag, RDMatPP &interpFact) const;
    
    std::string verbose(bool geographic, int wname, int wnet) const;
    
    const std::string &getName() const {return mName;};
    const std::string &getNetwork() const {return mNetwork;};
    
private:
    std::string mName;
    std::string mNetwork;
    double mTheta, mPhi;
    double mLat, mLon;
    double mDepth;
    double mBackAzimuth;
    bool mDumpStrain;
};

