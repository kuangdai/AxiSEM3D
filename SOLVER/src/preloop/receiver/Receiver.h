// Receiver.h
// created by Kuangdai on 13-May-2016 
// receiver   

#pragma once

#include <string>
#include "eigenp.h"

class Domain;
class Mesh;

class Receiver {
public:
    Receiver(const std::string &name, const std::string &network, 
        double theta_lat, double phi_lon, bool geographic, 
        double depth, double srcLat, double srcLon, double srcDep);
    
    void release(Domain &domain, const Mesh &mesh, 
        int recordInterval, int component,
        const std::string &path, bool binary, bool append, int bufferSize); 
    
    bool locate(const Mesh &mesh, int &elemTag, RDMatPP &interpFact) const;
    
    std::string verbose(bool geographic, int wname, int wnet) const;
    
private:
    std::string mName;
    std::string mNetwork;
    double mTheta, mPhi;
    double mLat, mLon;
    double mDepth;
    double mBackAzimuth;
};

