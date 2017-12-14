// Source.h
// created by Kuangdai on 8-May-2016 
// base class of source
// we only consider point sources located on the axis

#pragma once

#include "eigenp.h"
#include "eigenc.h"

class Quad;
class Mesh;
class Domain;
class Parameters;

class Source {
public:
    Source(double depth = 0., double lat = 90., double lon = 0.);
    
    virtual ~Source() {};
    
    void release(Domain &domain, const Mesh &mesh) const;
    
    virtual std::string verbose() const = 0;
        
    double getLatitude() const {return mLatitude;};
    double getLongitude() const {return mLongitude;};
    double getDepth() const {return mDepth;};
    
    static void parseLine(const std::string &line, const std::string &key, double &res);
    static void checkValue(const std::string &key, double res);
    static void buildInparam(Source *&src, const Parameters &par, int verbose);
    
protected:
    virtual void computeSourceFourier(const Quad &myQuad, const RDColP &interpFactZ,
        arPP_CMatX3 &fouriers) const = 0;
        
    double mDepth;
    double mLatitude;
    double mLongitude;
        
private:
    bool locate(const Mesh &mesh, int &locTag, RDColP &interpFactZ) const;
};

