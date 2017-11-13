// Earthquake.h
// created by Kuangdai on 8-May-2016 
// axial earthquake source

#pragma once
#include "Source.h"

class Earthquake: public Source {
public:
    // input: CMTSOLUTION components
    Earthquake(double depth = 0., double lat = 0., double lon = 0., 
        double Mrr = 0., double Mtt = 0., double Mpp = 0., 
        double Mrt = 0., double Mrp = 0., double Mtp = 0.);
        
    std::string verbose() const;
    
protected:    
    void computeSourceFourier(const Quad &myQuad, const RDColP &interpFactZ,
        arPP_CMatX3 &fouriers) const;
        
private:
    // store: Cartesian, paper components
    double mMxx;
    double mMyy;
    double mMzz;
    double mMxy;
    double mMxz;
    double mMyz;
};

