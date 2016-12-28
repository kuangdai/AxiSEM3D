// poinf.h
// created by Alex on 25-oct-2016
// earthquake source

#pragma once
#include "Source.h"

class pointf: public Source {
public:
    
    pointf(double depth = 0., double lat = 0., double lon = 0.,
        double f1 = 0., double f2 = 0., double f3 = 0.);

    std::string verbose() const;

protected:
    void computeSourceFourier(const Quad &myQuad, const RDColP &interpFactZ,
        arPP_CMatX3 &fouriers) const;

private:
    // store: Cartesian, paper components
    double px;
    double py;
    double pz;

};
