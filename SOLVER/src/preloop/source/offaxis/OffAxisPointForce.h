// OffAxisPointForce.h
// created by Kuangdai on 11-Nov-2017
// off-axis point-force source

#pragma once
#include "OffAxisSource.h"

class OffAxisPointForce: public OffAxisSource {
public:
    
    OffAxisPointForce(double depth, double lat, double lon,
        double srcLat, double srcLon, double srcDep,
        RDMatX3 q_sphiz);

    std::string verbose() const;

protected:
    void computeSourceFourier(const Quad &myQuad, 
        const RDColP &interpFactXii,
        const RDColP &interpFactEta,
        double phi,
        vec_arPP_CMatX3 &fouriers) const;

private:
    // store: cylindrical components (q_s, q_phi, q_z)
    RDMatX3 mQ_sphiz;
};
