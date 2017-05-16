// SpectralConstants.h
// created by Kuangdai on 3-May-2016 
// create spectral constants

#pragma once

#include "eigenp.h"

class SpectralConstants {
public:
    static void initialize(int nPol);
    static RDCol2 getXiEta(int xiPol, int etaPol, bool axial);    
    static RDCol2 getWeights(int xiPol, int etaPol, bool axial);
    static const RDColP &getP_GLL() {return sP_GLL;};
    static const RDColP &getP_GLJ() {return sP_GLJ;};
    static const RDMatPP &getG_GLL() {return sG_GLL;};
    static const RDMatPP &getG_GLJ() {return sG_GLJ;};
    
private:
    static RDColP sP_GLL;
    static RDColP sP_GLJ;
    static RDColP sW_GLL;
    static RDColP sW_GLJ;
    static RDMatPP sG_GLL;
    static RDMatPP sG_GLJ;
};
