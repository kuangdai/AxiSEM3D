// SurfaceInfo.cpp
// created by Kuangdai on 27-Nov-2017 
// surface element info

#include "SurfaceInfo.h"
#include "Element.h"
#include "Geodesy.h"

SurfaceInfo::SurfaceInfo(const Element *ele, int surfSide):
    mElement(ele), mSurfSide(surfSide) {
    const RDMatXX &sz = mElement->getCoordsOnSide(surfSide);
    const RDCol2 &sz0 = sz.col(0);
    const RDCol2 &sz1 = sz.col(nPol);
    mTheta0 = Geodesy::theta(sz0);
    mTheta1 = Geodesy::theta(sz1);
}
    
void SurfaceInfo::initBuffer(int bufferSize, CMatXX_RM &bufferDisp) {
    bufferDisp = CMatXX_RM::Zero(bufferSize, 
        nPntEdge * 3 * (mElement->getMaxNu() + 1));
}

void SurfaceInfo::feedBuffer(int bufferLine, CMatXX_RM &bufferDisp) {
    mElement->feedDispOnSide(mSurfSide, bufferDisp, bufferLine);
}

int SurfaceInfo::getMaxNu() const {
    return mElement->getMaxNu();
}
