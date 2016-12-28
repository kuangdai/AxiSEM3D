// NullSource.cpp
// created by Kuangdai on 1-Nov-2016 
// null source, for scaling test only

#include "NullSource.h"
#include <sstream>

std::string NullSource::verbose() const {
    std::stringstream ss;
    ss << "\n========================== Source ==========================" << std::endl;
    ss << "  I am a null source for test purposes." << std::endl;
    ss << "========================== Source ==========================\n" << std::endl;
    return ss.str();
}

void NullSource::computeSourceFourier(const Quad &myQuad, const RDColP &interpFactZ,
    arPP_CMatX3 &fouriers) const {
    // just initialize as zero
    for (int ipnt = 0; ipnt < nPntElem; ipnt++) {
        fouriers[ipnt] = CMatX3::Zero(0, 3);
    }
}