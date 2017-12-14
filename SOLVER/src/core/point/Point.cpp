// Point.h
// created by Kuangdai on 20-Apr-2016 
// base class of gll points

#include "Point.h"

Point::Point(int nr, bool axial, const RDCol2 &crds):
mNr(nr), mNu(nr / 2), mAxial(axial), mCoords(crds) {
    // nothing
}

void Point::scatterDisplToElement(vec_ar3_CMatPP &displ, int ipol, int jpol, int maxNu) const {
    throw std::runtime_error("Point::scatterDisplToElement || Incompatible point type.");
}

void Point::scatterDisplToElement(vec_CMatPP &displ, int ipol, int jpol, int maxNu) const {
    throw std::runtime_error("Point::scatterDisplToElement || Incompatible point type.");
}

void Point::gatherStiffFromElement(const vec_ar3_CMatPP &stiff, int ipol, int jpol) {
    throw std::runtime_error("Point::gatherStiffFromElement || Incompatible point type.");
}

void Point::gatherStiffFromElement(const vec_CMatPP &stiff, int ipol, int jpol) {
    throw std::runtime_error("Point::gatherStiffFromElement || Incompatible point type.");
}

void Point::addToStiff(const CMatX3 &source) {
    throw std::runtime_error("Point::addToStiff || Incompatible point type.");
}

std::string Point::costSignature() const {
    std::stringstream ss;
    ss << verbose() << "$DimAzimuth=" << mNr << "$Axial=" << (axial() ? "T" : "F");
    return ss.str();
}

const CMatX3 &Point::getDispFourierSolid() const {
    throw std::runtime_error("Point::getDispFourier || Incompatible point type.");
}

const CColX &Point::getDispFourierFluid() const {
    throw std::runtime_error("Point::getDispFourier || Incompatible point type.");
}