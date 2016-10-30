// Element.cpp
// created by Kuangdai on 27-Mar-2016 
// base class of axisymmetric spectral elements

#include "Element.h"
#include "Gradient.h"
#include "Point.h"

Element::Element(Gradient *grad, const std::array<Point *, nPntElem> &points):
mGradient(grad) {
    mMaxNr = mMaxNu = -1;
    for (int i = 0; i < nPntElem; i++) {
        mPoints[i] = points[i];
        mMaxNr = std::max(mMaxNr, mPoints[i]->getNr());
        mMaxNu = std::max(mMaxNu, mPoints[i]->getNu());
    }
}

Element::~Element() {
    delete mGradient;
}

void Element::addSourceTerm(const arPP_CMatX3 &source) const {
    for (int i = 0; i < nPntElem; i++) 
        mPoints[i]->addToStiff(source[i]);
}

int Element::sizeComm() const {
    int ipol, jpol;
    
    int size0 = 0;
    ipol = 0;
    for (jpol = 0; jpol <= nPol; jpol++) size0 += mPoints[ipol * nPntEdge + jpol]->sizeComm();
    
    int size1 = 0;
    ipol = nPol;
    for (jpol = 0; jpol <= nPol; jpol++) size1 += mPoints[ipol * nPntEdge + jpol]->sizeComm();
        
    int size2 = 0;
    jpol = 0;
    for (ipol = 0; ipol <= nPol; ipol++) size2 += mPoints[ipol * nPntEdge + jpol]->sizeComm();
    
    int size3 = 0;
    jpol = nPol;
    for (ipol = 0; ipol <= nPol; ipol++) size3 += mPoints[ipol * nPntEdge + jpol]->sizeComm();
    
    return std::max(size0, std::max(size1, std::max(size2, size3)));
}

bool Element::axial() const {
    return mPoints[0]->axial();
}

std::string Element::costSignature() const {
    std::stringstream ss;
    ss << verbose() << "$DimAzimuth=" << mMaxNr << "$Axial=" << (axial() ? "T" : "F");
    return ss.str();
}

