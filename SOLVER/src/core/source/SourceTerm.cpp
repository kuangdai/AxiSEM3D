// SourceTerm.cpp
// created by Kuangdai on 5-Apr-2016 
// base class of sources in the solver
// NOTE: there are many different types of seismic sources, such as earthquakes,
// explosions, surface forces (ocean), and point forces (adjoint). But to the 
// solver, all those sources are seen as a combination of force vectors on solid 
// points, or of scalars on fluid points if any.
// So, if we only consider sources in the solid domain, only one class is enough
// to implement all source types.

#include "SourceTerm.h"
#include "Element.h"
#include "Point.h"

SourceTerm::SourceTerm(Element *element, const arPP_CMatX3 &force):
mElement(element) {
    // make the order consistent
    for (int i = 0; i < nPntElem; i++) {
        int length = mElement->getPoint(i)->getNu() + 1;
        if (length < force[i].rows()) {
            mForce[i] = force[i].topRows(length);
        } else {
            mForce[i] = force[i];
        }
    }
    mForceXSTF = mForce;
}

void SourceTerm::apply(Real stf) {
    // using mForceXSTF avoids dynamic allocation
    for (int i = 0; i < nPntElem; i++) {
        mForceXSTF[i] = mForce[i] * stf;
    }
    mElement->addSourceTerm(mForceXSTF);
}