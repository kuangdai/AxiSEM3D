// OffAxisPointForce.h
// created by Kuangdai on 11-Nov-2017
// off-axis point-force source

#include "OffAxisPointForce.h"
#include "Quad.h"
#include "SpectralConstants.h"
#include "XMath.h"
#include <sstream>

#include "Relabelling.h"

OffAxisPointForce::OffAxisPointForce(double depth, double lat, double lon,
    double srcLat, double srcLon, double srcDep,
    RDMatX3 q_sphiz): OffAxisSource(depth, lat, lon, srcLat, srcLon, srcDep),
    mQ_sphiz(q_sphiz) {
    // nothing
}

void OffAxisPointForce::computeSourceFourier(const Quad &myQuad, 
    const RDColP &interpFactXii,
    const RDColP &interpFactEta,
    double phi, 
    vec_arPP_CMatX3 &fouriers) const {
    // number of steps
    int nstep = mQ_sphiz.rows();
    // Fourier order
    int nu = myQuad.getNu();
    // set zero
    for (int istep = 0; istep < nstep; istep++) {
        for (int ipnt = 0; ipnt < nPntElem; ipnt++) {
            fouriers[istep][ipnt] = CMatX3::Zero(nu + 1, 3);
        }
    }
    // particle relabelling
    RDRowN JPRT;
    if (myQuad.hasRelabelling()) {
        const RDMatXN &JJ = myQuad.getRelabelling().getStiffJacobian();
        JPRT = XMath::computeFourierAtPhi(JJ, phi);
    } else {
        JPRT = RDRowN::Ones();
    }
    // compute source pointwise
    for (int istep = 0; istep < nstep; istep++) {  // stf
        for (int ipol = 0; ipol <= nPol; ipol++) {
            for (int jpol = 0; jpol <= nPol; jpol++) {
                int ipnt = ipol * nPntEdge + jpol;
                double fact = interpFactXii(ipol) * interpFactEta(jpol);
                for (int beta = 0; beta <= nu; beta++) {
                    for (int idim = 0; idim < 3; idim++) {
                        fouriers[istep][ipnt](beta, idim) = Complex(
                            fact * exp(beta * phi * iid) * mQ_sphiz(istep, idim) * JPRT(ipnt));
                    }
                }
            }
        }
    }
}

std::string OffAxisPointForce::verbose() const {
    // TODO: should be verbosed collectively like receivers
    return "";
}
