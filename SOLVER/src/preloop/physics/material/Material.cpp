// Material.cpp
// created by Kuangdai on 17-May-2016 
// 3D seismic material properties

#include "Material.h"
#include "Quad.h"
#include "ExodusModel.h"
#include "SpectralConstants.h"
#include "XMath.h"

#include "Volumetric3D.h"
#include "Relabelling.h"
#include "Acoustic1D.h"
#include "Acoustic3X.h"

#include "AttBuilder.h"
#include "Attenuation1D_CG4.h"
#include "Attenuation1D_Full.h"
#include "Attenuation3D_CG4.h"
#include "Attenuation3D_Full.h"

#include "Isotropic1D.h"
#include "Isotropic3D.h"
#include "Isotropic3X.h"
#include "TransverselyIsotropic1D.h"
#include "TransverselyIsotropic3D.h"
#include "TransverselyIsotropic3X.h"

Material::Material(const Quad *myQuad, const ExodusModel &exModel): mMyQuad(myQuad) {
    // read Exodus model
    int quadTag = mMyQuad->getQuadTag();
    if (exModel.isIsotropic()) {
        for (int i = 0; i < 4; i++) {
            mVpv1D(i) = mVph1D(i) = exModel.getElementalVariables().at("VP_" + std::to_string(i))[quadTag];
            mVsv1D(i) = mVsh1D(i) = exModel.getElementalVariables().at("VS_" + std::to_string(i))[quadTag];
            mEta(i) = 1.;
        }
    } else {
        for (int i = 0; i < 4; i++) {
            mVpv1D(i) = exModel.getElementalVariables().at("VPV_" + std::to_string(i))[quadTag];
            mVph1D(i) = exModel.getElementalVariables().at("VPH_" + std::to_string(i))[quadTag];
            mVsv1D(i) = exModel.getElementalVariables().at("VSV_" + std::to_string(i))[quadTag];
            mVsh1D(i) = exModel.getElementalVariables().at("VSH_" + std::to_string(i))[quadTag];
            mEta(i) = exModel.getElementalVariables().at("ETA_" + std::to_string(i))[quadTag];
        }
    }
    for (int i = 0; i < 4; i++) {
        mRho1D(i) = exModel.getElementalVariables().at("RHO_" + std::to_string(i))[quadTag];
        mQmu(i) = exModel.getElementalVariables().at("QMU_" + std::to_string(i))[quadTag];
        mQkp(i) = exModel.getElementalVariables().at("QKAPPA_" + std::to_string(i))[quadTag];
    }
    
    // initialize 3D properties with 1D reference
    // RDMatXN mVpv3D, mVph3D;
    // RDMatXN mVsv3D, mVsh3D;
    // RDMatXN mRho3D;
    // arPP_RDColX mRhoMass3D;
    int Nr = mMyQuad->getNr();
    mVpv3D = RDMatXN::Zero(Nr, nPE);
    mVph3D = RDMatXN::Zero(Nr, nPE);
    mVsv3D = RDMatXN::Zero(Nr, nPE);
    mVsh3D = RDMatXN::Zero(Nr, nPE);
    mRho3D = RDMatXN::Zero(Nr, nPE);
    for (int ipol = 0; ipol <= nPol; ipol++) {
        for (int jpol = 0; jpol <= nPol; jpol++) {
            int ipnt = ipol * nPntEdge + jpol;
            const RDCol2 &xieta = SpectralConstants::getXiEta(ipol, jpol, mMyQuad->isAxial());
            // fill with 1D 
            mVpv3D.col(ipnt).fill(Mapping::interpolate(mVpv1D, xieta));
            mVph3D.col(ipnt).fill(Mapping::interpolate(mVph1D, xieta));
            mVsv3D.col(ipnt).fill(Mapping::interpolate(mVsv1D, xieta));
            mVsh3D.col(ipnt).fill(Mapping::interpolate(mVsh1D, xieta));
            mRho3D.col(ipnt).fill(Mapping::interpolate(mRho1D, xieta));
            // rho for mass
            int NrP = mMyQuad->getPointNr(ipol, jpol);
            mRhoMass3D[ipnt] = RDColX::Constant(NrP, mRho3D.col(ipnt)[0]);
        }
    }
}

void Material::addVolumetric3D(const Volumetric3D &m3D, double srcLat, double srcLon, double srcDep, double phi2D) {
    // radius at element center 
    double rElemCenter = mMyQuad->computeCenterRadius();
    
    // read 3D model
    int Nr = mMyQuad->getNr();
    for (int ipol = 0; ipol <= nPol; ipol++) {
        for (int jpol = 0; jpol <= nPol; jpol++) {
            const RDCol2 &xieta = SpectralConstants::getXiEta(ipol, jpol, mMyQuad->isAxial());
            // geographic oordinates of cardinal points
            const RDMatX3 &rtp = mMyQuad->computeGeocentricGlobal(srcLat, srcLon, srcDep, xieta, Nr, phi2D);
            // 1D reference values
            double vpv_ref = Mapping::interpolate(mVpv1D, xieta);
            double vph_ref = Mapping::interpolate(mVph1D, xieta);
            double vsv_ref = Mapping::interpolate(mVsv1D, xieta);
            double vsh_ref = Mapping::interpolate(mVsh1D, xieta);
            double rho_ref = Mapping::interpolate(mRho1D, xieta);
            int ipnt = ipol * nPntEdge + jpol;
            for (int alpha = 0; alpha < Nr; alpha++) {
                double dvpv, dvph, dvsv, dvsh, drho;
                if (m3D.get3dProperties(rtp(alpha, 0), rtp(alpha, 1), rtp(alpha, 2), rElemCenter,
                        dvpv, dvph, dvsv, dvsh, drho)) {
                    if (m3D.getReferenceType() == Volumetric3D::ReferenceTypes::Absolute) {
                        mVpv3D(alpha, ipnt) = dvpv;
                        mVph3D(alpha, ipnt) = dvph;
                        mVsv3D(alpha, ipnt) = dvsv;
                        mVsh3D(alpha, ipnt) = dvsh;
                        mRho3D(alpha, ipnt) = drho;
                    } else if (m3D.getReferenceType() == Volumetric3D::ReferenceTypes::Reference1D) {
                        mVpv3D(alpha, ipnt) = vpv_ref * (1. + dvpv);
                        mVph3D(alpha, ipnt) = vph_ref * (1. + dvph);
                        mVsv3D(alpha, ipnt) = vsv_ref * (1. + dvsv);
                        mVsh3D(alpha, ipnt) = vsh_ref * (1. + dvsh);
                        mRho3D(alpha, ipnt) = rho_ref * (1. + drho);
                    } else if (m3D.getReferenceType() == Volumetric3D::ReferenceTypes::Reference3D) {
                        mVpv3D(alpha, ipnt) *= 1. + dvpv;
                        mVph3D(alpha, ipnt) *= 1. + dvph;
                        mVsv3D(alpha, ipnt) *= 1. + dvsv;
                        mVsh3D(alpha, ipnt) *= 1. + dvsh;
                        mRho3D(alpha, ipnt) *= 1. + drho;
                    } else {
                        mVpv3D(alpha, ipnt) = (mVpv3D(alpha, ipnt) - vpv_ref) * (1. + dvpv) + vpv_ref;
                        mVph3D(alpha, ipnt) = (mVph3D(alpha, ipnt) - vph_ref) * (1. + dvph) + vph_ref;
                        mVsv3D(alpha, ipnt) = (mVsv3D(alpha, ipnt) - vsv_ref) * (1. + dvsv) + vsv_ref;
                        mVsh3D(alpha, ipnt) = (mVsh3D(alpha, ipnt) - vsh_ref) * (1. + dvsh) + vsh_ref;
                        mRho3D(alpha, ipnt) = (mRho3D(alpha, ipnt) - rho_ref) * (1. + drho) + rho_ref;
                    }
                } 
            }
            // rho for mass
            int NrP = mMyQuad->getPointNr(ipol, jpol);
            const RDMatX3 &rtpM = mMyQuad->computeGeocentricGlobal(srcLat, srcLon, srcDep, xieta, NrP, phi2D);
            for (int alpha = 0; alpha < NrP; alpha++) {
                double dvpv, dvph, dvsv, dvsh, drho;
                if (m3D.get3dProperties(rtpM(alpha, 0), rtpM(alpha, 1), rtpM(alpha, 2), rElemCenter,
                        dvpv, dvph, dvsv, dvsh, drho)) {
                    if (m3D.getReferenceType() == Volumetric3D::ReferenceTypes::Absolute) {
                        mRhoMass3D[ipnt](alpha) = drho;
                    } else if (m3D.getReferenceType() == Volumetric3D::ReferenceTypes::Reference1D) {
                        mRhoMass3D[ipnt](alpha) = rho_ref * (1. + drho);
                    } else if (m3D.getReferenceType() == Volumetric3D::ReferenceTypes::Reference3D) {
                        mRhoMass3D[ipnt](alpha) *= 1. + drho;
                    } else {
                        mRhoMass3D[ipnt](alpha) = (mRhoMass3D[ipnt](alpha) - vpv_ref) * (1. + drho) + rho_ref;
                    }
                }
            }
        }
    }
}

arPP_RDColX Material::computeElementalMass() const {
    arPP_RDColX mass, J;
    // Jacobian of topography
    if (mMyQuad->massRelabelling()) {
        J = mMyQuad->getRelabelling().getMassJacobian();
    } else {
        J = mRhoMass3D; // just to use the size
        for (int ipnt = 0; ipnt < nPE; ipnt++) J[ipnt].setOnes();
    }
    // general mass term
    const RDMatPP &iFact = mMyQuad->getIntegralFactor();
    for (int ipol = 0; ipol <= nPol; ipol++) {
        for (int jpol = 0; jpol <= nPol; jpol++) {
            int ipnt = ipol * nPntEdge + jpol;
            if (mMyQuad->isFluid()) {
                const RDCol2 &xieta = SpectralConstants::getXiEta(ipol, jpol, mMyQuad->isAxial());
                double vp = Mapping::interpolate(mVpv1D, xieta);
                double rho = Mapping::interpolate(mRho1D, xieta);
                mass[ipnt] = J[ipnt] * (1. / (rho * vp * vp) * iFact(ipol, jpol));
            } else {
                mass[ipnt] = J[ipnt].schur(mRhoMass3D[ipnt]) * iFact(ipol, jpol);
            }
        }
    }
    return mass;
}

Acoustic *Material::createAcoustic() const {
    const RDMatPP &iFact = mMyQuad->getIntegralFactor(); 
    if (mMyQuad->stiffRelabelling()) {
        RDRowN fluidK, theta;
        for (int ipol = 0; ipol <= nPol; ipol++) {
            for (int jpol = 0; jpol <= nPol; jpol++) {
                int ipnt = ipol * nPntEdge + jpol;
                const RDCol2 &xieta = SpectralConstants::getXiEta(ipol, jpol, mMyQuad->isAxial());
                double rho = Mapping::interpolate(mRho1D, xieta);
                fluidK(ipnt) = 1. / rho * iFact(ipol, jpol);
                theta(ipnt) = XMath::theta(mMyQuad->mapping(xieta));
            }
        }
        const RDMatXN &J = mMyQuad->getRelabelling().getStiffJacobian();
        const RDMatXN &KJ = fluidK.colwise().replicate(J.rows()).schur(J);
        return new Acoustic3X(theta, XMath::castToSolver<RDMatXN, RMatXN>(KJ),
            XMath::castToSolver<RDMatXN4, RMatXN4>(mMyQuad->getRelabelling().getStiffX()));
    } else {
        RDMatPP fluidK;
        for (int ipol = 0; ipol <= nPol; ipol++) {
            for (int jpol = 0; jpol <= nPol; jpol++) {
                const RDCol2 &xieta = SpectralConstants::getXiEta(ipol, jpol, mMyQuad->isAxial());
                double rho = Mapping::interpolate(mRho1D, xieta);
                fluidK(ipol, jpol) = 1. / rho * iFact(ipol, jpol);
            }
        }
        return new Acoustic1D(XMath::castToSolver(fluidK));
    }
}

Elastic *Material::createElastic(const AttBuilder *attBuild) const {
    if (isStiffness1D() && !mMyQuad->stiffRelabelling()) {
        return createElastic1D(attBuild);
    } else {
        return createElastic3D(attBuild);
    }
}

double Material::getVMaxRef() const {
    return std::max(mVph1D.maxCoeff(), mVpv1D.maxCoeff());
}

RDColX Material::getVMax() const {
    const RDColX &vpvMax = mVpv3D.rowwise().maxCoeff();
    const RDColX &vphMax = mVph3D.rowwise().maxCoeff();
    return (vpvMax.array().max(vphMax.array())).matrix();
}

bool Material::isStiffness1D() const {
    return XMath::equalRows(mVpv3D) && XMath::equalRows(mVph3D) &&
           XMath::equalRows(mVsv3D) && XMath::equalRows(mVsh3D) && 
           XMath::equalRows(mRho3D);
}

bool Material::isIsotropic() const {
    return (mVpv3D - mVph3D).norm() < tinyDouble * mVpv3D.norm() &&  
           (mVsv3D - mVsh3D).norm() < tinyDouble * mVsv3D.norm() && 
           (mEta - RDRow4::Ones()).norm() < tinyDouble;
}

Elastic *Material::createElastic1D(const AttBuilder *attBuild) const {
    const RDMatPP &iFact = mMyQuad->getIntegralFactor(); 
    RDMatPP theta, A, C, F, L, N;
    for (int ipol = 0; ipol <= nPol; ipol++) {
        for (int jpol = 0; jpol <= nPol; jpol++) {
            int ipnt = ipol * nPntEdge + jpol;
            const RDCol2 &xieta = SpectralConstants::getXiEta(ipol, jpol, mMyQuad->isAxial());
            double vpv = mVpv3D(0, ipnt);
            double vph = mVph3D(0, ipnt);
            double vsv = mVsv3D(0, ipnt);
            double vsh = mVsh3D(0, ipnt);
            double rho = mRho3D(0, ipnt);
            double eta = Mapping::interpolate(mEta, xieta);
            A(ipol, jpol) = rho * vph * vph * iFact(ipol, jpol);
            C(ipol, jpol) = rho * vpv * vpv * iFact(ipol, jpol);
            L(ipol, jpol) = rho * vsv * vsv * iFact(ipol, jpol);
            N(ipol, jpol) = rho * vsh * vsh * iFact(ipol, jpol);
            F(ipol, jpol) = eta * (A(ipol, jpol) - 2. * L(ipol, jpol));
            theta(ipol, jpol) = XMath::theta(mMyQuad->mapping(xieta));
        }
    }
    
    Attenuation1D *att = 0;
    if (attBuild) {
        // Voigt average
        RDMatPP kappa = (4. * A + C + 4. * F - 4. * N) / 9.;
        RDMatPP mu = (A + C - 2. * F + 6. * L + 5. * N) / 15.;
        A -= (kappa + 4. / 3. * mu);
        C -= (kappa + 4. / 3. * mu);
        F -= (kappa - 2. / 3. * mu);
        L -= mu;
        N -= mu;
        makeAttenuation1D(*attBuild, kappa, mu, att);
        A += (kappa + 4. / 3. * mu);
        C += (kappa + 4. / 3. * mu);
        F += (kappa - 2. / 3. * mu);
        L += mu;
        N += mu;
    } 
    if (isIsotropic()) 
        return new Isotropic1D(XMath::castToSolver(F), XMath::castToSolver(L), att);
    else 
        return new TransverselyIsotropic1D(theta, 
            XMath::castToSolver(A), XMath::castToSolver(C), XMath::castToSolver(F), 
            XMath::castToSolver(L), XMath::castToSolver(N), att);
}

Elastic *Material::createElastic3D(const AttBuilder *attBuild) const {
    const RDMatPP &iFact = mMyQuad->getIntegralFactor(); 
    RDMatXN A(mRho3D), C(mRho3D), F(mRho3D), L(mRho3D), N(mRho3D);
    RDRowN theta;
    for (int ipol = 0; ipol <= nPol; ipol++) {
        for (int jpol = 0; jpol <= nPol; jpol++) {
            // A C L N
            int ipnt = ipol * nPntEdge + jpol;
            A.col(ipnt).array() *= mVph3D.col(ipnt).array().pow(2.) * iFact(ipol, jpol);
            C.col(ipnt).array() *= mVpv3D.col(ipnt).array().pow(2.) * iFact(ipol, jpol);
            L.col(ipnt).array() *= mVsv3D.col(ipnt).array().pow(2.) * iFact(ipol, jpol);
            N.col(ipnt).array() *= mVsh3D.col(ipnt).array().pow(2.) * iFact(ipol, jpol);
            // F
            const RDCol2 &xieta = SpectralConstants::getXiEta(ipol, jpol, mMyQuad->isAxial());
            double eta = Mapping::interpolate(mEta, xieta);
            F.col(ipnt) = eta * (A.col(ipnt) - 2. * L.col(ipnt));
            // theta
            theta(ipnt) = XMath::theta(mMyQuad->mapping(xieta));
        }
    }
    // must do relabelling before attenuation
    if (mMyQuad->stiffRelabelling()) {
        const RDMatXN &J = mMyQuad->getRelabelling().getStiffJacobian();
        A = A.schur(J);
        C = C.schur(J);
        F = F.schur(J);
        L = L.schur(J);
        N = N.schur(J);
    }
    
    Attenuation3D *att = 0;
    if (attBuild) {
        // Voigt average
        RDMatXN kappa = (4. * A + C + 4. * F - 4. * N) / 9.;
        RDMatXN mu = (A + C - 2. * F + 6. * L + 5. * N) / 15.;
        A -= (kappa + 4. / 3. * mu);
        C -= (kappa + 4. / 3. * mu);
        F -= (kappa - 2. / 3. * mu);
        L -= mu;
        N -= mu;
        makeAttenuation3D(*attBuild, kappa, mu, att);
        A += (kappa + 4. / 3. * mu);
        C += (kappa + 4. / 3. * mu);
        F += (kappa - 2. / 3. * mu);
        L += mu;
        N += mu;
    }
    
    if (mMyQuad->stiffRelabelling()) {
        if (isIsotropic()) 
            return new Isotropic3X(theta, 
                XMath::castToSolver<RDMatXN, RMatXN>(F), 
                XMath::castToSolver<RDMatXN, RMatXN>(L), 
                XMath::castToSolver<RDMatXN4, RMatXN4>(mMyQuad->getRelabelling().getStiffX()), att);
        else 
            return new TransverselyIsotropic3X(theta, 
                XMath::castToSolver<RDMatXN, RMatXN>(A),
                XMath::castToSolver<RDMatXN, RMatXN>(C),
                XMath::castToSolver<RDMatXN, RMatXN>(F),
                XMath::castToSolver<RDMatXN, RMatXN>(L),
                XMath::castToSolver<RDMatXN, RMatXN>(N), 
                XMath::castToSolver<RDMatXN4, RMatXN4>(mMyQuad->getRelabelling().getStiffX()), att);
    } else {
        if (isIsotropic()) 
            return new Isotropic3D(
                XMath::castToSolver<RDMatXN, RMatXN>(F), 
                XMath::castToSolver<RDMatXN, RMatXN>(L), att);
        else 
            return new TransverselyIsotropic3D(theta, 
                XMath::castToSolver<RDMatXN, RMatXN>(A),
                XMath::castToSolver<RDMatXN, RMatXN>(C),
                XMath::castToSolver<RDMatXN, RMatXN>(F),
                XMath::castToSolver<RDMatXN, RMatXN>(L),
                XMath::castToSolver<RDMatXN, RMatXN>(N), att);
    }
}

double Material::getFieldVariable(const std::string &vname, int ipol, int jpol, int islice, int refType) {
    int ipnt = ipol * nPntEdge + jpol;
    const RDCol2 &xieta = SpectralConstants::getXiEta(ipol, jpol, mMyQuad->isAxial());
    double vpv_ref = Mapping::interpolate(mVpv1D, xieta);
    double vph_ref = Mapping::interpolate(mVph1D, xieta);
    double vsv_ref = Mapping::interpolate(mVsv1D, xieta);
    double vsh_ref = Mapping::interpolate(mVsh1D, xieta);
    double rho_ref = Mapping::interpolate(mRho1D, xieta);
    
    if (refType == Volumetric3D::ReferenceTypes::Absolute || refType == Volumetric3D::ReferenceTypes::Reference3D) {
        if (vname == "vpv" || vname == "vp") return mVpv3D(islice, ipnt);
        if (vname == "vsv" || vname == "vs") return mVsv3D(islice, ipnt);
        if (vname == "vph") return mVph3D(islice, ipnt);
        if (vname == "vsh") return mVsh3D(islice, ipnt);
        if (vname == "rho") return mRho3D(islice, ipnt);
    } else if (refType == Volumetric3D::ReferenceTypes::Reference1D) {
        if (vname == "vpv" || vname == "vp") return vpv_ref;
        if (vname == "vsv" || vname == "vs") return vsv_ref;
        if (vname == "vph") return vph_ref;
        if (vname == "vsh") return vsh_ref;
        if (vname == "rho") return rho_ref;
    } else {
        if (vname == "vpv" || vname == "vp") return (mVpv3D(islice, ipnt) - vpv_ref) / vpv_ref;
        if (vname == "vph") return (mVph3D(islice, ipnt) - vph_ref) / vph_ref;
        if (vname == "rho") return (mRho3D(islice, ipnt) - rho_ref) / rho_ref;
        if (mMyQuad->isFluid()) {
            if (vname == "vsv" || vname == "vs" || vname == "vsh") return 0.;
        } else {
            if (vname == "vsv" || vname == "vs") return (mVsv3D(islice, ipnt) - vsv_ref) / vsv_ref;
            if (vname == "vsh") return (mVsh3D(islice, ipnt) - vsh_ref) / vsh_ref;
        }
    }
    throw std::runtime_error("Material::getFieldVariable || Unknown field variable name: " + vname);
}

void Material::makeAttenuation1D(const AttBuilder &attBuild, 
    RDMatPP &kappa, RDMatPP &mu, Attenuation1D *&att) const {
    // get attenuation factors
    RDColX alpha, beta, gamma;
    double dKappaFact, dMuFact, kappaFactAtt, muFactAtt, kappaFactNoAtt, muFactNoAtt;
    bool doKappa;
    double Qmu = mQmu.sum() / 4.;
    double Qkp = mQkp.sum() / 4.;
    attBuild.computeFactors(Qmu, Qkp, alpha, beta, gamma, dKappaFact, dMuFact, 
        kappaFactAtt, muFactAtt, kappaFactNoAtt, muFactNoAtt, doKappa);
    // create attenuation
    if (attBuild.useCG4()) {
        const RDRow4 &weights_cg4 = mMyQuad->computeWeightsCG4();
        RDRow4 dkappa, dmu;
        dkappa(0) = weights_cg4(0) * dKappaFact * kappa(1, 1);
        dkappa(1) = weights_cg4(1) * dKappaFact * kappa(1, 3);
        dkappa(2) = weights_cg4(2) * dKappaFact * kappa(3, 1);
        dkappa(3) = weights_cg4(3) * dKappaFact * kappa(3, 3);
        dmu(0) = weights_cg4(0) * dMuFact * mu(1, 1);
        dmu(1) = weights_cg4(1) * dMuFact * mu(1, 3);
        dmu(2) = weights_cg4(2) * dMuFact * mu(3, 1);
        dmu(3) = weights_cg4(3) * dMuFact * mu(3, 3);
        kappa *= kappaFactNoAtt;
        mu *= muFactNoAtt;
        kappa(1, 1) *= 1. + weights_cg4(0) * (kappaFactAtt / kappaFactNoAtt - 1.);
        kappa(1, 3) *= 1. + weights_cg4(1) * (kappaFactAtt / kappaFactNoAtt - 1.);
        kappa(3, 1) *= 1. + weights_cg4(2) * (kappaFactAtt / kappaFactNoAtt - 1.);
        kappa(3, 3) *= 1. + weights_cg4(3) * (kappaFactAtt / kappaFactNoAtt - 1.);
        mu(1, 1) *= 1. + weights_cg4(0) * (muFactAtt / muFactNoAtt - 1.);
        mu(1, 3) *= 1. + weights_cg4(1) * (muFactAtt / muFactNoAtt - 1.);
        mu(3, 1) *= 1. + weights_cg4(2) * (muFactAtt / muFactNoAtt - 1.);
        mu(3, 3) *= 1. + weights_cg4(3) * (muFactAtt / muFactNoAtt - 1.);
        att = new Attenuation1D_CG4(attBuild.getNSLS(), 
            alpha.cast<Real>(), beta.cast<Real>(), gamma.cast<Real>(), mMyQuad->getNu(), 
            dkappa.cast<Real>(), dmu.cast<Real>(), doKappa);
    } else {
        const RDMatPP &dkappa = dKappaFact * kappa;
        const RDMatPP &dmu = dMuFact * mu;
        kappa *= kappaFactAtt;
        mu *= muFactAtt;
        att = new Attenuation1D_Full(attBuild.getNSLS(), 
            alpha.cast<Real>(), beta.cast<Real>(), gamma.cast<Real>(), mMyQuad->getNu(), 
            XMath::castToSolver(dkappa), XMath::castToSolver(dmu), doKappa);
    }    
}

void Material::makeAttenuation3D(const AttBuilder &attBuild, 
    RDMatXN &kappa, RDMatXN &mu, Attenuation3D *&att) const {
    // get attenuation factors
    RDColX alpha, beta, gamma;
    double dKappaFact, dMuFact, kappaFactAtt, muFactAtt, kappaFactNoAtt, muFactNoAtt;
    bool doKappa;
    double Qmu = mQmu.sum() / 4.;
    double Qkp = mQkp.sum() / 4.;
    attBuild.computeFactors(Qmu, Qkp, alpha, beta, gamma, dKappaFact, dMuFact, 
        kappaFactAtt, muFactAtt, kappaFactNoAtt, muFactNoAtt, doKappa);
    // create attenuation
    if (attBuild.useCG4()) {
        const RDRow4 &weights_cg4 = mMyQuad->computeWeightsCG4();
        int nr = mMyQuad->getNr();
        RDMatX4 dkappa(nr, 4), dmu(nr, 4);
        dkappa.col(0) = weights_cg4(0) * dKappaFact * kappa.col(nPntEdge * 1 + 1);
        dkappa.col(1) = weights_cg4(1) * dKappaFact * kappa.col(nPntEdge * 1 + 3);
        dkappa.col(2) = weights_cg4(2) * dKappaFact * kappa.col(nPntEdge * 3 + 1);
        dkappa.col(3) = weights_cg4(3) * dKappaFact * kappa.col(nPntEdge * 3 + 3);
        dmu.col(0) = weights_cg4(0) * dMuFact * mu.col(nPntEdge * 1 + 1);
        dmu.col(1) = weights_cg4(1) * dMuFact * mu.col(nPntEdge * 1 + 3);
        dmu.col(2) = weights_cg4(2) * dMuFact * mu.col(nPntEdge * 3 + 1);
        dmu.col(3) = weights_cg4(3) * dMuFact * mu.col(nPntEdge * 3 + 3);
        kappa *= kappaFactNoAtt;
        mu *= muFactNoAtt;
        kappa.col(nPntEdge * 1 + 1) *= 1. + weights_cg4(0) * (kappaFactAtt / kappaFactNoAtt - 1.);
        kappa.col(nPntEdge * 1 + 3) *= 1. + weights_cg4(1) * (kappaFactAtt / kappaFactNoAtt - 1.);
        kappa.col(nPntEdge * 3 + 1) *= 1. + weights_cg4(2) * (kappaFactAtt / kappaFactNoAtt - 1.);
        kappa.col(nPntEdge * 3 + 3) *= 1. + weights_cg4(3) * (kappaFactAtt / kappaFactNoAtt - 1.);
        mu.col(nPntEdge * 1 + 1) *= 1. + weights_cg4(0) * (muFactAtt / muFactNoAtt - 1.);
        mu.col(nPntEdge * 1 + 3) *= 1. + weights_cg4(1) * (muFactAtt / muFactNoAtt - 1.);
        mu.col(nPntEdge * 3 + 1) *= 1. + weights_cg4(2) * (muFactAtt / muFactNoAtt - 1.);
        mu.col(nPntEdge * 3 + 3) *= 1. + weights_cg4(3) * (muFactAtt / muFactNoAtt - 1.);
        att = new Attenuation3D_CG4(attBuild.getNSLS(), 
            alpha.cast<Real>(), beta.cast<Real>(), gamma.cast<Real>(), 
            XMath::castToSolver<RDMatX4, RMatX4>(dkappa), 
            XMath::castToSolver<RDMatX4, RMatX4>(dmu), doKappa);
    } else {
        const RDMatXN &dkappa = dKappaFact * kappa;
        const RDMatXN &dmu = dMuFact * mu;
        kappa *= kappaFactAtt;
        mu *= muFactAtt;
        att = new Attenuation3D_Full(attBuild.getNSLS(), 
            alpha.cast<Real>(), beta.cast<Real>(), gamma.cast<Real>(), 
            XMath::castToSolver<RDMatXN, RMatXN>(dkappa), 
            XMath::castToSolver<RDMatXN, RMatXN>(dmu), doKappa);
    }        
}

