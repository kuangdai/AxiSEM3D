// Material.cpp
// created by Kuangdai on 17-May-2016 
// 3D seismic material properties

#include "Material.h"
#include "Quad.h"
#include "ExodusModel.h"
#include "SpectralConstants.h"

#include "Volumetric3D.h"
#include "XMath.h"
#include "Relabelling.h"

#include "Acoustic1D.h"
#include "Acoustic3D.h"

#include "AttBuilder.h"
#include "Attenuation1D.h"
#include "Attenuation3D.h"

#include "Isotropic1D.h"
#include "Isotropic3D.h"
#include "TransverselyIsotropic1D.h"
#include "TransverselyIsotropic3D.h"

#include "Anisotropic1D.h"
#include "Anisotropic3D.h"
#include "Geodesy.h"

#include <boost/algorithm/string.hpp>
#include "SlicePlot.h"

Material::Material(const Quad *myQuad, const ExodusModel &exModel): mMyQuad(myQuad) {
    // read Exodus model
    int quadTag = mMyQuad->getQuadTag();
    if (exModel.isIsotropic()) {
        for (int i = 0; i < 4; i++) {
            mVpv1D(i) = mVph1D(i) = exModel.getElementalVariables("VP_" + std::to_string(i), quadTag);
            mVsv1D(i) = mVsh1D(i) = exModel.getElementalVariables("VS_" + std::to_string(i), quadTag);
            mEta1D(i) = 1.;
        }
    } else {
        for (int i = 0; i < 4; i++) {
            mVpv1D(i) = exModel.getElementalVariables("VPV_" + std::to_string(i), quadTag);
            mVph1D(i) = exModel.getElementalVariables("VPH_" + std::to_string(i), quadTag);
            mVsv1D(i) = exModel.getElementalVariables("VSV_" + std::to_string(i), quadTag);
            mVsh1D(i) = exModel.getElementalVariables("VSH_" + std::to_string(i), quadTag);
            mEta1D(i) = exModel.getElementalVariables("ETA_" + std::to_string(i), quadTag);
        }
    }
    for (int i = 0; i < 4; i++) {
        mRho1D(i) = exModel.getElementalVariables("RHO_" + std::to_string(i), quadTag);
    }
    if (exModel.hasAttenuation()) {
        for (int i = 0; i < 4; i++) {
            mQkp1D(i) = exModel.getElementalVariables("QKAPPA_" + std::to_string(i), quadTag);
            mQmu1D(i) = exModel.getElementalVariables("QMU_" + std::to_string(i), quadTag);
        }
    } else {
        mQkp1D.setZero();
        mQmu1D.setZero();
    }
    
    // initialize 3D properties with 1D reference
    int Nr = mMyQuad->getNr();
    mVpv3D = RDMatXN::Zero(1, nPE);
    mVph3D = RDMatXN::Zero(1, nPE);
    mVsv3D = RDMatXN::Zero(1, nPE);
    mVsh3D = RDMatXN::Zero(1, nPE);
    mRho3D = RDMatXN::Zero(1, nPE);
    mEta3D = RDMatXN::Zero(1, nPE);
    mQkp3D = RDMatXN::Zero(1, nPE);
    mQmu3D = RDMatXN::Zero(1, nPE);
    
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
            mEta3D.col(ipnt).fill(Mapping::interpolate(mEta1D, xieta));
            mQkp3D.col(ipnt).fill(Mapping::interpolate(mQkp1D, xieta));
            mQmu3D.col(ipnt).fill(Mapping::interpolate(mQmu1D, xieta));
            // rho for mass
            int NrP = mMyQuad->getPointNr(ipol, jpol);
            mRhoMass3D[ipnt] = RDColX::Constant(NrP, mRho3D.col(ipnt)(0));
            mVpFluid3D[ipnt] = RDColX::Constant(NrP, mVpv3D.col(ipnt)(0));
        }
    }
}

void Material::addVolumetric3D(const std::vector<Volumetric3D *> &m3D, 
    double srcLat, double srcLon, double srcDep, double phi2D) {
    if (m3D.size() == 0) {
        return;
    }    
        
    // pointers for fast access to material matrices
    std::vector<RDRow4 *>  prop1DPtr = {&mVpv1D, &mVph1D, &mVsv1D, &mVsh1D, &mRho1D, &mEta1D, &mQkp1D, &mQmu1D,
                                        &mVpv1D, &mVsv1D, // these two take no other effect than occupying the slots  
                                        &mC11_1D, &mC12_1D, &mC13_1D, &mC14_1D, &mC15_1D, &mC16_1D,
                                        &mC22_1D, &mC23_1D, &mC24_1D, &mC25_1D, &mC26_1D,
                                        &mC33_1D, &mC34_1D, &mC35_1D, &mC36_1D,
                                        &mC44_1D, &mC45_1D, &mC46_1D,
                                        &mC55_1D, &mC56_1D,
                                        &mC66_1D
    };
    std::vector<RDMatXN *> prop3DPtr = {&mVpv3D, &mVph3D, &mVsv3D, &mVsh3D, &mRho3D, &mEta3D, &mQkp3D, &mQmu3D,
                                        &mVpv3D, &mVsv3D, // these two take no other effect than occupying the slots  
                                        &mC11_3D, &mC12_3D, &mC13_3D, &mC14_3D, &mC15_3D, &mC16_3D,
                                        &mC22_3D, &mC23_3D, &mC24_3D, &mC25_3D, &mC26_3D,
                                        &mC33_3D, &mC34_3D, &mC35_3D, &mC36_3D,
                                        &mC44_3D, &mC45_3D, &mC46_3D,
                                        &mC55_3D, &mC56_3D,
                                        &mC66_3D
    };
    
    // radius at element center 
    double rElemCenter = mMyQuad->computeCenterRadius();
    
    // read 3D model
    int Nr = mMyQuad->getNr();
    for (int ipol = 0; ipol <= nPol; ipol++) {
        for (int jpol = 0; jpol <= nPol; jpol++) {
            // geographic oordinates of cardinal points
            const RDCol2 &xieta = SpectralConstants::getXiEta(ipol, jpol, mMyQuad->isAxial());
            const RDMatX3 &rtp = mMyQuad->computeGeocentricGlobal(srcLat, srcLon, srcDep, xieta, Nr, phi2D);
            int ipnt = ipol * nPntEdge + jpol;
            for (int alpha = 0; alpha < Nr; alpha++) {
                double r = rtp(alpha, 0);
                double t = rtp(alpha, 1);
                double p = rtp(alpha, 2);
                for (const auto &model: m3D) {
                    if (mMyQuad->isFluid() && !model->makeFluid3D()) {
                        continue;
                    }
                    std::vector<Volumetric3D::MaterialProperty> properties; 
                    std::vector<Volumetric3D::MaterialRefType> refTypes;
                    std::vector<double> values;
                    if (!model->get3dProperties(r, t, p, rElemCenter, properties, refTypes, values)) {
                        // point (r, t, p) not in model range
                        continue;
                    }
                    
                    if (!_3Dprepared()) {
                        prepare3D();
                    }
                    
                    // deal with VP and VS
                    std::vector<Volumetric3D::MaterialProperty> propertiesTIso; 
                    std::vector<Volumetric3D::MaterialRefType> refTypesTIso;
                    std::vector<double> valuesTIso;
                    for (int iprop = 0; iprop < properties.size(); iprop++) {
                        if (properties[iprop] == Volumetric3D::MaterialProperty::VP) {
                            propertiesTIso.push_back(Volumetric3D::MaterialProperty::VPV);
                            propertiesTIso.push_back(Volumetric3D::MaterialProperty::VPH);
                            refTypesTIso.push_back(refTypes[iprop]);
                            refTypesTIso.push_back(refTypes[iprop]);
                            valuesTIso.push_back(values[iprop]);
                            valuesTIso.push_back(values[iprop]);
                        } else if (properties[iprop] == Volumetric3D::MaterialProperty::VS) {
                            propertiesTIso.push_back(Volumetric3D::MaterialProperty::VSV);
                            propertiesTIso.push_back(Volumetric3D::MaterialProperty::VSH);
                            refTypesTIso.push_back(refTypes[iprop]);
                            refTypesTIso.push_back(refTypes[iprop]);
                            valuesTIso.push_back(values[iprop]);
                            valuesTIso.push_back(values[iprop]);
                        } else {
                            propertiesTIso.push_back(properties[iprop]);
                            refTypesTIso.push_back(refTypes[iprop]);
                            valuesTIso.push_back(values[iprop]);
                        }
                    }
                    // change values
                    for (int iprop = 0; iprop < propertiesTIso.size(); iprop++) {
                        // initialize anisotropy
                        if (!mFullAniso && propertiesTIso[iprop] >= Volumetric3D::MaterialProperty::C11) {
                            initAniso();
                        }
                        
                        // // check
                        // if (mFullAniso) {
                        //     if (propertiesTIso[iprop] <= Volumetric3D::MaterialProperty::ANIS_ETA) {
                        //         throw std::runtime_error("Material::addVolumetric3D || " 
                        //             "Velocity, density and eta can no longer be changed once "
                        //             "full anisotropy has been activated.");
                        //     }
                        // } 
                        
                        RDRow4 &row1D = *prop1DPtr[propertiesTIso[iprop]];
                        RDMatXN &mat3D = *prop3DPtr[propertiesTIso[iprop]];
                        Volumetric3D::MaterialRefType ref_type = refTypesTIso[iprop];
                        double value3D = valuesTIso[iprop];
                        if (ref_type == Volumetric3D::MaterialRefType::Absolute) {
                            mat3D(alpha, ipnt) = value3D;
                        } else if (ref_type == Volumetric3D::MaterialRefType::Reference1D) {
                            double ref1D = Mapping::interpolate(row1D, xieta);
                            mat3D(alpha, ipnt) = ref1D * (1. + value3D);
                        } else if (ref_type == Volumetric3D::MaterialRefType::Reference3D) {
                            mat3D(alpha, ipnt) *= 1. + value3D;
                        } else {
                            double ref1D = Mapping::interpolate(row1D, xieta);
                            mat3D(alpha, ipnt) = (mat3D(alpha, ipnt) - ref1D) * (1. + value3D) + ref1D;
                        }
                    }
                }
            }
        }
    }
    
    // form mass point sampling
    if (_3Dprepared()) {
        for (int ipol = 0; ipol <= nPol; ipol++) {
            for (int jpol = 0; jpol <= nPol; jpol++) {
                int ipnt = ipol * nPntEdge + jpol;
                int nr_mass = mMyQuad->getPointNr(ipol, jpol);
                mRhoMass3D[ipnt] = XMath::linearResampling(nr_mass, mRho3D.col(ipnt));
                mVpFluid3D[ipnt] = XMath::linearResampling(nr_mass, mVpv3D.col(ipnt));
            }
        }
    }
    
    // rotate anisotropy from geographic to source-centred
    if (mFullAniso) {
        rotateAniso(srcLat, srcLon, srcDep);
    }
}

arPP_RDColX Material::computeElementalMass() const {
    arPP_RDColX mass, J;
    // Jacobian of topography
    if (mMyQuad->hasRelabelling()) {
        J = mMyQuad->getRelabelling().getMassJacobian();
    } else {
        J = mRhoMass3D; // just to use the size
        for (int ipnt = 0; ipnt < nPntElem; ipnt++) {
            J[ipnt].setOnes();
        }
    }
    // general mass term
    const RDRowN &iFact = mMyQuad->getIntegralFactor();
    for (int ipnt = 0; ipnt < nPntElem; ipnt++) {
        if (mMyQuad->isFluid()) {
            mass[ipnt] = (mRhoMass3D[ipnt].array() * mVpFluid3D[ipnt].array().pow(2.)).pow(-1.).matrix();
        } else {
            mass[ipnt] = mRhoMass3D[ipnt];
        }
        mass[ipnt].array() *= iFact(ipnt) * J[ipnt].array();
    }
    return mass;
}

Acoustic *Material::createAcoustic(bool elem1D) const {
    const RDRowN &iFact = mMyQuad->getIntegralFactor();
    RDMatXN fluidK;
    if (_3Dprepared()) {
        fluidK = mRho3D.array().pow(-1.);    
    } else {
        fluidK = mRho3D.replicate(mMyQuad->getNr(), 1).array().pow(-1.);
    }
    for (int ipnt = 0; ipnt < nPntElem; ipnt++) {
       fluidK.col(ipnt) *= iFact(ipnt);
    }
    if (mMyQuad->hasRelabelling()) {
        fluidK.array() *= mMyQuad->getRelabelling().getStiffJacobian().array();
    }
    if (elem1D) {
        RDMatPP kstruct;
        XMath::structuredUseFirstRow(fluidK, kstruct);
        return new Acoustic1D(kstruct.cast<Real>());
    } else {
        return new Acoustic3D(fluidK.cast<Real>());
    }
}

Elastic *Material::createElastic(bool elem1D, const AttBuilder *attBuild) const {
    if (mFullAniso) {
        return createElasticAniso(elem1D, attBuild);
    }
    
    // elasticity tensor
    const RDRowN &iFact = mMyQuad->getIntegralFactor(); 
    RDMatXN vpv3D, vph3D;
    RDMatXN vsv3D, vsh3D;
    RDMatXN rho3D;
    RDMatXN eta3D;
    RDMatXN qkp3D, qmu3D;
    
    if (_3Dprepared()) {
        vpv3D = mVpv3D;
        vph3D = mVph3D;
        vsv3D = mVsv3D;
        vsh3D = mVsh3D;
        rho3D = mRho3D;
        eta3D = mEta3D;
        qkp3D = mQkp3D;
        qmu3D = mQmu3D;
    } else {
        vpv3D = mVpv3D.replicate(mMyQuad->getNr(), 1);
        vph3D = mVph3D.replicate(mMyQuad->getNr(), 1);
        vsv3D = mVsv3D.replicate(mMyQuad->getNr(), 1);
        vsh3D = mVsh3D.replicate(mMyQuad->getNr(), 1);
        rho3D = mRho3D.replicate(mMyQuad->getNr(), 1);
        eta3D = mEta3D.replicate(mMyQuad->getNr(), 1);
        qkp3D = mQkp3D.replicate(mMyQuad->getNr(), 1);
        qmu3D = mQmu3D.replicate(mMyQuad->getNr(), 1);
    }
    
    RDMatXN A(rho3D), C(rho3D), F(rho3D), L(rho3D), N(rho3D);
    for (int ipnt = 0; ipnt < nPntElem; ipnt++) {
        // A C L N
        A.col(ipnt).array() *= vph3D.col(ipnt).array().pow(2.) * iFact(ipnt);
        C.col(ipnt).array() *= vpv3D.col(ipnt).array().pow(2.) * iFact(ipnt);
        L.col(ipnt).array() *= vsv3D.col(ipnt).array().pow(2.) * iFact(ipnt);
        N.col(ipnt).array() *= vsh3D.col(ipnt).array().pow(2.) * iFact(ipnt);
    }
    // F
    F = eta3D.schur(A - 2. * L);
    // must do relabelling before attenuation
    if (mMyQuad->hasRelabelling()) {
        const RDMatXN &J = mMyQuad->getRelabelling().getStiffJacobian();
        A = A.schur(J);
        C = C.schur(J);
        F = F.schur(J);
        L = L.schur(J);
        N = N.schur(J);
    }
    
    // attenuation
    Attenuation1D *att1D = 0;
    Attenuation3D *att3D = 0;
    if (attBuild) {
        // Voigt average
        RDMatXN kappa = (4. * A + C + 4. * F - 4. * N) / 9.;
        RDMatXN mu = (A + C - 2. * F + 6. * L + 5. * N) / 15.;
        A -= (kappa + 4. / 3. * mu);
        C -= (kappa + 4. / 3. * mu);
        F -= (kappa - 2. / 3. * mu);
        L -= mu;
        N -= mu;
        if (elem1D) {
            att1D = attBuild->createAttenuation1D(qkp3D, qmu3D, kappa, mu, mMyQuad);
        } else {
            att3D = attBuild->createAttenuation3D(qkp3D, qmu3D, kappa, mu, mMyQuad);
        }
        A += (kappa + 4. / 3. * mu);
        C += (kappa + 4. / 3. * mu);
        F += (kappa - 2. / 3. * mu);
        L += mu;
        N += mu;
    }
    
    // Elastic pointers
    if (elem1D) {
        RDMatPP A0, C0, F0, L0, N0;
        XMath::structuredUseFirstRow(A, A0);
        XMath::structuredUseFirstRow(C, C0);
        XMath::structuredUseFirstRow(F, F0);
        XMath::structuredUseFirstRow(L, L0);
        XMath::structuredUseFirstRow(N, N0);
        if (isIsotropic()) {
            return new Isotropic1D(F0.cast<Real>(), L0.cast<Real>(), att1D);
        } else {
            return new TransverselyIsotropic1D(A0.cast<Real>(), C0.cast<Real>(), 
                F0.cast<Real>(), L0.cast<Real>(), N0.cast<Real>(), att1D);
        }    
    } else {
        if (isIsotropic()) {
            return new Isotropic3D(F.cast<Real>(), L.cast<Real>(), att3D);
        } else {
            return new TransverselyIsotropic3D(A.cast<Real>(), C.cast<Real>(), 
                F.cast<Real>(), L.cast<Real>(), N.cast<Real>(), att3D);
        }
    }
}

Elastic *Material::createElasticAniso(bool elem1D, const AttBuilder *attBuild) const {
    // elasticity tensor
    const RDRowN &iFact = mMyQuad->getIntegralFactor(); 
    RDMatXN C11_3D(mC11_3D), C12_3D(mC12_3D), C13_3D(mC13_3D), C14_3D(mC14_3D), C15_3D(mC15_3D), C16_3D(mC16_3D);
    RDMatXN C22_3D(mC22_3D), C23_3D(mC23_3D), C24_3D(mC24_3D), C25_3D(mC25_3D), C26_3D(mC26_3D);
    RDMatXN C33_3D(mC33_3D), C34_3D(mC34_3D), C35_3D(mC35_3D), C36_3D(mC36_3D);
    RDMatXN C44_3D(mC44_3D), C45_3D(mC45_3D), C46_3D(mC46_3D);
    RDMatXN C55_3D(mC55_3D), C56_3D(mC56_3D);
    RDMatXN C66_3D(mC66_3D);
    for (int ipnt = 0; ipnt < nPntElem; ipnt++) {
        C11_3D.col(ipnt) *= iFact(ipnt);
        C12_3D.col(ipnt) *= iFact(ipnt);
        C13_3D.col(ipnt) *= iFact(ipnt);
        C14_3D.col(ipnt) *= iFact(ipnt);
        C15_3D.col(ipnt) *= iFact(ipnt);
        C16_3D.col(ipnt) *= iFact(ipnt);
        C22_3D.col(ipnt) *= iFact(ipnt);
        C23_3D.col(ipnt) *= iFact(ipnt);
        C24_3D.col(ipnt) *= iFact(ipnt);
        C25_3D.col(ipnt) *= iFact(ipnt);
        C26_3D.col(ipnt) *= iFact(ipnt);
        C33_3D.col(ipnt) *= iFact(ipnt);
        C34_3D.col(ipnt) *= iFact(ipnt);
        C35_3D.col(ipnt) *= iFact(ipnt);
        C36_3D.col(ipnt) *= iFact(ipnt);
        C44_3D.col(ipnt) *= iFact(ipnt);
        C45_3D.col(ipnt) *= iFact(ipnt);
        C46_3D.col(ipnt) *= iFact(ipnt);
        C55_3D.col(ipnt) *= iFact(ipnt);
        C56_3D.col(ipnt) *= iFact(ipnt);
        C66_3D.col(ipnt) *= iFact(ipnt);
    }
    // must do relabelling before attenuation
    if (mMyQuad->hasRelabelling()) {
        const RDMatXN &J = mMyQuad->getRelabelling().getStiffJacobian();
        C11_3D = C11_3D.schur(J);
        C12_3D = C12_3D.schur(J);
        C13_3D = C13_3D.schur(J);
        C14_3D = C14_3D.schur(J);
        C15_3D = C15_3D.schur(J);
        C16_3D = C16_3D.schur(J);
        C22_3D = C22_3D.schur(J);
        C23_3D = C23_3D.schur(J);
        C24_3D = C24_3D.schur(J);
        C25_3D = C25_3D.schur(J);
        C26_3D = C26_3D.schur(J);
        C33_3D = C33_3D.schur(J);
        C34_3D = C34_3D.schur(J);
        C35_3D = C35_3D.schur(J);
        C36_3D = C36_3D.schur(J);
        C44_3D = C44_3D.schur(J);
        C45_3D = C45_3D.schur(J);
        C46_3D = C46_3D.schur(J);
        C55_3D = C55_3D.schur(J);
        C56_3D = C56_3D.schur(J);
        C66_3D = C66_3D.schur(J);
    }
    
    // attenuation
    Attenuation1D *att1D = 0;
    Attenuation3D *att3D = 0;
    if (attBuild) {
        // Voigt average
        // https://materialsproject.org/wiki/index.php/Elasticity_calculations
        RDMatXN kappa = (C11_3D + C22_3D + C33_3D + 2. * (C12_3D + C23_3D + C13_3D)) / 9.;
        RDMatXN mu = (C11_3D + C22_3D + C33_3D - (C12_3D + C23_3D + C13_3D) + 3. * (C44_3D + C55_3D + C66_3D)) / 15.;
        C11_3D -= (kappa + 4. / 3. * mu);
        C22_3D -= (kappa + 4. / 3. * mu);
        C33_3D -= (kappa + 4. / 3. * mu);
        C12_3D -= (kappa - 2. / 3. * mu);
        C23_3D -= (kappa - 2. / 3. * mu);
        C13_3D -= (kappa - 2. / 3. * mu);
        C44_3D -= mu;
        C55_3D -= mu;
        C66_3D -= mu;
        if (elem1D) {
            att1D = attBuild->createAttenuation1D(mQkp3D, mQmu3D, kappa, mu, mMyQuad);
        } else {
            att3D = attBuild->createAttenuation3D(mQkp3D, mQmu3D, kappa, mu, mMyQuad);
        }
        C11_3D += (kappa + 4. / 3. * mu);
        C22_3D += (kappa + 4. / 3. * mu);
        C33_3D += (kappa + 4. / 3. * mu);
        C12_3D += (kappa - 2. / 3. * mu);
        C23_3D += (kappa - 2. / 3. * mu);
        C13_3D += (kappa - 2. / 3. * mu);
        C44_3D += mu;
        C55_3D += mu;
        C66_3D += mu;
    }
    
    // Elastic pointers
    if (elem1D) {
        RDMatPP C11_1D, C12_1D, C13_1D, C14_1D, C15_1D, C16_1D;
        RDMatPP C22_1D, C23_1D, C24_1D, C25_1D, C26_1D;
        RDMatPP C33_1D, C34_1D, C35_1D, C36_1D;
        RDMatPP C44_1D, C45_1D, C46_1D;
        RDMatPP C55_1D, C56_1D;
        RDMatPP C66_1D;
        XMath::structuredUseFirstRow(C11_3D, C11_1D);
        XMath::structuredUseFirstRow(C12_3D, C12_1D);
        XMath::structuredUseFirstRow(C13_3D, C13_1D);
        XMath::structuredUseFirstRow(C14_3D, C14_1D);
        XMath::structuredUseFirstRow(C15_3D, C15_1D);
        XMath::structuredUseFirstRow(C16_3D, C16_1D);
        XMath::structuredUseFirstRow(C22_3D, C22_1D);
        XMath::structuredUseFirstRow(C23_3D, C23_1D);
        XMath::structuredUseFirstRow(C24_3D, C24_1D);
        XMath::structuredUseFirstRow(C25_3D, C25_1D);
        XMath::structuredUseFirstRow(C26_3D, C26_1D);
        XMath::structuredUseFirstRow(C33_3D, C33_1D);
        XMath::structuredUseFirstRow(C34_3D, C34_1D);
        XMath::structuredUseFirstRow(C35_3D, C35_1D);
        XMath::structuredUseFirstRow(C36_3D, C36_1D);
        XMath::structuredUseFirstRow(C44_3D, C44_1D);
        XMath::structuredUseFirstRow(C45_3D, C45_1D);
        XMath::structuredUseFirstRow(C46_3D, C46_1D);
        XMath::structuredUseFirstRow(C55_3D, C55_1D);
        XMath::structuredUseFirstRow(C56_3D, C56_1D);
        XMath::structuredUseFirstRow(C66_3D, C66_1D);
        return new Anisotropic1D(
            C11_1D.cast<Real>(), C12_1D.cast<Real>(), C13_1D.cast<Real>(), C14_1D.cast<Real>(), C15_1D.cast<Real>(), C16_1D.cast<Real>(),
            C22_1D.cast<Real>(), C23_1D.cast<Real>(), C24_1D.cast<Real>(), C25_1D.cast<Real>(), C26_1D.cast<Real>(),
            C33_1D.cast<Real>(), C34_1D.cast<Real>(), C35_1D.cast<Real>(), C36_1D.cast<Real>(),
            C44_1D.cast<Real>(), C45_1D.cast<Real>(), C46_1D.cast<Real>(),
            C55_1D.cast<Real>(), C56_1D.cast<Real>(), 
            C66_1D.cast<Real>(), 
            att1D);
    } else {
        return new Anisotropic3D(
            C11_3D.cast<Real>(), C12_3D.cast<Real>(), C13_3D.cast<Real>(), C14_3D.cast<Real>(), C15_3D.cast<Real>(), C16_3D.cast<Real>(),
            C22_3D.cast<Real>(), C23_3D.cast<Real>(), C24_3D.cast<Real>(), C25_3D.cast<Real>(), C26_3D.cast<Real>(),
            C33_3D.cast<Real>(), C34_3D.cast<Real>(), C35_3D.cast<Real>(), C36_3D.cast<Real>(),
            C44_3D.cast<Real>(), C45_3D.cast<Real>(), C46_3D.cast<Real>(),
            C55_3D.cast<Real>(), C56_3D.cast<Real>(), 
            C66_3D.cast<Real>(), 
            att3D);
    }
}

double Material::getVMaxRef() const {
    return std::max(mVph1D.maxCoeff(), mVpv1D.maxCoeff());
}

RDColX Material::getVMax() const {
    RDMatXN vpv3D, vph3D;
    if (_3Dprepared()) {
        vpv3D = mVpv3D;
        vph3D = mVph3D;
    } else {
        vpv3D = mVpv3D.replicate(mMyQuad->getNr(), 1);
        vph3D = mVph3D.replicate(mMyQuad->getNr(), 1);
    }
    const RDColX &vpvMax = vpv3D.rowwise().maxCoeff();
    const RDColX &vphMax = vph3D.rowwise().maxCoeff();
    return (vpvMax.array().max(vphMax.array())).matrix();
}

bool Material::isFluidPar1D() const {
    return XMath::equalRows(mVpv3D) && XMath::equalRows(mRho3D);
}

bool Material::isSolidPar1D(bool attenuation) const {
    bool result = false;
    if (mFullAniso) {
        result = XMath::equalRows(mC11_3D, 1e-7) && XMath::equalRows(mC12_3D, 1e-7) && XMath::equalRows(mC13_3D, 1e-7) && XMath::equalRows(mC14_3D, 1e-7) && XMath::equalRows(mC15_3D, 1e-7) && XMath::equalRows(mC16_3D, 1e-7) &&
                 XMath::equalRows(mC22_3D, 1e-7) && XMath::equalRows(mC23_3D, 1e-7) && XMath::equalRows(mC24_3D, 1e-7) && XMath::equalRows(mC25_3D, 1e-7) && XMath::equalRows(mC26_3D, 1e-7) &&
                 XMath::equalRows(mC33_3D, 1e-7) && XMath::equalRows(mC34_3D, 1e-7) && XMath::equalRows(mC35_3D, 1e-7) && XMath::equalRows(mC36_3D, 1e-7) &&
                 XMath::equalRows(mC44_3D, 1e-7) && XMath::equalRows(mC45_3D, 1e-7) && XMath::equalRows(mC46_3D, 1e-7) &&
                 XMath::equalRows(mC55_3D, 1e-7) && XMath::equalRows(mC56_3D, 1e-7) &&
                 XMath::equalRows(mC66_3D, 1e-7) && 
                 XMath::equalRows(mRho3D);
    } else {
        result = XMath::equalRows(mVpv3D) && XMath::equalRows(mVph3D) &&
                 XMath::equalRows(mVsv3D) && XMath::equalRows(mVsh3D) && 
                 XMath::equalRows(mRho3D) && XMath::equalRows(mEta3D);
    }
    if (attenuation) {
        result = result && XMath::equalRows(mQkp3D) && XMath::equalRows(mQmu3D);
    }              
    return result;
}

bool Material::isIsotropic() const {
    if (mFullAniso) {
        return false;
    }
    return (mVpv3D - mVph3D).norm() < tinyDouble * mVpv3D.norm() &&  
           (mVsv3D - mVsh3D).norm() < tinyDouble * mVsv3D.norm() && 
           (mEta3D - RDMatXN::Ones(mEta3D.rows(), mEta3D.cols())).norm() < tinyDouble;
}

RDMatXN Material::getProperty(const std::string &vname, int refType) {
    int Nr = mMyQuad->getNr();
    
    // name
    std::string varname = vname;
    if (boost::iequals(vname, "vp")) {
        varname = "vpv";
    }
    if (boost::iequals(vname, "vs")) {
        varname = "vsv";
    }
    RDRow4 data1D;
    RDMatXN data3D;
    if (boost::iequals(varname, "vpv")) {
        data1D = mVpv1D;
        data3D = mVpv3D;
    } else if (boost::iequals(varname, "vsv")) {
        data1D = mVsv1D;
        data3D = mVsv3D;
    } else if (boost::iequals(varname, "vph")) {
        data1D = mVph1D;
        data3D = mVph3D;
    } else if (boost::iequals(varname, "vsh")) {
        data1D = mVsh1D;
        data3D = mVsh3D;
    } else if (boost::iequals(varname, "rho")) {
        data1D = mRho1D;
        data3D = mRho3D;
    } else if (boost::iequals(varname, "eta")) {
        data1D = mEta1D;
        data3D = mEta3D;
    } else if (boost::iequals(varname, "qkappa")) {
        data1D = mQkp1D;
        data3D = mQkp3D;
    } else if (boost::iequals(varname, "qmu")) {
        data1D = mQmu1D;
        data3D = mQmu3D;
    } else if (mFullAniso) {
        if (boost::iequals(varname, "c11")) {
            data1D = mC11_1D;
            data3D = mC11_3D;
        } else if (boost::iequals(varname, "c12")) {
            data1D = mC12_1D;
            data3D = mC12_3D;
        } else if (boost::iequals(varname, "c13")) {
            data1D = mC13_1D;
            data3D = mC13_3D;
        } else if (boost::iequals(varname, "c14")) {
            data1D = mC14_1D;
            data3D = mC14_3D;
        } else if (boost::iequals(varname, "c15")) {
            data1D = mC15_1D;
            data3D = mC15_3D;
        } else if (boost::iequals(varname, "c16")) {
            data1D = mC16_1D;
            data3D = mC16_3D;
        } else if (boost::iequals(varname, "c22")) {
            data1D = mC22_1D;
            data3D = mC22_3D;
        } else if (boost::iequals(varname, "c23")) {
            data1D = mC23_1D;
            data3D = mC23_3D;
        } else if (boost::iequals(varname, "c24")) {
            data1D = mC24_1D;
            data3D = mC24_3D;
        } else if (boost::iequals(varname, "c25")) {
            data1D = mC25_1D;
            data3D = mC25_3D;
        } else if (boost::iequals(varname, "c26")) {
            data1D = mC26_1D;
            data3D = mC26_3D;
        } else if (boost::iequals(varname, "c33")) {
            data1D = mC33_1D;
            data3D = mC33_3D;
        } else if (boost::iequals(varname, "c34")) {
            data1D = mC34_1D;
            data3D = mC34_3D;
        } else if (boost::iequals(varname, "c35")) {
            data1D = mC35_1D;
            data3D = mC35_3D;
        } else if (boost::iequals(varname, "c36")) {
            data1D = mC36_1D;
            data3D = mC36_3D;
        } else if (boost::iequals(varname, "c44")) {
            data1D = mC44_1D;
            data3D = mC44_3D;
        } else if (boost::iequals(varname, "c45")) {
            data1D = mC45_1D;
            data3D = mC45_3D;
        } else if (boost::iequals(varname, "c46")) {
            data1D = mC46_1D;
            data3D = mC46_3D;
        } else if (boost::iequals(varname, "c55")) {
            data1D = mC55_1D;
            data3D = mC55_3D;
        } else if (boost::iequals(varname, "c56")) {
            data1D = mC56_1D;
            data3D = mC56_3D;
        } else if (boost::iequals(varname, "c66")) {
            data1D = mC66_1D;
            data3D = mC66_3D;
        } 
    } else {
        throw std::runtime_error("Material::getProperty || Unknown field variable name: " + vname);
    }
    
    if (data3D.rows() != Nr) {
        data3D = data3D.replicate(Nr, 1);
    }
    
    // 3D
    if (refType == SlicePlot::PropertyRefTypes::Property3D) {
        return data3D;
    }
    
    // fill 1D
    RDMatXN data1DXN(Nr, nPntElem);
    for (int ipol = 0; ipol <= nPol; ipol++) {
        for (int jpol = 0; jpol <= nPol; jpol++) {
            int ipnt = ipol * nPntEdge + jpol;
            const RDCol2 &xieta = SpectralConstants::getXiEta(ipol, jpol, mMyQuad->isAxial());
            data1DXN.col(ipnt).fill(Mapping::interpolate(data1D, xieta));
        }
    }
    
    // 1D
    if (refType == SlicePlot::PropertyRefTypes::Property1D) {
        return data1DXN;
    }
    
    // perturb
    RDMatXN data1DBase = data1DXN.array().max(tinyDouble).matrix(); // in fluid, vs = 0
    return ((data3D - data1DXN).array() / data1DBase.array()).matrix();
}

void Material::initAniso() {
    if (mMyQuad->isFluid()) {
        throw std::runtime_error("Material::initAniso || Cannot activate full anisotropy in fluid domain.");
    }
    
    // 1D elasticity tensor
    RDRow4 A_1D = mRho1D.schur(mVph1D).schur(mVph1D);
    RDRow4 C_1D = mRho1D.schur(mVpv1D).schur(mVpv1D);
    RDRow4 L_1D = mRho1D.schur(mVsv1D).schur(mVsv1D);
    RDRow4 N_1D = mRho1D.schur(mVsh1D).schur(mVsh1D);
    RDRow4 F_1D = mEta1D.schur(A_1D - 2. * L_1D);
    
    mC11_1D = mC22_1D = A_1D;
    mC33_1D = C_1D;
    mC44_1D = mC55_1D = L_1D;
    mC66_1D = N_1D;
    mC12_1D = A_1D - 2. * N_1D;
    mC13_1D = mC23_1D = F_1D;
    mC14_1D = mC15_1D = mC16_1D = RDRow4::Zero();
    mC24_1D = mC25_1D = mC26_1D = RDRow4::Zero();
    mC34_1D = mC35_1D = mC36_1D = RDRow4::Zero();
    mC45_1D = mC46_1D = mC56_1D = RDRow4::Zero();
    
    // 3D elasticity tensor
    RDMatXN A_3D = mRho3D.schur(mVph3D).schur(mVph3D);
    RDMatXN C_3D = mRho3D.schur(mVpv3D).schur(mVpv3D);
    RDMatXN L_3D = mRho3D.schur(mVsv3D).schur(mVsv3D);
    RDMatXN N_3D = mRho3D.schur(mVsh3D).schur(mVsh3D);
    RDMatXN F_3D = mEta3D.schur(A_3D - 2. * L_3D);
    
    mC11_3D = mC22_3D = A_3D;
    mC33_3D = C_3D;
    mC44_3D = mC55_3D = L_3D;
    mC66_3D = N_3D;
    mC12_3D = A_3D - 2. * N_3D;
    mC13_3D = mC23_3D = F_3D;
    mC14_3D = mC15_3D = mC16_3D = RDMatXN::Zero(A_3D.rows(), A_3D.cols());
    mC24_3D = mC25_3D = mC26_3D = RDMatXN::Zero(A_3D.rows(), A_3D.cols());
    mC34_3D = mC35_3D = mC36_3D = RDMatXN::Zero(A_3D.rows(), A_3D.cols());
    mC45_3D = mC46_3D = mC56_3D = RDMatXN::Zero(A_3D.rows(), A_3D.cols());

    // initialized
    mFullAniso = true;
}

void Material::rotateAniso(double srcLat, double srcLon, double srcDep) {
    RDMatXX inCijkl(6, 6);
    
    // 3D
    for (int alpha = 0; alpha < mC11_3D.rows(); alpha++) {
        // azimuth of the slice
        double phi = 2. * pi / mC11_3D.rows() * alpha;
        double recLon = Geodesy::phi2Lon(phi);
        // loop over GLL points
        for (int ipol = 0; ipol <= nPol; ipol++) {
            for (int jpol = 0; jpol <= nPol; jpol++) {
                int ipnt = ipol * nPntEdge + jpol;

                inCijkl(1 - 1, 1 - 1) = mC11_3D(alpha, ipnt);
                inCijkl(1 - 1, 2 - 1) = mC12_3D(alpha, ipnt);
                inCijkl(1 - 1, 3 - 1) = mC13_3D(alpha, ipnt);
                inCijkl(1 - 1, 4 - 1) = mC14_3D(alpha, ipnt);
                inCijkl(1 - 1, 5 - 1) = mC15_3D(alpha, ipnt);
                inCijkl(1 - 1, 6 - 1) = mC16_3D(alpha, ipnt);
                
                inCijkl(2 - 1, 1 - 1) = mC12_3D(alpha, ipnt);
                inCijkl(2 - 1, 2 - 1) = mC22_3D(alpha, ipnt);
                inCijkl(2 - 1, 3 - 1) = mC23_3D(alpha, ipnt);
                inCijkl(2 - 1, 4 - 1) = mC24_3D(alpha, ipnt);
                inCijkl(2 - 1, 5 - 1) = mC25_3D(alpha, ipnt);
                inCijkl(2 - 1, 6 - 1) = mC26_3D(alpha, ipnt);
                
                inCijkl(3 - 1, 1 - 1) = mC13_3D(alpha, ipnt);
                inCijkl(3 - 1, 2 - 1) = mC23_3D(alpha, ipnt);
                inCijkl(3 - 1, 3 - 1) = mC33_3D(alpha, ipnt);
                inCijkl(3 - 1, 4 - 1) = mC34_3D(alpha, ipnt);
                inCijkl(3 - 1, 5 - 1) = mC35_3D(alpha, ipnt);
                inCijkl(3 - 1, 6 - 1) = mC36_3D(alpha, ipnt);
                
                inCijkl(4 - 1, 1 - 1) = mC14_3D(alpha, ipnt);
                inCijkl(4 - 1, 2 - 1) = mC24_3D(alpha, ipnt);
                inCijkl(4 - 1, 3 - 1) = mC34_3D(alpha, ipnt);
                inCijkl(4 - 1, 4 - 1) = mC44_3D(alpha, ipnt);
                inCijkl(4 - 1, 5 - 1) = mC45_3D(alpha, ipnt);
                inCijkl(4 - 1, 6 - 1) = mC46_3D(alpha, ipnt);
                
                inCijkl(5 - 1, 1 - 1) = mC15_3D(alpha, ipnt);
                inCijkl(5 - 1, 2 - 1) = mC25_3D(alpha, ipnt);
                inCijkl(5 - 1, 3 - 1) = mC35_3D(alpha, ipnt);
                inCijkl(5 - 1, 4 - 1) = mC45_3D(alpha, ipnt);
                inCijkl(5 - 1, 5 - 1) = mC55_3D(alpha, ipnt);
                inCijkl(5 - 1, 6 - 1) = mC56_3D(alpha, ipnt);
                
                inCijkl(6 - 1, 1 - 1) = mC16_3D(alpha, ipnt);
                inCijkl(6 - 1, 2 - 1) = mC26_3D(alpha, ipnt);
                inCijkl(6 - 1, 3 - 1) = mC36_3D(alpha, ipnt);
                inCijkl(6 - 1, 4 - 1) = mC46_3D(alpha, ipnt);
                inCijkl(6 - 1, 5 - 1) = mC56_3D(alpha, ipnt);
                inCijkl(6 - 1, 6 - 1) = mC66_3D(alpha, ipnt);
            
                // compute backazimuth
                const RDCol2 &xieta = SpectralConstants::getXiEta(ipol, jpol, mMyQuad->isAxial());
                RDCol2 rtheta = Geodesy::rtheta(mMyQuad->mapping(xieta));
                double recLat = Geodesy::theta2Lat_r(rtheta(1), rtheta(0));
                double recDep = Geodesy::getROuter() - rtheta(0);
                double baz = Geodesy::backAzimuth(srcLat, srcLon, srcDep, recLat, recLon, recDep);
                // Bond Transformation
                const RDMatXX &outCijkl = bondTransformation(inCijkl, 0., 0., -baz);
            
                mC11_3D(alpha, ipnt) = outCijkl(1 - 1, 1 - 1);
                mC12_3D(alpha, ipnt) = outCijkl(1 - 1, 2 - 1);
                mC13_3D(alpha, ipnt) = outCijkl(1 - 1, 3 - 1);
                mC14_3D(alpha, ipnt) = outCijkl(1 - 1, 4 - 1);
                mC15_3D(alpha, ipnt) = outCijkl(1 - 1, 5 - 1);
                mC16_3D(alpha, ipnt) = outCijkl(1 - 1, 6 - 1);
                
                mC22_3D(alpha, ipnt) = outCijkl(2 - 1, 2 - 1);
                mC23_3D(alpha, ipnt) = outCijkl(2 - 1, 3 - 1);
                mC24_3D(alpha, ipnt) = outCijkl(2 - 1, 4 - 1);
                mC25_3D(alpha, ipnt) = outCijkl(2 - 1, 5 - 1);
                mC26_3D(alpha, ipnt) = outCijkl(2 - 1, 6 - 1);
                
                mC33_3D(alpha, ipnt) = outCijkl(3 - 1, 3 - 1);
                mC34_3D(alpha, ipnt) = outCijkl(3 - 1, 4 - 1);
                mC35_3D(alpha, ipnt) = outCijkl(3 - 1, 5 - 1);
                mC36_3D(alpha, ipnt) = outCijkl(3 - 1, 6 - 1);
                
                mC44_3D(alpha, ipnt) = outCijkl(4 - 1, 4 - 1);
                mC45_3D(alpha, ipnt) = outCijkl(4 - 1, 5 - 1);
                mC46_3D(alpha, ipnt) = outCijkl(4 - 1, 6 - 1);
                
                mC55_3D(alpha, ipnt) = outCijkl(5 - 1, 5 - 1);
                mC56_3D(alpha, ipnt) = outCijkl(5 - 1, 6 - 1);
                
                mC66_3D(alpha, ipnt) = outCijkl(6 - 1, 6 - 1);
            }
        }
    }
    
    // 1D
    for (int ip = 0; ip < 4; ip++) {
        inCijkl(1 - 1, 1 - 1) = mC11_1D(ip);
        inCijkl(1 - 1, 2 - 1) = mC12_1D(ip);
        inCijkl(1 - 1, 3 - 1) = mC13_1D(ip);
        inCijkl(1 - 1, 4 - 1) = mC14_1D(ip);
        inCijkl(1 - 1, 5 - 1) = mC15_1D(ip);
        inCijkl(1 - 1, 6 - 1) = mC16_1D(ip);
        
        inCijkl(2 - 1, 1 - 1) = mC12_1D(ip);
        inCijkl(2 - 1, 2 - 1) = mC22_1D(ip);
        inCijkl(2 - 1, 3 - 1) = mC23_1D(ip);
        inCijkl(2 - 1, 4 - 1) = mC24_1D(ip);
        inCijkl(2 - 1, 5 - 1) = mC25_1D(ip);
        inCijkl(2 - 1, 6 - 1) = mC26_1D(ip);
        
        inCijkl(3 - 1, 1 - 1) = mC13_1D(ip);
        inCijkl(3 - 1, 2 - 1) = mC23_1D(ip);
        inCijkl(3 - 1, 3 - 1) = mC33_1D(ip);
        inCijkl(3 - 1, 4 - 1) = mC34_1D(ip);
        inCijkl(3 - 1, 5 - 1) = mC35_1D(ip);
        inCijkl(3 - 1, 6 - 1) = mC36_1D(ip);
        
        inCijkl(4 - 1, 1 - 1) = mC14_1D(ip);
        inCijkl(4 - 1, 2 - 1) = mC24_1D(ip);
        inCijkl(4 - 1, 3 - 1) = mC34_1D(ip);
        inCijkl(4 - 1, 4 - 1) = mC44_1D(ip);
        inCijkl(4 - 1, 5 - 1) = mC45_1D(ip);
        inCijkl(4 - 1, 6 - 1) = mC46_1D(ip);
        
        inCijkl(5 - 1, 1 - 1) = mC15_1D(ip);
        inCijkl(5 - 1, 2 - 1) = mC25_1D(ip);
        inCijkl(5 - 1, 3 - 1) = mC35_1D(ip);
        inCijkl(5 - 1, 4 - 1) = mC45_1D(ip);
        inCijkl(5 - 1, 5 - 1) = mC55_1D(ip);
        inCijkl(5 - 1, 6 - 1) = mC56_1D(ip);
        
        inCijkl(6 - 1, 1 - 1) = mC16_1D(ip);
        inCijkl(6 - 1, 2 - 1) = mC26_1D(ip);
        inCijkl(6 - 1, 3 - 1) = mC36_1D(ip);
        inCijkl(6 - 1, 4 - 1) = mC46_1D(ip);
        inCijkl(6 - 1, 5 - 1) = mC56_1D(ip);
        inCijkl(6 - 1, 6 - 1) = mC66_1D(ip);
        
        // compute backazimuth
        RDCol2 rtheta = Geodesy::rtheta(mMyQuad->getNodalCoords().col(ip));
        double recLat = Geodesy::theta2Lat_r(rtheta(1), rtheta(0));
        double recDep = Geodesy::getROuter() - rtheta(0);
        double baz = Geodesy::backAzimuth(srcLat, srcLon, srcDep, recLat, 0., recDep);
        // Bond Transformation
        const RDMatXX &outCijkl = bondTransformation(inCijkl, 0., 0., -baz);
        
        mC11_1D(ip) = outCijkl(1 - 1, 1 - 1);
        mC12_1D(ip) = outCijkl(1 - 1, 2 - 1);
        mC13_1D(ip) = outCijkl(1 - 1, 3 - 1);
        mC14_1D(ip) = outCijkl(1 - 1, 4 - 1);
        mC15_1D(ip) = outCijkl(1 - 1, 5 - 1);
        mC16_1D(ip) = outCijkl(1 - 1, 6 - 1);
        
        mC22_1D(ip) = outCijkl(2 - 1, 2 - 1);
        mC23_1D(ip) = outCijkl(2 - 1, 3 - 1);
        mC24_1D(ip) = outCijkl(2 - 1, 4 - 1);
        mC25_1D(ip) = outCijkl(2 - 1, 5 - 1);
        mC26_1D(ip) = outCijkl(2 - 1, 6 - 1);
        
        mC33_1D(ip) = outCijkl(3 - 1, 3 - 1);
        mC34_1D(ip) = outCijkl(3 - 1, 4 - 1);
        mC35_1D(ip) = outCijkl(3 - 1, 5 - 1);
        mC36_1D(ip) = outCijkl(3 - 1, 6 - 1);
        
        mC44_1D(ip) = outCijkl(4 - 1, 4 - 1);
        mC45_1D(ip) = outCijkl(4 - 1, 5 - 1);
        mC46_1D(ip) = outCijkl(4 - 1, 6 - 1);
        
        mC55_1D(ip) = outCijkl(5 - 1, 5 - 1);
        mC56_1D(ip) = outCijkl(5 - 1, 6 - 1);
        
        mC66_1D(ip) = outCijkl(6 - 1, 6 - 1);
    }
}

RDMatXX Material::bondTransformation(RDMatXX inCijkl, double alpha, double beta, double gamma) {
    RDMat33 R1, R2, R3, R;
    R1 << 1., 0., 0.,
          0., cos(alpha), sin(alpha),
          0., -sin(alpha), cos(alpha);
    R2 << cos(beta), 0., -sin(beta), 
          0., 1., 0.,
          sin(beta), 0, cos(beta);
    R3 << cos(gamma), sin(gamma), 0., 
          -sin(gamma), cos(gamma), 0.,
          0., 0., 1.;
    
    R = R1 * R2 * R3;
    
    RDMat33 K1, K2, K3, K4;
    K1.array() = R.array().pow(2.);
    K2 << R(0, 1) * R(0, 2), R(0, 2) * R(0, 0), R(0, 0) * R(0, 1),
          R(1, 1) * R(1, 2), R(1, 2) * R(1, 0), R(1, 0) * R(1, 1),
          R(2, 1) * R(2, 2), R(2, 2) * R(2, 0), R(2, 0) * R(2, 1);
    K3 << R(1, 0) * R(2, 0), R(1, 1) * R(2, 1), R(1, 2) * R(2, 2),
          R(2, 0) * R(0, 0), R(2, 1) * R(0, 1), R(2, 2) * R(0, 2),
          R(0, 0) * R(1, 0), R(0, 1) * R(1, 1), R(0, 2) * R(1, 2);
    K4 << R(1, 1) * R(2, 2) + R(1, 2) * R(2, 1), 
          R(1, 2) * R(2, 0) + R(1, 0) * R(2, 2), 
          R(1, 0) * R(2, 1) + R(1, 1) * R(2, 0),
          R(2, 1) * R(0, 2) + R(2, 2) * R(0, 1), 
          R(2, 2) * R(0, 0) + R(2, 0) * R(0, 2), 
          R(2, 0) * R(0, 1) + R(2, 1) * R(0, 0),
          R(0, 1) * R(1, 2) + R(0, 2) * R(1, 1), 
          R(0, 2) * R(1, 0) + R(0, 0) * R(1, 2), 
          R(0, 0) * R(1, 1) + R(0, 1) * R(1, 0);
    
    RDMatXX K(6, 6);
    K.block(0, 0, 3, 3) = K1;
    K.block(0, 3, 3, 3) = 2. * K2;
    K.block(3, 0, 3, 3) = K3;
    K.block(3, 3, 3, 3) = K4;
    
    RDMatXX outCijkl(K * inCijkl * K.transpose());

    // deal with numerical errors
    double tol = std::abs(outCijkl(4, 4)) * 1e-7;
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            if (std::abs(outCijkl(i, j)) < tol) {
                outCijkl(i, j) = 0.;
            }
        }
    }
    return outCijkl;
}

void Material::prepare3D() {
    int Nr = mMyQuad->getNr();
    mVpv3D = RDMatXN::Zero(Nr, nPE);
    mVph3D = RDMatXN::Zero(Nr, nPE);
    mVsv3D = RDMatXN::Zero(Nr, nPE);
    mVsh3D = RDMatXN::Zero(Nr, nPE);
    mRho3D = RDMatXN::Zero(Nr, nPE);
    mEta3D = RDMatXN::Zero(Nr, nPE);
    mQkp3D = RDMatXN::Zero(Nr, nPE);
    mQmu3D = RDMatXN::Zero(Nr, nPE);
    
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
            mEta3D.col(ipnt).fill(Mapping::interpolate(mEta1D, xieta));
            mQkp3D.col(ipnt).fill(Mapping::interpolate(mQkp1D, xieta));
            mQmu3D.col(ipnt).fill(Mapping::interpolate(mQmu1D, xieta));
            // rho for mass
            int NrP = mMyQuad->getPointNr(ipol, jpol);
            mRhoMass3D[ipnt] = RDColX::Constant(NrP, mRho3D.col(ipnt)(0));
            mVpFluid3D[ipnt] = RDColX::Constant(NrP, mVpv3D.col(ipnt)(0));
        }
    }
}

bool Material::_3Dprepared() const {
    return mVpv3D.rows() == mMyQuad->getNr();
}


