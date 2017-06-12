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
        mQkp1D(i) = exModel.getElementalVariables("QKAPPA_" + std::to_string(i), quadTag);
        mQmu1D(i) = exModel.getElementalVariables("QMU_" + std::to_string(i), quadTag);
    }
    
    // initialize 3D properties with 1D reference
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
            mRhoMass3D[ipnt] = RDColX::Constant(NrP, mRho3D.col(ipnt)[0]);
            mVpFluid3D[ipnt] = RDColX::Constant(NrP, mVpv3D.col(ipnt)[0]);
        }
    }
}

void Material::addVolumetric3D(const std::vector<Volumetric3D *> &m3D, 
    double srcLat, double srcLon, double srcDep, double phi2D) {
    if (m3D.size() == 0) {
        return;
    }    
        
    // pointers for fast access to material matrices
    std::vector<RDRow4 *>  prop1DPtr = {&mVpv1D, &mVph1D, &mVsv1D, &mVsh1D, &mRho1D, &mEta1D, &mQkp1D, &mQmu1D};
    std::vector<RDMatXN *> prop3DPtr = {&mVpv3D, &mVph3D, &mVsv3D, &mVsh3D, &mRho3D, &mEta3D, &mQkp3D, &mQmu3D};
    
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
                    std::vector<Volumetric3D::MaterialProperty> properties; 
                    std::vector<Volumetric3D::MaterialRefType> refTypes;
                    std::vector<double> values;
                    if (!model->get3dProperties(r, t, p, rElemCenter, properties, refTypes, values)) {
                        // point (r, t, p) not in model range
                        continue;
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
    for (int ipol = 0; ipol <= nPol; ipol++) {
        for (int jpol = 0; jpol <= nPol; jpol++) {
            int ipnt = ipol * nPntEdge + jpol;
            int nr_mass = mMyQuad->getPointNr(ipol, jpol);
            mRhoMass3D[ipnt] = XMath::linearResampling(nr_mass, mRho3D.col(ipnt));
            mVpFluid3D[ipnt] = XMath::linearResampling(nr_mass, mVpv3D.col(ipnt));
        }
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
    RDMatXN fluidK = mRho3D.array().pow(-1.);
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
    // elasticity tensor
    const RDRowN &iFact = mMyQuad->getIntegralFactor(); 
    RDMatXN A(mRho3D), C(mRho3D), F(mRho3D), L(mRho3D), N(mRho3D);
    for (int ipnt = 0; ipnt < nPntElem; ipnt++) {
        // A C L N
        A.col(ipnt).array() *= mVph3D.col(ipnt).array().pow(2.) * iFact(ipnt);
        C.col(ipnt).array() *= mVpv3D.col(ipnt).array().pow(2.) * iFact(ipnt);
        L.col(ipnt).array() *= mVsv3D.col(ipnt).array().pow(2.) * iFact(ipnt);
        N.col(ipnt).array() *= mVsh3D.col(ipnt).array().pow(2.) * iFact(ipnt);
    }
    // F
    F = mEta3D.schur(A - 2. * L);
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
            att1D = attBuild->createAttenuation1D(mQkp3D, mQmu3D, kappa, mu, mMyQuad);
        } else {
            att3D = attBuild->createAttenuation3D(mQkp3D, mQmu3D, kappa, mu, mMyQuad);
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

double Material::getVMaxRef() const {
    return std::max(mVph1D.maxCoeff(), mVpv1D.maxCoeff());
}

RDColX Material::getVMax() const {
    const RDColX &vpvMax = mVpv3D.rowwise().maxCoeff();
    const RDColX &vphMax = mVph3D.rowwise().maxCoeff();
    return (vpvMax.array().max(vphMax.array())).matrix();
}

bool Material::isFluidPar1D() const {
    return XMath::equalRows(mVpv3D) && XMath::equalRows(mRho3D);
}

bool Material::isSolidPar1D(bool attenuation) const {
    bool result = XMath::equalRows(mVpv3D) && XMath::equalRows(mVph3D) &&
                  XMath::equalRows(mVsv3D) && XMath::equalRows(mVsh3D) && 
                  XMath::equalRows(mRho3D) && XMath::equalRows(mEta3D);
    if (attenuation) {
        result = result && XMath::equalRows(mQkp3D) && XMath::equalRows(mQmu3D);
    }              
    return result;
}

bool Material::isIsotropic() const {
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
    } else {
        throw std::runtime_error("Material::getProperty || Unknown field variable name: " + vname);
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
