// Volumetric3D.h
// created by Kuangdai on 16-May-2016 
// base class of volumetric 3D models, including 
// velocity, density, attenuation

#pragma once
#include <string>
#include <vector>

class Parameters;
class ExodusModel;

class Volumetric3D {
public:
    // reference type of material properties
    enum MaterialRefType {Absolute, Reference1D, Reference3D, ReferencePerturb};
    const std::vector<std::string> MaterialRefTypeStringFull = {"Absolute", 
        "Reference1D", "Reference3D", "ReferencePerturb"};
    const std::vector<std::string> MaterialRefTypeStringShort = {"Abs", 
        "Ref1D", "Ref3D", "RefPtb"};
        
    // allowed keys for material properties
    enum MaterialProperty {VP, VPV, VPH, VS, VSV, VSH, RHO, ANIS_ETA, QMU, QKAPPA};
    const std::vector<std::string> MaterialPropertyString = {"VP", "VPV", "VPH", 
        "VS", "VSV", "VSH", "RHO", "ANIS_ETA", "QMU", "QKAPPA"};
    
    virtual ~Volumetric3D() {finalize();};
    
    // initialize internal variables if needed
    virtual void initialize() {};
    virtual void initialize(const std::vector<std::string> &params) {initialize();};
    
    // finalize internal variables if needed
    virtual void finalize() {};

    // get perturbations or absolute values at location r/theta/phi 
    // IMPORTANT NOTES: 
    // a) This function should be realized such that r/theta/phi are the geocentric 
    //    coordinates, without rotating the source to the north pole.  
    //    For models given in geographic coordinates, geocentric-to-geographic  
    //    conversions from (theta, phi) to (lat, lon) has to be performed internally.
    // b) rElemCenter: radius of the element center, used to decide whether the given 
    //    location is within the model range. We use rElemCenter instead of r itself 
    //    to honour vertical discontinuities.
    // c) output: propNames, propValues, propRefTypes
    // d) return: a "false" return means the outputs are all ignored, 
    //    e.g., the input location is out of the model range
    virtual bool get3dProperties(double r, double theta, double phi, double rElemCenter,
        std::vector<MaterialProperty> &propNames, 
        std::vector<double> &propValues, 
        std::vector<MaterialRefType> &propRefTypes) const {return false;};
        
    // get perturbations or absolute values at location r/theta/phi 
    // IMPORTANT: this function should be realized such that r/theta/phi are 
    // the geocentric coordinates in the SPHERICAL reference configuration, that is,
    // the perfect sphere before topography and ellipticity
    // RETURN: a "false" return means the input location is out of model range
    // rElemCenter: radius of the element center, used to decide whether the given 
    // location is within model range. We use rElemCenter instead of r itself 
    // to honour vertical discontinuities
    virtual bool get3dProperties(double r, double theta, double phi, double rElemCenter,
        double &vpv, double &vph, double &vsv, double &vsh, double &rho) const {return false;};
    
    // reference type
    virtual MaterialRefType getReferenceType() const {return MaterialRefType::Absolute;};
    
    // verbose 
    virtual std::string verbose() const = 0;
    
    // set source location, if needed
    virtual void setSourceLocation(double srcLat, double srcLon, double srcDep) {};
    
    // obtain additional mesh information from ExodusModel, if needed
    virtual void setupExodusModel(const ExodusModel *exModel) {};
    
    // build from input parameters
    static void buildInparam(std::vector<Volumetric3D *> &models, 
        const Parameters &par, const ExodusModel *exModel, 
        double srcLat, double srcLon, double srcDep, int verbose);

};
