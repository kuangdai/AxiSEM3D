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
    const std::vector<std::string> MaterialRefTypeString = {"Absolute", 
        "Reference1D", "Reference3D", "ReferencePerturb"};
    const std::vector<std::string> MaterialRefTypeStringShort = {"Abs", 
        "Ref1D", "Ref3D", "RefPtb"};
        
    // allowed keys for material properties
    enum MaterialProperty {VPV, VPH, VSV, VSH, RHO, ANIS_ETA, QKAPPA, QMU, VP, VS, 
        C11, C12, C13, C14, C15, C16,
        C22, C23, C24, C25, C26,
        C33, C34, C35, C36,
        C44, C45, C46,
        C55, C56,
        C66
    };
    const std::vector<std::string> MaterialPropertyString = {"VPV", "VPH", "VSV", "VSH", 
        "RHO", "ANIS_ETA", "QKAPPA", "QMU", "VP", "VS",
        "C11", "C12", "C13", "C14", "C15", "C16",
        "C22", "C23", "C24", "C25", "C26",
        "C33", "C34", "C35", "C36",
        "C44", "C45", "C46",
        "C55", "C56",
        "C66"
    };
    const std::vector<double> MaterialPropertyAbsSI = {1e3, 1e3, 1e3, 1e3, 
        1e3, 1., 1., 1., 1e3, 1e3, 
        // Cijkl should be given in GPa
        1e9, 1e9, 1e9, 1e9, 1e9, 1e9, 
        1e9, 1e9, 1e9, 1e9, 1e9, 
        1e9, 1e9, 1e9, 1e9,
        1e9, 1e9, 1e9,
        1e9, 1e9, 
        1e9
    };
    
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
    // c) output: properties, refTypes, values
    // d) return: a "false" return means the outputs are all ignored, 
    //    e.g., the input location is out of the model range
    virtual bool get3dProperties(double r, double theta, double phi, double rElemCenter,
        std::vector<MaterialProperty> &properties, 
        std::vector<MaterialRefType> &refTypes,
        std::vector<double> &values) const = 0;
        
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
        
    // some models may extent into the fluid core
    // but may not create 3D fluid actually
    virtual bool makeFluid3D() const {return false;};    

};
