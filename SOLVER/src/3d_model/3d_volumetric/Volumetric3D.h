// Volumetric3D.h
// created by Kuangdai on 16-May-2016 
// base class of Volumetric 3D models, including both mantle and crust
// currently we assume 1D Q-factors (Q_mu and Q_kappa)

#pragma once
#include <string>
#include <vector>

class Parameters;

class Volumetric3D {
public:
    // reference type
    enum ReferenceTypes {Absolute, Reference1D, Reference3D, ReferencePerturb};
    const std::string ReferenceTypesString[4] 
        {"Absolute", "Reference1D", "Reference3D", "ReferencePerturb"};
    
    virtual ~Volumetric3D() {finalize();};
    
    // initialize internal variables if needed
    virtual void initialize() {};
    virtual void initialize(const std::vector<std::string> &params) {initialize();};
    
    // finalize internal variables if needed
    virtual void finalize() {};
    
    // get perturbations or absolute values at location r/theta/phi 
    // IMPORTANT: this function should be realized such that r/theta/phi are 
    // the geocentric coordinates in the SPHERICAL reference configuration, that is,
    // the perfect sphere before topography and ellipticity
    // RETURN: a "false" return means the input location is out of model range
    // rElemCenter: radius of the element center, used to decide whether the given 
    // location is within model range. We use rElemCenter instead of r itself 
    // to honour vertical discontinuities
    virtual bool get3dProperties(double r, double theta, double phi, double rElemCenter,
        double &vpv, double &vph, double &vsv, double &vsh, double &rho) const = 0;
    
    // reference type
    virtual ReferenceTypes getReferenceType() const = 0;
    
    // verbose 
    virtual std::string verbose() const = 0;
    
    // set outer radius
    virtual void setROuter(double router) {};
    
    // build from input parameters
    static void buildInparam(std::vector<Volumetric3D *> &models, 
        const Parameters &par, int verbose);

};
