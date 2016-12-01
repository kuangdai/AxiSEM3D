// Geometric3D.h
// created by Kuangdai on 4-Jun-2016 
// base class of geometric 3D models, including topography and ellipticity

#pragma once
#include <string>
#include <vector>

class Parameters;

class Geometric3D {
public:

    virtual ~Geometric3D() {finalize();};
    
    // initialize internal variables if needed
    virtual void initialize() {};
    virtual void initialize(const std::vector<std::string> &params) {initialize();};
    
    // finalize internal variables if needed
    virtual void finalize() {};
    
    // get undulation (deltaR)
    // IMPORTANT: this function should be realized such that r/theta/phi are 
    // the geocentric coordinates in the SPHERICAL reference configuration, that is,
    // the perfect sphere before topography and ellipticity
    // For topography models given in geographic coordinates, one needs to perform
    // the conversions internally
    virtual double getDeltaR(double r, double theta, double phi, double rElemCenter) const = 0;
    
    // // get gradient of deltaR
    // // deltaR_r     = d(deltaR)/d(r)
    // // deltaR_theta = d(deltaR)/d(theta) / r
    // // deltaR_phi   = d(deltaR)/d(phi) / (r * sin(theta))
    // // though the general method based on variation is provided in Geometric3D.cpp, we 
    // // strongly recommend explicit implementation of getNablaDeltaR() for each model
    // virtual bool getNablaDeltaR(double r, double theta, double phi, double rElemCenter,
    //     double &deltaR_r, double &deltaR_theta, double &deltaR_phi) const;
    
    // verbose 
    virtual std::string verbose() const = 0;
    
    // set outer radius
    virtual void setROuter(double router) {};
    
    // build from input parameters
    static void buildInparam(std::vector<Geometric3D *> &models, 
        const Parameters &par, int verbose);

};

