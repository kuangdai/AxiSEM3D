// Geometric3D.h
// created by Kuangdai on 4-Jun-2016 
// base class of geometric 3D models, such as 
// topography, ellipticity and undulating Moho

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
    // IMPORTANT NOTES: 
    // a) This function should be realized such that r/theta/phi are the geocentric 
    //    coordinates, without rotating the source to the north pole.  
    //    For models given in geographic coordinates, geocentric-to-geographic  
    //    conversions from (theta, phi) to (lat, lon) has to be performed internally.
    // b) All models should be defined independently with respect to the perfect sphere.
    virtual double getDeltaR(double r, double theta, double phi, double rElemCenter) const = 0;
    
    // verbose 
    virtual std::string verbose() const = 0;
    
    // build from input parameters
    static void buildInparam(std::vector<Geometric3D *> &models, 
        const Parameters &par, int verbose);

};

