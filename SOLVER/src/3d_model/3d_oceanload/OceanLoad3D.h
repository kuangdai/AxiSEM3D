// OceanLoad3D.h
// created by Kuangdai on 8-Oct-2016 
// base class of Ocean Load models

#pragma once
#include <string>
#include <vector>

class Parameters;

class OceanLoad3D {
public:

    virtual ~OceanLoad3D() {finalize();};
    
    // initialize internal variables if needed
    virtual void initialize() {};
    virtual void initialize(const std::vector<std::string> &params) {initialize();};
    
    // finalize internal variables if needed
    virtual void finalize() {};
    
    // get water depth at location theta/phi 
    virtual double getOceanDepth(double theta, double phi) const = 0;
    
    // verbose 
    virtual std::string verbose() const = 0;
    
    // build from input parameters
    static void buildInparam(OceanLoad3D *&model, 
        const Parameters &par, int verbose);

    static const double mWaterDensity;
};
