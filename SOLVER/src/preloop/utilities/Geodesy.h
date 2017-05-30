// Geodesy.h
// created by Kuangdai on 15-May-2017 
// geodetic tools

#include "eigenp.h"

class Geodesy {
public:
    
    // 2D, (s, z) -> (r, theta)
    static void rtheta(const RDCol2 &sz, double &r, double &theta);
    static RDCol2 rtheta(const RDCol2 &sz);
    static double theta(const RDCol2 &sz) {return rtheta(sz)(1);};
    
    // phi = atan(y / x), phi in [0, 2 pi)
    static double atan4(double y, double x, bool &defined);
    // spherical to Cartesian
    static RDCol3 toCartesian(const RDCol3 &rtp);
    // Cartesian to spherical
    static RDCol3 toSpherical(const RDCol3 &xyz, bool &defined);
    // rotation matrix   
    static RDMat33 rotationMatrix(double theta, double phi);
        
    // geographic to geocentric 
    static double lat2Theta_r(double lat, double radius);
    static double lat2Theta_d(double lat, double depth) {
        return lat2Theta_r(lat, sROuter - depth);
    };
    static double lon2Phi(double lon);
    
    // geocentric to geographic 
    static double theta2Lat_r(double theta, double radius);
    static double theta2Lat_d(double theta, double depth) {
        return theta2Lat_r(theta, sROuter - depth);
    };
    static double phi2Lon(double lon);
    
    // source-centered to globe
    static RDCol3 rotateSrc2Glob(const RDCol3 &rtpS, double srclat, double srclon, double srcdep);
    
    // globe to source-centered
    static RDCol3 rotateGlob2Src(const RDCol3 &rtpG, double srclat, double srclon, double srcdep);
    
    // compute back azimuth (copied from specfem)
    static double backAzimuth(double srclat, double srclon, double srcdep,
                              double reclat, double reclon, double recdep);
                              
    
    // setup
    static void setup(double router, double flattening, 
        const RDColX &ellip_knots, 
        const RDColX &ellip_coeffs);
        
    // compute flattening at a radius
    static double getFlattening(double r);
    
    // get
    static double getROuter() {return sROuter;};
    static double getFlattening() {return sFlattening;};
    
private:
    static double sROuter;
    static double sFlattening;
    static RDColX sEllipKnots;
    static RDColX sEllipCoeffs;
};



