// XMath.h
// created by Kuangdai on 9-May-2016 
// small math tools

#include "eigenp.h"
#include "eigenc.h"

#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <typeinfo>

class XMath {
public:
    ////////////////// coordinate related //////////////////
    // 2D, (s, z) -> (r, theta)
    static void rtheta(const RDCol2 &sz, double &r, double &theta);
    static RDCol2 rtheta(const RDCol2 &sz);
    static double theta(const RDCol2 &sz) {return rtheta(sz)(1);};
    // phi = atan(y / x), phi in [0, 2 pi)
    static double atan4(double y, double x, bool &defined);
    // rotation matrix   
    static RDMat33 rotationMatrix(double theta, double phi);
    // spherical vector basis
    static RDMat33 sphericalBasis(double theta, double phi);
    // spherical to Cartesian
    static RDCol3 toCartesian(const RDCol3 &rtp);
    // Cartesian to spherical
    static RDCol3 toSpherical(const RDCol3 &xyz, bool &defined);
    // geographic to geocentric 
    static double lat2Theta(double lat, double depth);
    static double lon2Phi(double lon);
    // geocentric to geographic 
    static double theta2Lat(double theta, double depth);
    static double phi2Lon(double lon);
    // source-centered to globe
    static RDCol3 rotateSrc2Glob(const RDCol3 &rtpS, double srclat, double srclon, double srcdep);
    // globe to source-centered
    static RDCol3 rotateGlob2Src(const RDCol3 &rtpG, double srclat, double srclon, double srcdep);
    // compute back azimuth (copied from specfem)
    static double backAzimuth(double srclat, double srclon, double srcdep,
                              double reclat, double reclon, double recdep);
    // a and b are two angles close to each other  
    // but their values may be very different, e.g., a = 0.01, b = 2 * pi - 0.01,
    // which is problematic when linear interpolation is carried out
    static void makeClose(double &a, double &b);     
    // find closest distance among 2D points
    static double findClosestDist(const std::vector<RDCol2> &crds);         
    
    // Lagrange interpolation
    static void interpLagrange(double target, int nbases, const double *bases, double *results);
    
    // Gaussian smoothing
    static void gaussianSmoothing(RDColX &data, int order, double dev, bool period);
    static void gaussianSmoothing(RDMatXX &data, 
        IColX orderRow, RDColX devRow, bool periodRow, 
        IColX orderCol, RDColX devCol, bool periodCol);
    
    ////////////////// fftw related //////////////////
    // http://www.fftw.org/fftw2_doc/fftw_3.html
    // FFTW is best at handling sizes of the form 2^a 3^b 5^c 7^d 11^e 13^f, 
    // where e+f is either 0 or 1, and the other exponents are arbitrary.
    // We call numbers of the form lucky numbers.
    static bool isLuckyNumber(int n, bool forceOdd = false);
    static int nextLuckyNumber(int n, bool forceOdd = false);
    
    
    ////////////////// preproccessor-to-solver cast //////////////////
    // structured
    static RMatPP castToSolver(const RDMatPP &mp);
    
    // flattened
    template<class TIN, class TOUT>
    static TOUT castToSolver(const TIN &mat) {
        return mat.template cast<Real>();
    };
    
    ////////////////// Fourier methods //////////////////
    // check if a matrix has identical rows
    template<class TMat>
    static bool equalRows(const TMat &mat, double tol = tinyDouble) {
        bool equal = true;
        for (int i = 1; i < mat.rows(); i++) {
            equal = equal && (mat.row(i) - mat.row(0)).norm() < tol * mat.row(0).norm();
            if (!equal) break;
        }
        return equal;
    };
    
    // resampling
    static RDColX trigonResampling(int newSize, const RDColX &original);
    static RDColX linearResampling(int newSize, const RDColX &original);
    
    // Fourier
    static RDRowN computeFourierAtPhi(const RDMatXN &data, double phi);
    
    /////////////// string cast ////////////////
    template<typename parType>
    static void castValue(parType &result, const std::string &val_in, const std::string &source) {
        try {
            // special care of bool
            std::string val = val_in;
            if (typeid(parType) == typeid(bool)) {
                boost::to_upper<std::string>(val);
                if (val == "TRUE" || val == "YES" || val == "ON") val = "1";
                if (val == "FALSE" || val == "NO" || val == "OFF") val = "0";
            }
            result = boost::lexical_cast<parType>(val);
        } catch(std::exception) {
            throw std::runtime_error("XMath::castValue || "
                "Invalid argument encountered in " + source + ", arg = " + val_in + ".");
        }
    };
    
    //////////////// Ellipticity ////////////////
    static void setEllipticity(double flattening, double router, 
        const std::vector<double> &ellip_knots, 
        const std::vector<double> &ellip_coeffs);
    static double getFlattening(double r);
    static double getFlattening() {return sFlattening;};
    static double getROuter() {return sROuter;};
private:
    static double sFlattening;
    static double sROuter;
    static RDColX sEllipKnots;
    static RDColX sEllipCoeffs;
};

