// XMath.h
// created by Kuangdai on 9-May-2016 
// miscellaneous math tools

#include "eigenp.h"
#include "eigenc.h"

class XMath {
public:
    // a and b are two angles close to each other  
    // but their values may be very different, e.g., a = 0.01, b = 2 * pi - 0.01,
    // which is problematic when linear interpolation is carried out
    static void makeClose(double &a, double &b);     
    
    // find closest distance among 2D points
    static double findClosestDist(const std::vector<RDCol2> &crds);         
    
    // Lagrange interpolation
    static void interpLagrange(double target, int nbases, 
        const double *bases, double *results);
    
    // Gaussian smoothing
    static void gaussianSmoothing(RDColX &data, int order, double dev, bool period);
    static void gaussianSmoothing(RDMatXX &data, 
        IColX orderRow, RDColX devRow, bool periodRow, 
        IColX orderCol, RDColX devCol, bool periodCol);
    
    ////////////////// preproccessor-to-solver cast //////////////////
    // structured
    static RMatPP castToSolver(const RDMatPP &mp);
    
    // flattened
    template<class TIN, class TOUT>
    static TOUT castToSolver(const TIN &mat) {
        return mat.template cast<Real>();
    };
    
    ////////////////// Fourier-related methods //////////////////
    // check if a matrix has identical rows
    template<class TMat>
    static bool equalRows(const TMat &mat, double tol = tinyDouble) {
        for (int i = 1; i < mat.rows(); i++) {
            bool equal = (mat.row(i) - mat.row(0)).norm() <= tol * mat.row(0).norm();
            if (!equal) return false;
        }
        return true;
    };
    
    // resampling
    static RDColX trigonResampling(int newSize, const RDColX &original);
    static RDColX linearResampling(int newSize, const RDColX &original);
    
    // compute value at phi 
    static RDRowN computeFourierAtPhi(const RDMatXN &data, double phi);
};

