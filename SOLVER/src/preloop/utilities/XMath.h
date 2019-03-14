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
        
    // linear interpolation
    static void interpLinear(double target, const RDColX &bases, int &loc, double &weight); 
    
    // check sorted
    template<class TIN>
    static bool sortedAscending(const TIN &bases) {
        for (int i = 0; i < bases.size() - 1; i++) {
            if (bases(i) > bases(i + 1)) {
                return false;
            }
        }
        return true;
    }
    
    // check limits
    static void checkLimits(double &value, double low, double up, double tol = tinyDouble);
    
    // Gaussian smoothing
    static void gaussianSmoothing(RDColX &data, int order, double dev, bool period);
    static void gaussianSmoothing(RDMatXX &data, 
        IColX orderRow, RDColX devRow, bool periodRow, 
        IColX orderCol, RDColX devCol, bool periodCol);
    
    ////////////////// flatten-structured cast //////////////////
    // structured
    template<class TIN>
    static void structuredUseFirstRow(const TIN &flat, RDMatPP &strct) {
        for (int ipol = 0; ipol < nPntEdge; ipol++) {
            strct.block(ipol, 0, 1, nPntEdge) 
            = flat.block(0, nPntEdge * ipol, 1, nPntEdge);
        }
    };
    
    // flatten
    template<class TOUT>
    static void flattenFillWithFirstRow(const RDMatPP &strct, TOUT &flat) {
        int nr = flat.rows();
        for (int ipol = 0; ipol < nPntEdge; ipol++) {
            flat.block(0, nPntEdge * ipol, nr, nPntEdge) 
            = strct.block(ipol, 0, 1, nPntEdge).colwise().replicate(nr); 
        }
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
    
    
    // memory info for eigen
    template<class EigenMat>
    static std::string eigenMemoryInfo(const std::string &title, const EigenMat &mat) {
        // scalar size
        size_t scalar = sizeof(typename EigenMat::Scalar);
        // total size
        double memGB = (mat.size() * scalar) / 1e9;
        // format
        std::stringstream ss;
        ss << title << ": ";
        ss << "dimensions = " << mat.rows() << " x " << mat.cols() << "; ";
        ss << "scalar bytes = " << scalar << "; ";
        ss << "memory = " << memGB << " GB";
        return ss.str();
    }


};

