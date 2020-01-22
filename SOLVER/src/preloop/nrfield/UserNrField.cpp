// UserNrField.cpp
// created by Kuangdai on 11-Jun-2017 
// user-defined nr integer field

#include "UserNrField.h"
#include <sstream>
#include "Geodesy.h"

UserNrField::UserNrField(bool useLucky, const std::vector<double> &params): 
NrField(useLucky), mParameters(params) {
    // nothing
}

int UserNrField::getNrAtPoint(const RDCol2 &coords) const {
    // input: coordinates
    // s, z, r, theta, depth
    double s = coords(0);
    double z = coords(1);
    double r, theta;
    Geodesy::rtheta(coords, r, theta);
    double depth = Geodesy::getROuter() - r;
    
    // output: Fourier Expansion order nu
    double nu = 2.;
    
    ////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////
    // TODO: Compute nu from s, z, r, theta, depth and mParameters
    // NOTE: Only edit within this box!
    //       If this is the first time you are looking at this part
    //       of code, the following example implements a Nu field
    //       that depends on depth only. The depth-nu profile is given
    //       by NU_USER_PARAMETER_LIST as (depth in meters): 
    //       NU_USER_PARAMETER_LIST   depth0 nu0 depth1 nu1 depth2 nu2 ...
    
    // check size
    int ndepth = mParameters.size() / 2;
    if (ndepth == 0 || mParameters.size() % 2 == 1) {
        throw std::runtime_error("UserNrField::getNrAtPoint || "
        "Using the default UserNrField.cpp, "
        "the number of the user-defined parameters must be even and positive.");
    }
    
    if (depth <= mParameters[0]) {
        // shallower than shallowest
        nu = mParameters[1];
    } else if (depth >= mParameters[(ndepth - 1) * 2]) {
        // deeper than deepest
        nu = mParameters[(ndepth - 1) * 2 + 1];
    } else {
        // middle
        for (int i = 1; i < ndepth; i++) {
            double depth0 = mParameters[i * 2 - 2];
            double depth1 = mParameters[i * 2];
            double nu0 = mParameters[i * 2 - 1];
            double nu1 = mParameters[i * 2 + 1];
            if (depth <= depth1) {
                nu = (nu1 - nu0) / (depth1 - depth0) * (depth - depth0) + nu0;
                break;
            }
        }
    }
    
    // NOTE: Only edit within this box!
    ////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////
    
    // get nr from nu
    int nr = 2 * nu + 1;
    if (nr <= 0) {
        throw std::runtime_error("UserNrField::getNrAtPoint || Non-positive Nr.");
    }
    return nr;
}

std::string UserNrField::verbose() const {
    std::stringstream ss;
    ss << "\n================= Fourier Expansion Order ==================" << std::endl;
    ss << "  Type                     =   User-defined" << std::endl;
    if (mParameters.size() > 0) {
        ss << "  Parameter List           =   ";
        for (int ipar = 0; ipar < mParameters.size(); ipar++) {
            ss << mParameters[ipar] << " ";
        }
        ss << std::endl;
    }
    ss << "  Use FFTW Lucky Numbers   =   " << (mUseLuckyNumber ? "YES" : "NO") << std::endl;
    ss << "================= Fourier Expansion Order ==================\n" << std::endl;
    return ss.str();
}

