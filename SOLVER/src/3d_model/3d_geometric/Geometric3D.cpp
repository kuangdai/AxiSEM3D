// Geometric3D.cpp
// created by Kuangdai on 4-Jun-2016 
// base class of geometric 3D models, including topography and ellipticity

#include "Geometric3D.h"
#include "global.h"

#include "XMPI.h"
#include "XMath.h"
#include "Parameters.h"
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

/////////////////////////////// user-defined models here
#include "Ellipticity.h"
#include "Geometric3D_crust1.h"
/////////////////////////////// user-defined models here

void Geometric3D::buildInparam(std::vector<Geometric3D *> &models, 
    const Parameters &par, int verbose) {
    for (const auto &m: models) delete m;    
    models.clear();
    // first check size
    int nmodels = par.getValue<int>("MODEL_3D_GEOMETRIC_NUM");
    int nsize = par.getSize("MODEL_3D_GEOMETRIC_LIST");
    if (nmodels > nsize) throw std::runtime_error("Geometric3D::buildInparam || "
        "Not enough model names provided in MODEL_3D_GEOMETRIC_LIST ||"
        "MODEL_3D_GEOMETRIC_NUM = " + par.getValue<std::string>("MODEL_3D_GEOMETRIC_NUM") + ".");
    
    for (int i = 0; i < nmodels; i++) {
        // split model name and parameters
        std::string mstr = par.getValue<std::string>("MODEL_3D_GEOMETRIC_LIST", i);
        std::vector<std::string> strs;
        boost::trim_if(mstr, boost::is_any_of("\t "));
        boost::split(strs, mstr, boost::is_any_of("$"), boost::token_compress_on);
        std::string name = strs[0];
        std::vector<double> params;
        try { 
            for (int i = 1; i < strs.size(); i++) 
                params.push_back(boost::lexical_cast<double>(strs[i]));
        } catch (std::exception) {
            throw std::runtime_error("Geometric3D::buildInparam || "
                "Invalid parameter following geometric model " + name + ".");
        }
        
        // create model
        Geometric3D * m;
        if (boost::iequals(name, "crust1")) {
            m = new Geometric3D_crust1();
            
        /////////////////////////////// 
        // user-defined models here
        /////////////////////////////// 
            
        } else {
            throw std::runtime_error("Geometric3D::buildInparam || "
                "Unknown geometric model name " + name + ".");
        }
        
        // initialize
        m->setROuter(XMath::getROuter());
        m->initialize(params);
        if (verbose) XMPI::cout << m->verbose();
        models.push_back(m);
    }
    
    // ELLIPTICITY
    std::string emode = par.getValue<std::string>("MODEL_3D_ELLIPTICITY_MODE");
    if (boost::iequals(emode, "full")) {
        Geometric3D *m = new Ellipticity();
        if (verbose) XMPI::cout << m->verbose();
        models.push_back(m);
    }
}

// bool Geometric3D::getNablaDeltaR(double r, double theta, double phi, double rElemCenter,
//     double &deltaR_r, double &deltaR_theta, double &deltaR_phi) const {
//     // original value
//     double deltaR, deltaR_p, deltaR_m;
//     bool within = getDeltaR(r, theta, phi, rElemCenter, deltaR);
//     if (!within) {
//         deltaR_r = deltaR_theta = deltaR_phi = 0.;
//         return false;
//     }
//     
//     // avoid doing finite difference at origin
//     if (r < tinyDouble) throw std::runtime_error("Geometric3D::getNablaDeltaR || Undefined gradient. " 
//         "Check model or implement getNablaDeltaR() explicitly.");
//             
//     // finite difference in r-direction
//     double dr = 100.;
//     if (rElemCenter >= r) {
//         double r_p = r + dr;
//         bool within_p = getDeltaR(r_p, theta, phi, rElemCenter, deltaR_p);
//         if (within_p) deltaR_r = (deltaR_p - deltaR) / dr;
//         else throw std::runtime_error("Geometric3D::getNablaDeltaR || Undefined gradient. " 
//             "Check model or implement getNablaDeltaR() explicitly.");
//     } else {
//         double r_m = r - dr;
//         bool within_m = getDeltaR(r_m, theta, phi, rElemCenter, deltaR_m);
//         if (within_m) deltaR_r = (deltaR - deltaR_m) / dr;
//         else throw std::runtime_error("Geometric3D::getNablaDeltaR || Undefined gradient. " 
//             "Check model or implement getNablaDeltaR() explicitly.");
//     }
//     
//     double dtheta = .001 * degree;
//     double dphi = .001 * degree;
//     if (theta < dtheta || theta > pi - dtheta) {
//         // on the axis, deltaR_theta and deltaR_phi are not defined (as phi is indeterminant)
//         // so we compute Cartesian components deltaR_x and deltaR_y 
//         // and rotate the result to the orientation of phi
//         bool north = theta < dtheta;
//         double theta_pole = north ? dtheta : pi - dtheta;
//         double drdx = 0., drdy = 0.;
//         for (int iphi = 0; iphi < 360; iphi++) {
//             double phi_ring = iphi * degree;
//             double drdt = 0., drdp = 0.;
//             
//             // finite difference in theta-direction
//             double theta_p = theta_pole + dtheta;
//             double theta_m = theta_pole - dtheta;
//             bool within_p = getDeltaR(r, theta_p, phi_ring, rElemCenter, deltaR_p);
//             bool within_m = getDeltaR(r, theta_m, phi_ring, rElemCenter, deltaR_m);
//             if (within_p && within_m) drdt = (deltaR_p - deltaR_m) / (2. * dtheta);
//             else throw std::runtime_error("Geometric3D::getNablaDeltaR || Undefined gradient. " 
//                 "Check model or implement getNablaDeltaR() explicitly.");
//             drdt /= r;
//             
//             // finite difference in phi-direction
//             double sint = sin(theta_pole);
//             double phi_p = phi_ring + dphi;
//             double phi_m = phi_ring - dphi;
//             if (phi_p > 2. * pi) phi_p -= 2. * pi;
//             if (phi_m < 0.) phi_m += 2. * pi;
//             within_p = getDeltaR(r, theta_pole, phi_p, rElemCenter, deltaR_p);
//             within_m = getDeltaR(r, theta_pole, phi_m, rElemCenter, deltaR_m);
//             if (within_p && within_m) drdp = (deltaR_p - deltaR_m) / (2. * dphi);
//             else throw std::runtime_error("Geometric3D::getNablaDeltaR || Undefined gradient. " 
//                 "Check model or implement getNablaDeltaR() explicitly.");
//             drdp /= r * sint;
//             
//             // convert to Cartesian
//             // the relation between theta-axis and x-/y-axis differs for north and south
//             if (!north) drdt *= -1.;
//             drdx += drdt * cos(phi_ring) - drdp * sin(phi_ring);
//             drdy += drdt * sin(phi_ring) + drdp * cos(phi_ring);
//         }
//         drdx /= 360.;
//         drdy /= 360.;
//         
//         // convert back to spherical
//         deltaR_theta = drdx * cos(phi) + drdy * sin(phi);
//         deltaR_phi = -drdx * sin(phi) + drdy * cos(phi);
//         if (!north) deltaR_theta *= -1.;
//         
//     } else {
//         // finite difference in theta-direction
//         double theta_p = theta + dtheta;
//         double theta_m = theta - dtheta;
//         bool within_p = getDeltaR(r, theta_p, phi, rElemCenter, deltaR_p);
//         bool within_m = getDeltaR(r, theta_m, phi, rElemCenter, deltaR_m);
//         if (within_p && within_m) deltaR_theta = (deltaR_p - deltaR_m) / (2. * dtheta);
//         else throw std::runtime_error("Geometric3D::getNablaDeltaR || Undefined gradient. " 
//             "Check model or implement getNablaDeltaR() explicitly.");
//         deltaR_theta /= r;
//         
//         // finite difference in phi-direction
//         double sint = sin(theta);
//         double phi_p = phi + dphi;
//         double phi_m = phi - dphi;
//         if (phi_p > 2. * pi) phi_p -= 2. * pi;
//         if (phi_m < 0.) phi_m += 2. * pi;
//         within_p = getDeltaR(r, theta, phi_p, rElemCenter, deltaR_p);
//         within_m = getDeltaR(r, theta, phi_m, rElemCenter, deltaR_m);
//         if (within_p && within_m) deltaR_phi = (deltaR_p - deltaR_m) / (2. * dphi);
//         else throw std::runtime_error("Geometric3D::getNablaDeltaR || Undefined gradient. " 
//             "Check model or implement getNablaDeltaR() explicitly.");
//         deltaR_phi /= r * sint;    
//     }
//             
//     return true;
// }


