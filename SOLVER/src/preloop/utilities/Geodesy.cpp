// Geodesy.cpp
// created by Kuangdai on 15-May-2017 
// geodetic tools

#include "Geodesy.h"

double Geodesy::sFlattening = 0.;
double Geodesy::sROuter = 6371e3;
RDColX Geodesy::sEllipKnots = RDColX::Zero(0);
RDColX Geodesy::sEllipCoeffs = RDColX::Zero(0);


void Geodesy::rtheta(const RDCol2 &sz, double &r, double &theta) {
    r = sz.norm();
    theta = (r < tinyDouble) ? 0. : acos(sz(1) / r);
}

RDCol2 Geodesy::rtheta(const RDCol2 &sz) {
    RDCol2 rt;
    rtheta(sz, rt(0), rt(1));
    return rt;
}

double Geodesy::atan4(double y, double x, bool &defined) {
    if (sqrt(x * x + y * y) < tinyDouble) {
        defined = false;
        return 0.;
    }
    double t = atan2(y, x);
    if (t < 0.) {
        t += pi;
    }
    if (y < 0.) {
        t += pi;
    }
    defined = true;
    return t;
}

RDCol3 Geodesy::toCartesian(const RDCol3 &rtp) {
    RDCol3 xyz;
    xyz(0) = rtp(0) * sin(rtp(1)) * cos(rtp(2));
    xyz(1) = rtp(0) * sin(rtp(1)) * sin(rtp(2));
    xyz(2) = rtp(0) * cos(rtp(1));
    return xyz;
}

RDCol3 Geodesy::toSpherical(const RDCol3 &xyz, bool &defined) {
    RDCol3 rtp;
    rtp(0) = xyz.norm();
    rtp(1) = rtp(0) < tinyDouble ? 0. : acos(xyz(2) / rtp(0));
    rtp(2) = atan4(xyz(1), xyz(0), defined);
    return rtp;
}

RDMat33 Geodesy::rotationMatrix(double theta, double phi) {
    RDMat33 Q;
    Q(0, 0) = cos(theta) * cos(phi);
    Q(0, 1) = -sin(phi);
    Q(0, 2) = sin(theta) * cos(phi);
    Q(1, 0) = cos(theta) * sin(phi);
    Q(1, 1) = cos(phi);
    Q(1, 2) = sin(theta) * sin(phi);
    Q(2, 0) = -sin(theta);
    Q(2, 1) = 0.;
    Q(2, 2) = cos(theta);
    return Q;
}

double Geodesy::lat2Theta_r(double lat, double radius) {
    double flattening = getFlattening(radius);
    double one_minus_f_squared = (1. - flattening) * (1. - flattening);
    // take care of poles
    double limit_round = 90. - tinyDouble;
    double limit_bound = 90.1;
    if (lat > limit_round) {
        if (lat > limit_bound) {
            throw std::runtime_error("Geodesy::lat2Theta || Latitude out of range [-90, 90]");
        }
        lat = limit_round;
    }
    if (lat < -limit_round) {
        if (lat < -limit_bound) {
            throw std::runtime_error("Geodesy::lat2Theta || Latitude out of range [-90, 90]");
        }
        lat = -limit_round;
    }
    return pi / 2. - atan(one_minus_f_squared * tan(lat * degree));
}

double Geodesy::lon2Phi(double lon) {
    if (lon < 0.) {
        lon += 360.;
    }
    return lon * degree;
}

double Geodesy::theta2Lat_r(double theta, double radius) {
    double flattening = getFlattening(radius);
    double inv_one_minus_f_squared = 1. / (1. - flattening) / (1. - flattening);
    // take care of poles
    double lat = pi / 2. - theta;
    double limit_round = (90. - tinyDouble) * degree;
    double limit_bound = 90.1 * degree;
    if (lat > limit_round) {
        if (lat > limit_bound) {
            throw std::runtime_error("Geodesy::theta2Lat || Theta out of range [0, pi]");
        }
        lat = limit_round;
    }
    if (lat < -limit_round) {
        if (lat < -limit_bound) {
            throw std::runtime_error("Geodesy::theta2Lat || Theta out of range [0, pi]");
        }
        lat = -limit_round;
    }
    return atan(inv_one_minus_f_squared * tan(lat)) / degree;
}

double Geodesy::phi2Lon(double phi) {
    if (phi > pi) {
        phi -= 2. * pi;
    }
    return phi / degree;
}

RDCol3 Geodesy::rotateSrc2Glob(const RDCol3 &rtpS, double srclat, double srclon, double srcdep) {
    const RDCol3 &xyzS = toCartesian(rtpS);
    const RDCol3 &xyzG = rotationMatrix(lat2Theta_d(srclat, srcdep), lon2Phi(srclon)) * xyzS;
    bool defined = true;
    RDCol3 rtpG = toSpherical(xyzG, defined);
    if (!defined) {
        rtpG(2) = rtpS(2);
    }
    return rtpG;
}

RDCol3 Geodesy::rotateGlob2Src(const RDCol3 &rtpG, double srclat, double srclon, double srcdep) {
    const RDCol3 &xyzG = toCartesian(rtpG);
    const RDCol3 &xyzS = rotationMatrix(lat2Theta_d(srclat, srcdep), lon2Phi(srclon)).transpose() * xyzG;
    bool defined = true;
    RDCol3 rtpS = toSpherical(xyzS, defined);
    if (!defined) {
        rtpS(2) = rtpG(2);
    }
    return rtpS;
}

double Geodesy::backAzimuth(double srclat, double srclon, double srcdep,
    double reclat, double reclon, double recdep) {    
    // event
    double srcTheta = lat2Theta_d(srclat, srcdep);
    double srcPhi = lon2Phi(srclon);
    double d = sin(srcPhi);
    double e = -cos(srcPhi);
    double f = -sin(srcTheta);
    double c = cos(srcTheta);
    double a = f * e;
    double b = -f * d;
    // station
    double recTheta = lat2Theta_d(reclat, recdep);
    double recPhi = lon2Phi(reclon);
    double d1 = sin(recPhi);
    double e1 = -cos(recPhi);
    double f1 = -sin(recTheta);
    double c1 = cos(recTheta);
    double g1 = -c1 * e1;
    double h1 = c1 * d1;
    // baz
    double ss = (a - d1) * (a - d1) + (b - e1) * (b - e1) + c * c - 2.;
    double sc = (a - g1) * (a - g1) + (b - h1) * (b - h1) + (c - f1) * (c - f1) - 2.;
    bool defined; // useless
    return atan4(ss, sc, defined);
}

void Geodesy::setup(double router, double flattening, 
                     const RDColX &ellip_knots, 
                     const RDColX &ellip_coeffs) {
    sROuter = router;
    sFlattening = flattening;
    sEllipKnots = ellip_knots;
    sEllipCoeffs = ellip_coeffs;
}

double Geodesy::getFlattening(double r) {
    double f = 1.;
    double r_ref = r / sROuter;
    int nknots = sEllipKnots.size();
    for (int i = 1; i < nknots; i++) {
        if (r_ref <= sEllipKnots(i)) {
            f = (sEllipCoeffs(i) - sEllipCoeffs(i - 1)) 
                / (sEllipKnots(i) - sEllipKnots(i - 1))
                * (r_ref - sEllipKnots(i - 1)) 
                + sEllipCoeffs(i - 1);
            break;
        }
    }
    return f * sFlattening;
}

