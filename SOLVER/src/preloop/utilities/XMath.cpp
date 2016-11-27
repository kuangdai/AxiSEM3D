// XMath.cpp
// created by Kuangdai on 9-May-2016 
// small math tools

#include "XMath.h"
#include "PreloopFFTW.h"

void XMath::rtheta(const RDCol2 &sz, double &r, double &theta) {
    r = sz.norm();
    theta = (r < tinyDouble) ? 0. : acos(sz(1) / r);
}

RDCol2 XMath::rtheta(const RDCol2 &sz) {
    RDCol2 rt;
    rtheta(sz, rt(0), rt(1));
    return rt;
}

double XMath::atan4(double y, double x, bool &defined) {
    if (sqrt(x * x + y * y) < tinyDouble) {
        defined = false;
        return 0.;
    }
    double t = atan2(y, x);
    if (t < 0.) t += pi;
    if (y < 0.) t += pi;
    defined = true;
    return t;
}

RDMat33 XMath::rotationMatrix(double theta, double phi) {
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

RDMat33 XMath::sphericalBasis(double theta, double phi) {
    RDMat33 Q;
    Q(0, 0) = sin(theta) * cos(phi);
    Q(0, 1) = sin(theta) * sin(phi);
    Q(0, 2) = cos(theta);
    Q(1, 0) = cos(theta) * cos(phi);
    Q(1, 1) = cos(theta) * sin(phi);
    Q(1, 2) = -sin(theta);
    Q(2, 0) = -sin(phi);
    Q(2, 1) = cos(phi);
    Q(2, 2) = 0.;
    return Q;
}

RDCol3 XMath::toCartesian(const RDCol3 &rtp) {
    RDCol3 xyz;
    xyz(0) = rtp(0) * sin(rtp(1)) * cos(rtp(2));
    xyz(1) = rtp(0) * sin(rtp(1)) * sin(rtp(2));
    xyz(2) = rtp(0) * cos(rtp(1));
    return xyz;
}

RDCol3 XMath::toSpherical(const RDCol3 &xyz, bool &defined) {
    RDCol3 rtp;
    rtp(0) = xyz.norm();
    rtp(1) = rtp(0) < tinyDouble ? 0. : acos(xyz(2) / rtp(0));
    rtp(2) = atan4(xyz(1), xyz(0), defined);
    return rtp;
}

double XMath::lat2Theta(double lat, double depth) {
    double flattening = getFlattening(sROuter - depth);
    double one_minus_f_squared = (1. - flattening) * (1. - flattening);
    return pi / 2. - atan(one_minus_f_squared * tan(lat * degree));
}

double XMath::lon2Phi(double lon) {
    if (lon < 0.) lon += 360.;
    return lon * degree;
}

double XMath::theta2Lat(double theta, double depth) {
    double flattening = getFlattening(sROuter - depth);
    double inv_one_minus_f_squared = 1. / (1. - flattening) / (1. - flattening);
    return atan(inv_one_minus_f_squared * tan(pi / 2. - theta)) / degree;
}

double XMath::phi2Lon(double phi) {
    if (phi > pi) phi -= 2. * pi;
    return phi / degree;
}

RDCol3 XMath::rotateSrc2Glob(const RDCol3 &rtpS, double srclat, double srclon, double srcdep) {
    const RDCol3 &xyzS = toCartesian(rtpS);
    const RDCol3 &xyzG = rotationMatrix(lat2Theta(srclat, srcdep), lon2Phi(srclon)) * xyzS;
    bool defined;
    RDCol3 rtpG = toSpherical(xyzG, defined);
    if (!defined) rtpG(2) = rtpS(2);
    return rtpG;
}

RDCol3 XMath::rotateGlob2Src(const RDCol3 &rtpG, double srclat, double srclon, double srcdep) {
    const RDCol3 &xyzG = toCartesian(rtpG);
    const RDCol3 &xyzS = rotationMatrix(lat2Theta(srclat, srcdep), lon2Phi(srclon)).transpose() * xyzG;
    bool defined;
    RDCol3 rtpS = toSpherical(xyzS, defined);
    if (!defined) rtpS(2) = rtpG(2);
    return rtpS;
}

double XMath::backAzimuth(double srclat, double srclon, double srcdep,
    double reclat, double reclon, double recdep) {    
    // event
    double srcTheta = lat2Theta(srclat, srcdep);
    double srcPhi = lon2Phi(srclon);
    double d = sin(srcPhi);
    double e = -cos(srcPhi);
    double f = -sin(srcTheta);
    double c = cos(srcTheta);
    double a = f * e;
    double b = -f * d;
    // station
    double recTheta = lat2Theta(reclat, recdep);
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

void XMath::makeClose(double &a, double &b) {
    if (a - b > pi) b += 2. * pi;
    if (b - a > pi) a += 2. * pi;
}

double XMath::findClosestDist(const std::vector<RDCol2> &crds) {
    int npoint = crds.size();
    double distMin = 1e100;
    for (int i = 0; i < npoint; i++) {
        for (int j = i + 1; j < npoint; j++) {
            double dist = (crds[i] - crds[j]).norm();
            if (dist < distMin) distMin = dist;
        }
    }
    return distMin;
}

void XMath::interpLagrange(double target, int nbases, const double *bases, double *results) {
    for (int dgr = 0; dgr < nbases; dgr++) {
        double prod1 = 1.;
        double prod2 = 1.;
        double x0 = bases[dgr];
        for (int i = 0; i < nbases; i++) {
            if (i != dgr) {
                double x = bases[i];
                prod1 *= target - x;
                prod2 *= x0 - x;
            }
        }
        results[dgr] = prod1 / prod2;
    }
}

void XMath::gaussianSmoothing(RDColX &data, int order, double dev, bool period) {
    if (order <= 0) return;
    order = std::min(order, (int)((data.size() + 1) / 2 - 1));
    dev *= order;
    
    // gaussian kernel
    RDColX gaussian(order * 2 + 1);
    for (int i = 0; i <= order; i++) {
        gaussian(order + i) = exp(- .5 * i * i / (dev * dev));
        gaussian(order - i) = gaussian(order + i);
    }
    gaussian /= gaussian.sum();
    
    // convolve
    RDColX result = RDColX::Zero(data.size());
    for (int i = 0; i < data.size(); i++) {
        for (int j = -order; j <= order; j++) {
            int k = i + j;
            if (period) {
                while (k < 0) k += data.size();
                while (k > data.size() - 1) k -= data.size();
            } else {
                // using fixed padding
                if (k < 0) k = 0;
                if (k > data.size() - 1) k = data.size() - 1;
            }
            result(i) += gaussian(j + order) * data(k);
        }
    }
    // assign
    data = result;
}

void XMath::gaussianSmoothing(RDMatXX &data, 
    IColX orderRow, RDColX devRow, bool periodRow, 
    IColX orderCol, RDColX devCol, bool periodCol) {
    for (int i = 0; i < data.rows(); i++) {
        RDColX temp = data.row(i).transpose();
        gaussianSmoothing(temp, orderRow(i), devRow(i), periodRow);
        data.row(i) = temp.transpose();
    }
    for (int i = 0; i < data.cols(); i++) {
        RDColX temp = data.col(i);
        gaussianSmoothing(temp, orderCol(i), devCol(i), periodCol);
        data.col(i) = temp;
    }
}

bool XMath::isLuckyNumber(int n, bool forceOdd)
{
    int num = n;
    
    // We always hope to use even numbers that are generally faster,
    // but the Nyquist frequency sometimes causes trouble.
    // force odd
    if (forceOdd && num % 2 == 0) return false;
    
    // use even when n > 10
    if (!forceOdd && num % 2 != 0 && num > 10) return false;
    
    for (int i = 2; i <= num; i++) {  
        while(num % i == 0) {
            num /= i;
            if (i > 13) return false;
        }
    }
    num = n;
    int e = 0;
    while(num % 11 == 0) {
        num /= 11;
        e++;
    }
    num = n;
    int f = 0;
    while(num % 13 == 0) {
        num /= 13;
        f++;
    }
    if (e + f > 1) return false;
    return true;
}

int XMath::nextLuckyNumber(int n, bool forceOdd)
{
    while(true) {
        if (isLuckyNumber(n, forceOdd)) return n;
        n++;
    }
}

RMatPP XMath::castToSolver(const RDMatPP &mp) {
    RMatPP ms = RMatPP::Zero();
    ms.topLeftCorner(nPntEdge, nPntEdge) = mp.cast<Real>();
    return ms;
}

RDColX XMath::trigonResampling(int newSize, const RDColX &original) {
    int nslices = original.rows();
    if (newSize == nslices) return original;
    if (equalRows(original)) return RDColX::Constant(newSize, original(0));
    
    // fft
    PreloopFFTW::getR2C_RMat(nslices) = original;
    PreloopFFTW::computeR2C(nslices);
    CDColX &fourier = PreloopFFTW::getR2C_CMat(nslices);
    
    // densed sampling
    double dphi = 2. * pi / newSize;
    RDColX densed(newSize);
    for (int islice = 0; islice < newSize; islice++) {
        double phi = islice * dphi;
        double value_phi = fourier(0).real();
        for (int alpha = 1; alpha < fourier.size(); alpha++) {
            double factor = (nslices % 2 == 0 && alpha == fourier.size() - 1) ? 1. : 2.;
            value_phi += factor * (fourier(alpha) * exp(alpha * phi * iid)).real();
        }
        densed(islice) = value_phi;
    }
    return densed; 
}

RDColX XMath::linearResampling(int newSize, const RDColX &original) {
    int nslices = original.rows();
    if (newSize == nslices) return original;
    if (equalRows(original)) return RDColX::Constant(newSize, original(0));
    
    // densed sampling
    double dphi = 2. * pi / newSize;
    double dphi_orig = 2. * pi / nslices;
    RDColX densed(newSize);
    for (int islice = 0; islice < newSize; islice++) {
        double phi = islice * dphi;
        int loc0 = (int)(phi / dphi_orig);
        double phi0 = loc0 * dphi_orig;
        int loc1 = loc0 + 1;
        if (loc1 == nslices) loc1 = 0;
        densed(islice) = (original(loc1) - original(loc0)) / dphi_orig * (phi - phi0) + original(loc0);
    }
    return densed; 
}

RDRowN XMath::computeFourierAtPhi(const RDMatXN &data, double phi) {
    int nslices = data.rows();
    RDRowN result;
    for (int col = 0; col < nPntElem; col++) {
        PreloopFFTW::getR2C_RMat(nslices) = data.col(col);
        PreloopFFTW::computeR2C(nslices);
        CDColX &fourier = PreloopFFTW::getR2C_CMat(nslices);
        double value_phi = fourier(0).real();
        for (int alpha = 1; alpha < fourier.size(); alpha++) {
            double factor = (nslices % 2 == 0 && alpha == fourier.size() - 1) ? 1. : 2.;
            value_phi += factor * (fourier(alpha) * exp(alpha * phi * iid)).real();
        }
        result(col) = value_phi;
    }
    return result;
}

#include <boost/timer/timer.hpp>
double XMath::getClockResolution(bool user) {
    boost::timer::cpu_timer cpu;
    boost::timer::cpu_times start_time, current_time;
    if (user) {
        cpu.start();
        start_time = cpu.elapsed();
        current_time.user = start_time.user;
        while (current_time.user == start_time.user)
            current_time = cpu.elapsed();
        return (current_time.user - start_time.user) * 1.;
    } else {
        cpu.start();
        start_time = cpu.elapsed();
        current_time.wall = start_time.wall;
        while (current_time.wall == start_time.wall)
            current_time = cpu.elapsed();
        return (current_time.wall - start_time.wall) * 1.;
    }
}

double XMath::sFlattening = 0.;
double XMath::sROuter = 6371e3;
RDColX XMath::sEllipKnots = RDColX::Zero(0);
RDColX XMath::sEllipCoeffs = RDColX::Zero(0);
void XMath::setEllipticity(double flattening, double router, 
    const std::vector<double> &ellip_knots, 
    const std::vector<double> &ellip_coeffs) {
    sFlattening = flattening;
    sROuter = router;
    int nknots = ellip_knots.size();
    sEllipKnots = sEllipCoeffs = RDColX(nknots);
    for (int i = 0; i < nknots; i++) {
        sEllipKnots(i) = ellip_knots[i];
        sEllipCoeffs(i) = ellip_coeffs[i];
    }
}

double XMath::getFlattening(double r) {
    double f = 1.;
    double r1 = r / sROuter;
    int nknots = sEllipKnots.size();
    for (int i = 1; i < nknots; i++) {
        if (r1 <= sEllipKnots(i)) {
            f = (sEllipCoeffs(i) - sEllipCoeffs(i - 1)) / (sEllipKnots(i) - sEllipKnots(i - 1))
                * (r1 - sEllipKnots(i - 1)) + sEllipCoeffs(i - 1);
            break;
        }
    }
    return f * sFlattening;
}

