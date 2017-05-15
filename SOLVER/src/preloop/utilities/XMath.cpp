// XMath.cpp
// created by Kuangdai on 9-May-2016 
// miscellaneous math tools

#include "XMath.h"
#include "PreloopFFTW.h"

void XMath::makeClose(double &a, double &b) {
    if (a - b > pi) {
        b += 2. * pi;
    }
    if (b - a > pi) {
        a += 2. * pi;
    }
}

double XMath::findClosestDist(const std::vector<RDCol2> &crds) {
    size_t npoint = crds.size();
    double distMin = DBL_MAX;
    for (size_t i = 0; i < npoint; i++) {
        for (size_t j = i + 1; j < npoint; j++) {
            double dist = (crds[i] - crds[j]).norm();
            if (dist < distMin) {
                distMin = dist;
            }
        }
    }
    return distMin;
}

void XMath::interpLagrange(double target, size_t nbases, 
    const std::vector<double> &bases, std::vector<double> &results) {
    for (size_t dgr = 0; dgr < nbases; dgr++) {
        double prod1 = 1.;
        double prod2 = 1.;
        double x0 = bases[dgr];
        for (size_t i = 0; i < nbases; i++) {
            if (i != dgr) {
                double x = bases[i];
                prod1 *= target - x;
                prod2 *= x0 - x;
            }
        }
        results[dgr] = prod1 / prod2;
    }
}

void XMath::gaussianSmoothing(RDColX &data, size_t order, double dev, bool period) {
    if (data.size() == 0) return;
    order = std::min(order, (data.size() + 1) / 2 - 1);
    if (order == 0) return;
    dev *= order;
    
    // gaussian kernel
    RDColX gaussian(order * 2 + 1);
    for (size_t i = 0; i <= order; i++) {
        gaussian(order + i) = exp(- .5 * i * i / (dev * dev));
        gaussian(order - i) = gaussian(order + i);
    }
    gaussian /= gaussian.sum();
    
    // convolve
    RDColX result = RDColX::Zero(data.size());
    for (size_t i = 0; i < data.size(); i++) {
        for (int j = -order; j <= order; j++) {
            int k = i + j;
            if (period) {
                while (k < 0) {
                    k += data.size();
                }
                while (k > data.size() - 1) {
                    k -= data.size();
                }
            } else {
                // using fixed padding
                if (k < 0) {
                    k = 0;
                }
                if (k > data.size() - 1) {
                    k = data.size() - 1;
                }
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
    for (size_t i = 0; i < data.rows(); i++) {
        RDColX temp = data.row(i).transpose();
        gaussianSmoothing(temp, orderRow(i), devRow(i), periodRow);
        data.row(i) = temp.transpose();
    }
    for (size_t i = 0; i < data.cols(); i++) {
        RDColX temp = data.col(i);
        gaussianSmoothing(temp, orderCol(i), devCol(i), periodCol);
        data.col(i) = temp;
    }
}

RMatPP XMath::castToSolver(const RDMatPP &mp) {
    RMatPP ms = RMatPP::Zero();
    ms.topLeftCorner(nPntEdge, nPntEdge) = mp.cast<Real>();
    return ms;
}

RDColX XMath::trigonResampling(size_t newSize, const RDColX &original) {
    size_t nslices = original.rows();
    if (newSize == nslices) {
        return original;
    }
    if (equalRows(original)) {
        return RDColX::Constant(newSize, original(0));
    }
    
    // fft
    PreloopFFTW::getR2C_RMat(nslices) = original;
    PreloopFFTW::computeR2C(nslices);
    CDColX &fourier = PreloopFFTW::getR2C_CMat(nslices);
    
    // densed sampling
    double dphi = 2. * pi / newSize;
    RDColX densed(newSize);
    for (size_t islice = 0; islice < newSize; islice++) {
        double phi = islice * dphi;
        double value_phi = fourier(0).real();
        for (size_t alpha = 1; alpha < fourier.size(); alpha++) {
            double factor = (nslices % 2 == 0 && alpha == fourier.size() - 1) ? 1. : 2.;
            value_phi += factor * (fourier(alpha) * exp(alpha * phi * iid)).real();
        }
        densed(islice) = value_phi;
    }
    return densed; 
}

RDColX XMath::linearResampling(size_t newSize, const RDColX &original) {
    size_t nslices = original.rows();
    if (newSize == nslices) {
        return original;
    }
    if (equalRows(original)) {
        return RDColX::Constant(newSize, original(0));
    }
    
    // densed sampling
    double dphi = 2. * pi / newSize;
    double dphi_orig = 2. * pi / nslices;
    RDColX densed(newSize);
    for (size_t islice = 0; islice < newSize; islice++) {
        double phi = islice * dphi;
        size_t loc0 = (size_t)(phi / dphi_orig);
        double phi0 = loc0 * dphi_orig;
        size_t loc1 = loc0 + 1;
        if (loc1 == nslices) {
            loc1 = 0;
        }
        densed(islice) = (original(loc1) - original(loc0)) 
            / dphi_orig * (phi - phi0) + original(loc0);
    }
    return densed; 
}

RDRowN XMath::computeFourierAtPhi(const RDMatXN &data, double phi) {
    size_t nslices = data.rows();
    RDRowN result;
    for (size_t col = 0; col < nPntElem; col++) {
        PreloopFFTW::getR2C_RMat(nslices) = data.col(col);
        PreloopFFTW::computeR2C(nslices);
        CDColX &fourier = PreloopFFTW::getR2C_CMat(nslices);
        double value_phi = fourier(0).real();
        for (size_t alpha = 1; alpha < fourier.size(); alpha++) {
            double factor = (nslices % 2 == 0 && alpha == fourier.size() - 1) ? 1. : 2.;
            value_phi += factor * (fourier(alpha) * exp(alpha * phi * iid)).real();
        }
        result(col) = value_phi;
    }
    return result;
}
