// NuWisdom.cpp
// created by Kuangdai on 9-Oct-2016 
// Wisdom of Nu(s, z) field

#include "NuWisdom.h"
#include "Parameters.h"
#include "XMPI.h"
#include "NetCDF_Reader.h"
#include "NetCDF_ReaderAscii.h"
#include "NetCDF_Writer.h"
#include "eigenp.h"

void NuWisdom::insert(double s, double z, int nu_learn, int nu_orign) {
    RTreePoint newPoint(s, z);
    std::array<int, 2> nu = {nu_learn, nu_orign};
    mRTree.insert(std::make_pair(newPoint, nu));
}

void NuWisdom::writeToFile(const std::string &fname) const {
    if (XMPI::root()) {
        #ifndef NDEBUG
            Eigen::internal::set_is_malloc_allowed(true);
        #endif
        RDMatXX_RM data(mRTree.size(), 4);
        #ifndef NDEBUG
            Eigen::internal::set_is_malloc_allowed(false);
        #endif
        int row = 0;
        for(auto const &value: mRTree) {
            data(row, 0) = value.first.get<0>();
            data(row, 1) = value.first.get<1>();
            data(row, 2) = value.second[0] * 1.;
            data(row, 3) = value.second[1] * 1.;
            row++;
        }
        NetCDF_Writer ncw;
        ncw.open(fname, true);
        std::vector<size_t> dims;
        dims.push_back(mRTree.size());
        dims.push_back(4);
        ncw.defineVariable<double>("axisem3d_wisdom", dims);
        ncw.writeVariableWhole("axisem3d_wisdom", data);
        ncw.close();
    }
}

void NuWisdom::readFromFile(const std::string &fname) {
    RDMatXX data;
    if (XMPI::root()) {
        RDMatXX dataRead;
        if (NetCDF_Reader::checkNetCDF_isAscii(fname)) {
            NetCDF_ReaderAscii reader;
            reader.open(fname);
            reader.read2D("axisem3d_wisdom", dataRead);
            reader.close();
        } else {
            NetCDF_Reader reader;
            reader.open(fname);
            reader.read2D("axisem3d_wisdom", dataRead);
            reader.close();
        }
        if (dataRead.cols() == 3) {
            data = RDMatXX(dataRead.rows(), 4);
            data.leftCols(3) = dataRead;
            data.col(3) = data.col(2);
        } else if (dataRead.cols() == 4) {
            data = dataRead;
        } else {
            throw std::runtime_error("NuWisdom::readFromFile || "
                "Inconsistent dimensions for wisdom data, must be a matrix of 3 or 4 columns "
                "|| File = " + fname);
        }
    }
    XMPI::bcastEigen(data);
    
    // pop rtree
    mRTree.clear();
    for (int i = 0; i < data.rows(); i++) {
        // no need to check duplicated here
        RTreePoint newPoint(data(i, 0), data(i, 1));
        int nu_learn = round(data(i, 2));
        int nu_orign = round(data(i, 3));
        std::array<int, 2> nu = {nu_learn, nu_orign};
        mRTree.insert(std::make_pair(newPoint, nu));
    }
}

int NuWisdom::getNu(double s, double z, int numSamples) const {
    RTreePoint target(s, z);
    const std::vector<RTreeValue> &nearest = queryKNN(target, numSamples);
    if (nearest.size() == 0) {
        throw std::runtime_error("NuWisdom::getNu || Wisdom is empty.");
    }
    
    double nuTarget = 0.;
    double distTotal = 0.;
    for (int i = 0; i < nearest.size(); i++) {
        double dist = boost::geometry::distance(target, nearest[i].first);
        if (dist < tinyDouble) {
            return nearest[i].second[0];
        }
        distTotal += 1. / dist;
        nuTarget += nearest[i].second[0] / dist;
    }
    return round(nuTarget / distTotal);
}

int NuWisdom::getMaxNu() const {
    int nu_max = -1;
    for (RTreeValue const &v: mRTree) {
        int nu = v.second[0];
        nu_max = std::max(nu_max, nu);
    }
    return nu_max;
}

double NuWisdom::getCompressionRatio() const {
    double nu_learn = 0.;
    double nu_orign = 0.;
    for (RTreeValue const &v: mRTree) {
        nu_learn += v.second[0] * 1.;
        nu_orign += v.second[1] * 1.;
    }
    return nu_learn / nu_orign;
}

std::vector<RTreeValue> NuWisdom::queryKNN(const RTreePoint &target, int number) const {
    std::vector<RTreeValue> returned_values; 
    mRTree.query(boost::geometry::index::nearest(target, number),
        std::back_inserter(returned_values));
    return returned_values;    
}

LearnParameters::LearnParameters(const Parameters &par) {
    mInvoked = par.getValue<bool>("NU_WISDOM_LEARN");
    mCutoff = par.getValue<double>("NU_WISDOM_LEARN_EPSILON");
    if (mCutoff > .1) {
        mCutoff = .1;
    }
    if (mCutoff < 1e-5) {
        mCutoff = 1e-5;
    }
    mInterval = par.getValue<int>("NU_WISDOM_LEARN_INTERVAL");
    if (mInterval <= 0) {
        mInterval = 5;
    }
    mFileName = par.getValue<std::string>("NU_WISDOM_LEARN_OUTPUT");
    mFileName = Parameters::sOutputDirectory + "/" + mFileName;
}
