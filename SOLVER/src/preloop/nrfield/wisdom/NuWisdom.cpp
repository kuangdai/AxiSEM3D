// NuWisdom.cpp
// created by Kuangdai on 9-Oct-2016 
// Wisdom of Nu(s, z) field

#include "NuWisdom.h"
#include "Parameters.h"
#include <boost/algorithm/string.hpp>
#include <fstream>
#include "XMPI.h"

void NuWisdom::insert(double s, double z, int nu_learn, int nu_orign) {
    RTreePoint newPoint(s, z);
    std::array<int, 2> nu = {nu_learn, nu_orign};
    mRTree.insert(std::make_pair(newPoint, nu));
}

void NuWisdom::writeToFile(const std::string &fname, bool binary) const {
    std::fstream::openmode mode = std::fstream::out;
    if (binary) mode = mode | std::fstream::binary;
    std::fstream fs(fname, mode);
    if (binary) {
        for(auto const &value: mRTree) {
            float s = (float)value.first.get<0>();
            float z = (float)value.first.get<1>();
            fs.write((char *) &s, sizeof(float));
            fs.write((char *) &z, sizeof(float));
            fs.write((char *) &(value.second[0]), sizeof(int));
            fs.write((char *) &(value.second[1]), sizeof(int));
        }
    } else {
        for(auto const &value: mRTree) {
            fs << (float)value.first.get<0>() << " ";
            fs << (float)value.first.get<1>() << " ";
            fs << value.second[0] << " " << value.second[1] << std::endl;
        }
    }
    fs.close();
}

void NuWisdom::readFromFile(const std::string &fname, bool binary) {
    // read data
    std::vector<float> buffer;
    if (XMPI::root()) {
        std::fstream::openmode mode = std::fstream::in;
        if (binary) mode = mode | std::fstream::binary;
        std::fstream fs(fname, mode);
        if (!fs) throw std::runtime_error("NuWisdom::readFromFile || Error opening Wisdom file: ||" + fname);
        float s, z; int nul, nuo;
        if (binary) {
            while (true) {
                fs.read((char *) &s, sizeof(float));
                fs.read((char *) &z, sizeof(float));
                fs.read((char *) &nul, sizeof(int));
                fs.read((char *) &nuo, sizeof(int));
                if (!fs) break;
                buffer.push_back(s);
                buffer.push_back(z);
                buffer.push_back(nul);
                buffer.push_back(nuo);
            }
        } else {
            while (fs >> s >> z >> nul >> nuo) {
                buffer.push_back(s);
                buffer.push_back(z);
                buffer.push_back(nul);
                buffer.push_back(nuo);
            }
        }
        fs.close();
    }
    XMPI::bcast(buffer);
    
    // pop rtree
    mRTree.clear();
    for (int i = 0; i < buffer.size() / 4; i++) {
        // no need to check duplicated here
        RTreePoint newPoint(buffer[i * 4], buffer[i * 4 + 1]);
        int nu_learn = round(buffer[i * 4 + 2]);
        int nu_orign = round(buffer[i * 4 + 3]);
        std::array<int, 2> nu = {nu_learn, nu_orign};
        mRTree.insert(std::make_pair(newPoint, nu));
    }
}

int NuWisdom::getNu(double s, double z, int numSamples) const {
    RTreePoint target(s, z);
    const std::vector<RTreeValue> &nearest = queryKNN(target, numSamples);
    if (nearest.size() == 0) 
        throw std::runtime_error("NuWisdom::getNu || Wisdom is empty.");
    
    double nuTarget = 0.;
    double distTotal = 0.;
    for (int i = 0; i < nearest.size(); i++) {
        double dist = boost::geometry::distance(target, nearest[i].first);
        if (dist < tinyDouble) return nearest[i].second[0];
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
    int nu_learn = 0;
    int nu_orign = 0;
    for (RTreeValue const &v: mRTree) {
        nu_learn += v.second[0];
        nu_orign += v.second[1];
    }
    return 1. * nu_learn / nu_orign;
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
    if (mCutoff > .1) mCutoff = .1;
    if (mCutoff < 1e-5) mCutoff = 1e-5;
    mInterval = par.getValue<int>("NU_WISDOM_LEARN_INTERVAL");
    if (mInterval <= 0) mInterval = 10;
    mFileName = par.getValue<std::string>("NU_WISDOM_LEARN_OUTPUT");
    mFileName = Parameters::sOutputDirectory + "/" + mFileName;
}
