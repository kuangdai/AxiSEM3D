// NuWisdom.h
// created by Kuangdai on 9-Oct-2016 
// Wisdom of Nu(s, z) field

#pragma once 

class Parameters;

#include <string>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/rtree.hpp>

typedef boost::geometry::model::point<double, 2, boost::geometry::cs::cartesian> RTreePoint;
typedef std::pair<RTreePoint, std::array<int, 2>> RTreeValue;

class NuWisdom {
public:
    NuWisdom() {};
    ~NuWisdom() {};
    
    void insert(double s, double z, int nu_learn, int nu_orign);
    void writeToFile(const std::string &fname) const;
    void readFromFile(const std::string &fname);
    int getNu(double s, double z, int numSamples) const;
    int getMaxNu() const;
    double getCompressionRatio() const;
    
private:
    // KNN query
    std::vector<RTreeValue> queryKNN(const RTreePoint &target, int number) const;
    // rtree
    boost::geometry::index::rtree<RTreeValue, boost::geometry::index::quadratic<16>> mRTree;   
};

// learning options
struct LearnParameters {
    LearnParameters(const Parameters &par);
    bool mInvoked;
    double mCutoff;
    int mInterval;
    std::string mFileName;
};

