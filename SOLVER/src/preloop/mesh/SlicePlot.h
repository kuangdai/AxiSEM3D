// SlicePlot.h
// created by Kuangdai on 26-Nov-2016 
// plot a slice

#pragma once

#include <string>
#include <vector>
#include "eigenp.h"
class Parameters;
class Mesh;
class Domain;

class SlicePlot {
    
public:
    enum PropertyRefTypes {Property1D, Property3D, PropertyPerturb};
    enum SampleTypes {Center, Vertex, GLLPnt};
    
    SlicePlot(const std::string &params, const Mesh *mesh);
    
    void plotUnweighted() const;
    void plotWeighted() const;
    void plotEleType(const Domain &domain) const;
    void plotMeasured(const RDColX &cost) const;
        
    // build from input parameters
    static void buildInparam(std::vector<SlicePlot *> &splots, 
        const Parameters &par, const Mesh *mesh, int verbose);   
    
private:
    void plotPhysicalProperties() const;
    void plotNu() const;
    void plotRank(bool weighted) const;
    
    int numArgs(const std::string &parName) const;
    bool isPhysical(const std::string &parName) const;
    std::string verbose(char sep = ' ') const; 
    
    std::string mParName;
    SampleTypes mSampleType = SampleTypes::Center;
    double mLat = 0., mLon = 0., mPhi = 0.;
    PropertyRefTypes mRefType = PropertyRefTypes::Property1D;
    
    const Mesh *mMesh;
};