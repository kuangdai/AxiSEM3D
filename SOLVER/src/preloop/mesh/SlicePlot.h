// SlicePlot.h
// created by Kuangdai on 26-Nov-2016 
// plot a slice

#pragma once

#include <string>
#include <vector>
class Parameters;
class Mesh;

class SlicePlot {
    
public:
    enum PropertyTypes {Property1D, Property3D, PropertyPerturb};
    enum SampleTypes {Center, Vertex, GLLPnt};
    
    SlicePlot(const std::string &params, const Mesh *mesh);
    
    std::string verbose() const; 
    
    void plotUnweighted() const;
    void plotWeighted() const;
        
    // build from input parameters
    static void buildInparam(std::vector<SlicePlot *> &splots, 
        const Parameters &par, const Mesh *mesh);   
    
private:
    bool isSlicePar(const std::string &parName) const;
    
    std::string mParName;
    double mLat, mLon, mPhi;
    PropertyTypes mPropertyType;
    SampleTypes mSampleType;
    double mRMin, mRMax;
    
    const Mesh *mMesh;
};