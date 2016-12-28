// DualGraph.h
// created by Kuangdai on 18-Jun-2016 
// metis dual graph

#pragma once

#include <eigenp.h>

// domain decomposition option
struct DecomposeOption {
    
    DecomposeOption() {
        mDoubleConstrants = false;
        mElemWeights1 = mElemWeights2 = std::vector<double>();
        mElemCommSize = std::vector<int>();
        mImbalance1 = mImbalance2 = 0.;
        mNPartition = 0;
        mCommVol = false;
    };
    
    DecomposeOption(int nelem, bool balance_sf, double weight, int csize, 
        int npart, bool cvol, float imbal = 0.001f) {
        mDoubleConstrants = balance_sf;
        mElemWeights1 = mElemWeights2 = std::vector<double>(nelem, weight);
        mElemCommSize = std::vector<int>(nelem, csize);
        mImbalance1 = mImbalance2 = imbal;
        mNPartition = npart;
        mCommVol = cvol;
    };
    
    // balance two kinds of weights individually
    // e.g., solid/fluid, point/element
    bool mDoubleConstrants;
    // element weights [nelem]
    std::vector<double> mElemWeights1;
    std::vector<double> mElemWeights2;
    // load imbalance  
    float mImbalance1;
    float mImbalance2;
    // element communication size [nelem]
    std::vector<int> mElemCommSize;
    int mNPartition;
    bool mCommVol;
};

class DualGraph {
    
public:    
    // form neighbourhood of connectivity
    static void formNeighbourhood(const IMatX4 &connectivity, int ncommon, 
        std::vector<IColX> &neighbours);
    
    // domain decomposition
    static void decompose(const IMatX4 &connectivity, const DecomposeOption &option, int nproc, 
        IColX &elemToProc);
    
private:
    static void formAdjacency(const IMatX4 &connectivity, int ncommon, int *&xadj, int *&adjncy);
    static void freeAdjacency(int *&xadj, int *&adjncy);
    static void metisError(const int retval, const std::string &func_name);
    static void check_idx_t();
};


