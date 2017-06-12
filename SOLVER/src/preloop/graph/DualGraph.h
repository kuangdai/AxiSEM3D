// DualGraph.h
// created by Kuangdai on 18-Jun-2016 
// metis dual graph

#pragma once

#include <eigenp.h>

// domain decomposition option
struct DecomposeOption {
    RDColX mElemWeights = RDColX::Zero(0);
    double mImbalance = 0.01;
    int mProcInterval = 1;
    int mNCutsPerProc = 1;
};

class DualGraph {
    
public:    
    // form neighbourhood of connectivity
    static void formNeighbourhood(const IMatX4 &connectivity, int ncommon, 
        std::vector<IColX> &neighbours);
    
    // domain decomposition
    static void decompose(const IMatX4 &connectivity, const DecomposeOption &option, 
        IColX &elemToProc);
    
private:
    static void formAdjacency(const IMatX4 &connectivity, int ncommon, int *&xadj, int *&adjncy);
    static void freeAdjacency(int *&xadj, int *&adjncy);
    static void metisError(const int retval, const std::string &func_name);
    static void check_idx_t();
};


