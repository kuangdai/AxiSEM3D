// DualGraph.cpp
// created by Kuangdai on 18-Jun-2016 
// metis dual graph

#include "DualGraph.h"
#include "XMPI.h"
#include <metis.h>
#include <algorithm>
#include <limits>
#include <numeric>

void DualGraph::formNeighbourhood(const IMatX4 &connectivity, int ncommon, 
    std::vector<IColX> &neighbours) {
    // form graph
    int *xadj, *adjncy;
    formAdjacency(connectivity, ncommon, xadj, adjncy);
    
    // convert CRS format to neighbours
    int nelem = connectivity.rows();
    neighbours.clear();
    neighbours.reserve(nelem);
    for (int i = 0; i < nelem; i++) {
        int nNeighb = xadj[i + 1] - xadj[i];
        IColX neighb(nNeighb);
        for (int j = 0; j < nNeighb; j++) {
            neighb(j) = adjncy[xadj[i] + j];
        }
        neighbours.push_back(neighb);
    }

    // delete graph 
    freeAdjacency(xadj, adjncy);
}

void DualGraph::decompose(const IMatX4 &connectivity, const DecomposeOption &option, 
    IColX &elemToProc) {
    // size
    int nelem = connectivity.rows();    
    elemToProc = IColX::Zero(nelem);
    int nproc = XMPI::nproc();
    if (nproc == 1) {
        return;
    }
    
    int objval = std::numeric_limits<int>::max();
    if (XMPI::rank() % option.mProcInterval == 0) {
        // form graph
        int *xadj, *adjncy;
        formAdjacency(connectivity, 2, xadj, adjncy);
        
        // check size of weights
        bool welem = option.mElemWeights.size() > 0;
        if (welem && option.mElemWeights.size() != nelem) {
            throw std::runtime_error("DualGraph::decompose || "
                "Incompatible size of element weights.");
        }
        
        // metis options
        int metis_option[METIS_NOPTIONS];
        metisError(METIS_SetDefaultOptions(metis_option), "METIS_SetDefaultOptions");
        metis_option[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;
        metis_option[METIS_OPTION_CONTIG] = 1;
        metis_option[METIS_OPTION_NCUTS] = option.mNCutsPerProc;
        metis_option[METIS_OPTION_SEED] = XMPI::rank(); // generate different partitions
        
        // prepare input
        int ncon = 1;
        float ubvec = (float)(1. + option.mImbalance);
        int *vwgt = NULL;
        IColX elemWeightsInt;
        double imax = std::numeric_limits<int>::max() * .9;
        double sum = 0.;
        if (welem) {
            sum += option.mElemWeights.sum();
            elemWeightsInt = (option.mElemWeights / sum * imax).array().round().matrix().cast<int>();
            vwgt = elemWeightsInt.data();
        }
        
        // run
        metisError(METIS_PartGraphKway(&nelem, &ncon, xadj, adjncy, 
            vwgt, NULL, NULL, &nproc, NULL, &ubvec, 
            metis_option, &objval, elemToProc.data()), 
            "METIS_PartGraphKway");
         
        // free memory 
        freeAdjacency(xadj, adjncy);
    }
    
    // find best result
    std::vector<int> objall;
    XMPI::gather(objval, objall, true);
    int proc_min = std::min_element(objall.begin(), objall.end()) - objall.begin();
    XMPI::bcastEigen(elemToProc, proc_min);
}

void DualGraph::formAdjacency(const IMatX4 &connectivity, int ncommon, int *&xadj, int *&adjncy) {
    DualGraph::check_idx_t();
    
    // get unique node list
    int nelem = connectivity.rows();
    std::vector<int> nodes(connectivity.data(), connectivity.data() + nelem * 4);
    std::sort(nodes.begin(), nodes.end());
    auto itu = std::unique(nodes.begin(), nodes.end());
    int nnode = std::distance(nodes.begin(), itu);
    
    // convert connectivity to CSR format
    IColX eptr = IColX(nelem + 1);
    for (int i = 0; i < nelem + 1; i++) {
        eptr(i) = i * 4;
    }
    IColX eind = IColX(nelem * 4);
    for (int i = 0; i < nelem; i++) {
        for (int j = 0; j < 4; j++) {
            eind(i * 4 + j) = std::distance(nodes.begin(), 
                std::lower_bound(nodes.begin(), itu, connectivity(i, j)));
        }
    }

    // convert to dual graph
    int numflag = 0;
    metisError(METIS_MeshToDual(&nelem, &nnode, eptr.data(), eind.data(), 
        &ncommon, &numflag, &xadj, &adjncy), "METIS_MeshToDual"); 
}

void DualGraph::freeAdjacency(int *&xadj, int *&adjncy) {
    // free CRS created by metis 
    metisError(METIS_Free(xadj), "METIS_Free"); 
    metisError(METIS_Free(adjncy), "METIS_Free"); 
}

void DualGraph::metisError(const int retval, const std::string &func_name) {
    if (retval != METIS_OK) {
        throw std::runtime_error("DualGraph::metisError || "
            "Error in metis function: " + func_name);
    }
}

void DualGraph::check_idx_t() {
    if (IDXTYPEWIDTH != 32 || REALTYPEWIDTH != 32) {
        throw std::runtime_error("DualGraph::check_idx_t || Incompatible METIS build ||"
            "Please re-install METIS with 32-bit configuration. ||"
            "Or edit METIS_ROOT in CMakeLists.txt to locate a 32-bit build.");    
    }
}
