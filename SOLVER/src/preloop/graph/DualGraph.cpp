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
        for (int j = 0; j < nNeighb; j++) neighb(j) = adjncy[xadj[i] + j];
        neighbours.push_back(neighb);
    }

    // delete graph 
    freeAdjacency(xadj, adjncy);
}

void DualGraph::decompose(const IMatX4 &connectivity, const DecomposeOption &option, int nproc, 
    IColX &elemToProc) {
    int nelem = connectivity.rows();    
    elemToProc = IColX::Zero(nelem);
    if (nproc == 1) return;
    
    if (XMPI::root()) {    
        // form graph
        int *xadj, *adjncy;
        formAdjacency(connectivity, 2, xadj, adjncy);
        
        // weights
        double imax = std::numeric_limits<int>::max() * .9;
        int ncon = 0;
        std::vector<int> vweight;
        float ubvec[2] = {1.f, 1.f};
        std::vector<int> vsize = option.mElemCommSize;
        if (option.mDoubleConstrants) {
            ncon = 2;
            vweight.reserve(nelem * 2);
            double sum1 = std::accumulate(option.mElemWeights1.begin(), option.mElemWeights1.end(), 0.);
            double sum2 = std::accumulate(option.mElemWeights2.begin(), option.mElemWeights2.end(), 0.);
            for (int i = 0; i < nelem; i++) {
                vweight.push_back((int)round(option.mElemWeights1[i] / sum1 * imax));
                vweight.push_back((int)round(option.mElemWeights2[i] / sum2 * imax));
            }
            ubvec[0] = {1.f + option.mImbalance1};
            ubvec[1] = {1.f + option.mImbalance2};
        } else {
            ncon = 1;
            vweight.reserve(nelem);
            double sum1 = std::accumulate(option.mElemWeights1.begin(), option.mElemWeights1.end(), 0.);
            for (int i = 0; i < nelem; i++) {
                vweight.push_back((int)round(option.mElemWeights1[i] / sum1 * imax));
            }
            ubvec[0] = {1.f + option.mImbalance1};
        }
        
        // metis options
        int metis_option[METIS_NOPTIONS];
        metisError(METIS_SetDefaultOptions(metis_option), "METIS_SetDefaultOptions");
        metis_option[METIS_OPTION_NCUTS] = option.mNPartition;
        metis_option[METIS_OPTION_UFACTOR] = round(option.mImbalance1 * 1000);
        if (option.mCommVol) {
            metis_option[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_VOL;
            metis_option[METIS_OPTION_CONTIG] = 1;
        }
        
        int objval;
        if (option.mCommVol) {
            metisError(METIS_PartGraphKway(&nelem, &ncon, xadj, adjncy, 
                vweight.data(), vsize.data(), NULL, &nproc, NULL, ubvec, 
                metis_option, &objval, elemToProc.data()), "METIS_PartGraphKway");
        } else {
            metisError(METIS_PartGraphRecursive(&nelem, &ncon, xadj, adjncy, 
                vweight.data(), vsize.data(), NULL, &nproc, NULL, ubvec, 
                metis_option, &objval, elemToProc.data()), "METIS_PartGraphRecursive");
        }
         
        // free memory 
        freeAdjacency(xadj, adjncy);
    }
    
    XMPI::bcastEigen(elemToProc);
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
    for (int i = 0; i < nelem + 1; i++) eptr(i) = i * 4;
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
    if (retval != METIS_OK) throw std::runtime_error("DualGraph::metisError || "
        "Error in metis function: " + func_name);
}

void DualGraph::check_idx_t() {
    if (IDXTYPEWIDTH != 32)
        throw std::runtime_error("DualGraph::check_idx_t || Incompatible METIS build ||"
            "Please re-install METIS with 32-bit configuration. ||"
            "Or edit METIS_ROOT in CMakeLists.txt to locate the 32-bit build.");
}
