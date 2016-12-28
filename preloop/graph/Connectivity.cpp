// Connectivity.cpp
// created by Kuangdai on 18-Jun-2016 
// mesh connectivity

#include "Connectivity.h"
#include "DualGraph.h"
#include "XMPI.h"
#include <algorithm>
#include <map>

#include <XTimer.h>

std::array<std::vector<IRow2>, 4> Connectivity::sNodeIJPol;
std::array<std::vector<IRow2>, 4> Connectivity::sEdgeIJPol;

Connectivity::Connectivity(const std::vector<std::array<int, 4>> &excon) {
    int nelem = excon.size();
    mGlobalQuadID = IColX(nelem);
    mConnectivity = IMatX4(nelem, 4);
    for (int i = 0; i < nelem; i++) {
        mGlobalQuadID(i) = i;
        mConnectivity(i, 0) = excon[i][0];
        mConnectivity(i, 1) = excon[i][1];
        mConnectivity(i, 2) = excon[i][2];
        mConnectivity(i, 3) = excon[i][3];
    }
}

Connectivity::Connectivity(const Connectivity &super, const IColX &mask) {
    int nsup = super.size();
    int nsub = 0;
    for (int i = 0; i < nsup; i++) nsub += (mask(i) ? 1 : 0);
    mGlobalQuadID = IColX(nsub);
    mConnectivity = IMatX4(nsub, 4);
    int pos = 0;
    for (int i = 0; i < nsup; i++) {
        if (mask(i)) {
            mGlobalQuadID(pos) = super.mGlobalQuadID(i);
            mConnectivity.row(pos) = super.mConnectivity.row(i);
            pos++;
        }
    }
}

void Connectivity::formElemToGLL(int &ngll, std::vector<IMatPP> &elemToGLL, std::vector<IColX> &neighbours) const {
    // form static
    if (sNodeIJPol[0].size() == 0) formNodeEdge();
    
    // get neighbours from metis
    DualGraph::formNeighbourhood(mConnectivity, 1, neighbours);

    // build local to global mapping
    ngll = 0;
    int nelem = size();
    elemToGLL = std::vector<IMatPP>(nelem, IMatPP::Constant(-1));

    // loop over all elements
    for (int ielem = 0; ielem < nelem; ielem++) {

        // the mask will be True for all points, that are not shared with an
        // element with lower element id
        IMatPP mask = IMatPP::Ones();

        // loop over all neighbouring elements with smaller element id
        for (int in = 0; in < neighbours[ielem].rows(); in++) {
            int ineighbour = neighbours[ielem](in);
            if (ineighbour >= ielem) continue;
            
            // get ipol and jpol of shared points in both elements
            std::vector<IRow2> map1, map2;
            get_shared_DOF_quad(mConnectivity.row(ineighbour), mConnectivity.row(ielem), map1, map2, ielem);
            
            // for global points that already have a number, get the global
            // number from the other element
            for (int i = 0; i < map2.size(); i++) {
                mask(map2[i](0), map2[i](1)) = 0;
                elemToGLL[ielem](map2[i](0), map2[i](1)) = 
                    elemToGLL[ineighbour](map1[i](0), map1[i](1));
            }
        }

        // compute number of new global points in this element and label them
        for (int ipol = 0; ipol <= nPol; ipol++) {
            for (int jpol = 0; jpol <= nPol; jpol++) {
                if (mask(ipol, jpol)) elemToGLL[ielem](ipol, jpol) = ngll++;
            }
        }        
    }
}

void Connectivity::decompose(const DecomposeOption &option, 
    int &nGllLocal, std::vector<IMatPP> &elemToGllLocal, 
    MessagingInfo &msg, IColX &procMask) const {
    // domain decomposition
    XTimer::begin("Metis Part", 3);
    IColX elemToProc;
    DualGraph::decompose(mConnectivity, option, XMPI::nproc(), elemToProc);
    XTimer::end("Metis Part", 3);
    
    // global element-gll mapping
    // neighbourhood with ncommon = 1 (NOT 2)
    // NOTE: Though we decompose with ncommon = 2, metis may still (but rarely) yields 
    //       a decomposition where two processors only share one single point. 
    //       This could happen when a large nproc is used on a relatively small mesh.
    XTimer::begin("global element-gll", 3);
    int nElemGlobal = size();
    int nGllGlobal = 0;
    std::vector<IMatPP> elemToGllGlobal;
    std::vector<IColX> neighbours;
    formElemToGLL(nGllGlobal, elemToGllGlobal, neighbours);
    XTimer::end("global element-gll", 3);
    
    // map of to-be-communicated global gll points 
    XTimer::begin("to-be-communicated global", 3);
    // key: proc_id
    // value: map<global_gll, array_of_3(elem_id, ipol, jpol)>
    std::map<int, std::map<int, std::array<int, 3>>> gllCommGlb;
    for (int ielem = 0; ielem < nElemGlobal; ielem++) {
        if (elemToProc(ielem) != XMPI::rank()) continue;
        for (int in = 0; in < neighbours[ielem].rows(); in++) {
            int ineighbour = neighbours[ielem](in);
            // within the same proc
            if (elemToProc(ineighbour) == XMPI::rank()) continue;
            // add proc 
            int rankOther = elemToProc(ineighbour);
            gllCommGlb.insert(std::pair<int, std::map<int, std::array<int, 3>>>
                (rankOther, std::map<int, std::array<int, 3>>()));
                
            // create temp vector of neighbour gll's for fast search
            std::vector<int> gllOther;
            for (int ipol = 0; ipol <= nPol; ipol++) 
                for (int jpol = 0; jpol <= nPol; jpol++) 
                    if (onEdge(ipol, jpol)) gllOther.push_back(elemToGllGlobal[ineighbour](ipol, jpol));
            // find common gll
            int nfound = 0;
            for (int ipol = 0; ipol <= nPol; ipol++) {
                for (int jpol = 0; jpol <= nPol; jpol++) {
                    if (!onEdge(ipol, jpol)) continue;
                    int targetGll = elemToGllGlobal[ielem](ipol, jpol);
                    if (std::find(gllOther.begin(), gllOther.end(), targetGll) != gllOther.end()) {
                        std::array<int, 3> ielem_ipol_jpol;
                        ielem_ipol_jpol[0] = ielem;
                        ielem_ipol_jpol[1] = ipol;
                        ielem_ipol_jpol[2] = jpol;
                        // NOTE: the components are inherently sorted by global gll-point tags
                        //       we do not care how they are sorted by std::map since all processors
                        //       should use the same sorting rule 
                        gllCommGlb.at(rankOther).insert(std::pair<int, std::array<int, 3>>
                            (targetGll, ielem_ipol_jpol));
                        nfound++;
                    }
                }
            }
            // either a common edge or a common point
            if (nfound != nPntEdge && nfound != 1) 
                throw std::runtime_error("Connectivity::decompose || Domain decomposition failed.");
        }
    }
    XTimer::end("to-be-communicated global", 3);
    
    // global-to-local element map and local mask
    XTimer::begin("global-to-local element", 3);
    IColX elemGlbToLoc = IColX::Constant(nElemGlobal, -1);
    procMask = IColX::Zero(nElemGlobal);
    int nElemLocal = 0;
    for (int ielem = 0; ielem < nElemGlobal; ielem++) {
        if (elemToProc(ielem) == XMPI::rank()) {
            elemGlbToLoc(ielem) = nElemLocal++;
            procMask(ielem) = 1;
        }
    } 
    XTimer::end("global-to-local element", 3);
    
    // local element-gll mapping
    XTimer::begin("local element-gll", 3);
    Connectivity(*this, procMask).formElemToGLL(nGllLocal, elemToGllLocal, neighbours);
    XTimer::end("local element-gll", 3);
    
    // form local messaging
    XTimer::begin("local messaging", 3);
    msg.mIProcComm.clear();
    msg.mNLocalPoints.clear();
    msg.mILocalPoints.clear();
    for (auto it_proc = gllCommGlb.begin(); it_proc != gllCommGlb.end(); it_proc++) {
        msg.mIProcComm.push_back(it_proc->first);
        std::vector<int> gll_loc;
        for (auto it_gll = it_proc->second.begin(); it_gll != it_proc->second.end(); it_gll++) {
            std::array<int, 3> ielem_ipol_jpol = it_gll->second;
            int ielem_glb = ielem_ipol_jpol[0];
            int ipol = ielem_ipol_jpol[1];
            int jpol = ielem_ipol_jpol[2];
            int ielem_loc = elemGlbToLoc(ielem_glb);
            gll_loc.push_back(elemToGllLocal[ielem_loc](ipol, jpol));
        }
        msg.mNLocalPoints.push_back(gll_loc.size());
        msg.mILocalPoints.push_back(gll_loc);
    }
    msg.mNProcComm = msg.mIProcComm.size();
    msg.mReqSend = std::vector<XMPI::Request>(msg.mNProcComm, XMPI::Request());
    msg.mReqRecv = std::vector<XMPI::Request>(msg.mNProcComm, XMPI::Request());
    XTimer::end("local messaging", 3);
}

void Connectivity::get_shared_DOF_quad(const IRow4 &connectivity1, const IRow4 &connectivity2, 
    std::vector<IRow2> &map1, std::vector<IRow2> &map2, int ielem) {
    int ncommon = 0;
    int index1 = -1;
    int index2 = -1;
    common_nodes(connectivity1, connectivity2, ncommon, index1, index2);
    
    if (ncommon == 1) {
        // a single common point
        map1 = sNodeIJPol[index1];
        map2 = sNodeIJPol[index2];
    } else { 
        // a common edge
        map1 = sEdgeIJPol[index1];
        map2 = sEdgeIJPol[index2];
        std::reverse(map2.begin(), map2.end());
    }    
}

void Connectivity::common_nodes(const IRow4 &a, const IRow4 &b, 
    int &ncommon, int &aindex, int &bindex) {    
    IRow2 ac, bc;
    ncommon = 0;
    for (int i = 0; i < 4; i++) {
        int ai = a(i);
        for (int j = 0; j < 4; j++) {
            if (ai == b(j)) {
                ac(ncommon) = i;
                bc(ncommon) = j;
                ncommon++;
                break;
            }
        }
    }
    if (ncommon == 1) {
        aindex = ac(0);
        bindex = bc(0);
    } else if (ncommon == 2) {
        aindex = ac.minCoeff();
        bindex = bc.minCoeff();
        if (aindex == 0 && ac.maxCoeff() == 3) aindex = 3;
        if (bindex == 0 && bc.maxCoeff() == 3) bindex = 3;
    } else {
        throw std::runtime_error("Connectivity::common_nodes || Error topology in mesh.");
    }
}

void Connectivity::formNodeEdge() {
    // node to ipol jpol
    IRow2 ijpol;
    ijpol(0) = 0;
    ijpol(1) = 0;
    sNodeIJPol[0].push_back(ijpol);
    ijpol(0) = nPol;
    ijpol(1) = 0;
    sNodeIJPol[1].push_back(ijpol);
    ijpol(0) = nPol;
    ijpol(1) = nPol;
    sNodeIJPol[2].push_back(ijpol);
    ijpol(0) = 0;
    ijpol(1) = nPol;
    sNodeIJPol[3].push_back(ijpol);
    // edge to ipol jpol
    for (int i = 0; i < nPntEdge; i++) {
        ijpol(0) = i;
        ijpol(1) = 0;
        sEdgeIJPol[0].push_back(ijpol);
    }
    for (int i = 0; i < nPntEdge; i++) {
        ijpol(0) = nPol;
        ijpol(1) = i;
        sEdgeIJPol[1].push_back(ijpol);
    }
    for (int i = 0; i < nPntEdge; i++) {
        ijpol(0) = nPol - i;
        ijpol(1) = nPol;
        sEdgeIJPol[2].push_back(ijpol);
    }
    for (int i = 0; i < nPntEdge; i++) {
        ijpol(0) = 0;
        ijpol(1) = nPol - i;
        sEdgeIJPol[3].push_back(ijpol);
    }
}

bool Connectivity::onEdge(int ipol, int jpol) {
    if (ipol > 0 && ipol < nPol && jpol > 0 && jpol < nPol) return false;
    return true;
}




