// Connectivity.h
// created by Kuangdai on 18-Jun-2016 
// mesh connectivity

#pragma once

#include <eigenp.h>

struct DecomposeOption;
struct MessagingInfo; 

class Connectivity {
    
public:    
    // constructor using exodus connectivity 
    Connectivity(const IMatX4 &excon);
    
    // constructor of a subset
    Connectivity(const Connectivity &super, const IColX &mask);
    
    // size, nelem
    int size() const {return mGlobalQuadID.rows();};
    
    // domain decomposition
    void decompose(const DecomposeOption &option, 
        int &nlocalGLL, std::vector<IMatPP> &localElemToGLL, 
        MessagingInfo &msgInfo, IColX &procMask) const;
    
private:
    IColX mGlobalQuadID;
    IMatX4 mConnectivity;
    
private:
    // form element-to-gll mapping 
    void formElemToGLL(int &ngll, std::vector<IMatPP> &elemToGLL, 
        std::vector<IColX> &neighbours, int ncommon) const;
    static void get_shared_DOF_quad(const IRow4 &connectivity1, const IRow4 &connectivity2, 
        std::vector<IRow2> &map1, std::vector<IRow2> &map2, int ielem);
    static void common_nodes(const IRow4 &a, const IRow4 &b, 
        int &ncommon, int &aindex, int &bindex);
    static void formNodeEdge();
    static std::array<std::vector<IRow2>, 4> sNodeIJPol;
    static std::array<std::vector<IRow2>, 4> sEdgeIJPol;
    static bool onEdge(int ipol, int jpol);
};


