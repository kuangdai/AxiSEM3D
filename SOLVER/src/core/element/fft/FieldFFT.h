// FieldFFT.h
// created by Kuangdai on 19-May-2017 
// FFT between Fourier and physical spaces

#pragma once

#include "eigenc.h"

class FieldFFT {
public:
    
    static void transformF2P_3(const vec_ar3_CMatPP &uc, int Nr);
    static void transformF2P_6(const vec_ar6_CMatPP &uc, int Nr);
    static void transformF2P_9(const vec_ar9_CMatPP &uc, int Nr);
    
    static void transformP2F_3(vec_ar3_CMatPP &uc, int Nr);
    static void transformP2F_6(vec_ar6_CMatPP &uc, int Nr);
    static void transformP2F_9(vec_ar9_CMatPP &uc, int Nr);
    
    template<class vec_arY_CMatPP, class CMatXNY>
    static void makeFlat(const vec_arY_CMatPP &ucStruct, CMatXNY &ucFlat, int Nu) {
        for (int alpha = 0; alpha <= Nu; alpha++) {
            for (int i = 0; i < ucStruct[0].size(); i++) {
                for (int j = 0; j < nPntEdge; j++) {
                    ucFlat.block(alpha, nPE * i + nPntEdge * j, 1, nPntEdge) 
                        = ucStruct[alpha][i].block(j, 0, 1, nPntEdge);
                }
            } 
        } 
    };
    
    template<class vec_arY_CMatPP, class CMatXNY>
    static void makeStruct(vec_arY_CMatPP &ucStruct, const CMatXNY &ucFlat, int Nu) {
        for (int alpha = 0; alpha <= Nu; alpha++) {
            for (int i = 0; i < ucStruct[0].size(); i++) {
                for (int j = 0; j < nPntEdge; j++) {
                    ucStruct[alpha][i].block(j, 0, 1, nPntEdge)
                        = ucFlat.block(alpha, nPE * i + nPntEdge * j, 1, nPntEdge);
                }
            } 
        } 
    };
};

