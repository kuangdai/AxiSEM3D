// Elastic3D.h
// created by Kuangdai on 29-Apr-2016 
// base class of 3D elasticity 

#include "Elastic3D.h"
#include "Attenuation3D.h"

Elastic3D::Elastic3D(Attenuation3D *att):
mAttenuation(att) {
    // nothing
}

Elastic3D::~Elastic3D() {
    if (mAttenuation) delete mAttenuation;
} 

void Elastic3D::checkCompatibility(int Nr, bool isVoigt) const {
    if (mAttenuation) mAttenuation->checkCompatibility(Nr);
}

void Elastic3D::resetZero() {
    if (mAttenuation) mAttenuation->resetZero();
}

// data structure convertors
void Elastic3D::flattenVector(const vec_ar9_CMatPP &mat, CMatXN9 &row, int Nu) {
    for (int alpha = 0; alpha <= Nu; alpha++) 
        for (int i = 0; i < 9; i++) 
            for (int j = 0; j < nPntEdge; j++)
                row.block(alpha, nPE * i + nPntEdge * j, 1, nPntEdge) 
                    = mat[alpha][i].block(j, 0, 1, nPntEdge);
}

void Elastic3D::flattenVectorVoigt(const vec_ar9_CMatPP &mat, CMatXN6 &row, int Nu) {
    for (int alpha = 0; alpha <= Nu; alpha++) 
        for (int i = 0; i < 6; i++) 
            for (int j = 0; j < nPntEdge; j++)
                row.block(alpha, nPE * i + nPntEdge * j, 1, nPntEdge) 
                    = mat[alpha][i].block(j, 0, 1, nPntEdge);
}

void Elastic3D::stackupVector(const CMatXN9 &row, vec_ar9_CMatPP &mat, int Nu) {
    for (int alpha = 0; alpha <= Nu; alpha++) 
        for (int i = 0; i < 9; i++) 
            for (int j = 0; j < nPntEdge; j++)
                mat[alpha][i].block(j, 0, 1, nPntEdge) 
                    = row.block(alpha, nPE * i + nPntEdge * j, 1, nPntEdge);
}

void Elastic3D::stackupVectorVoigt(const CMatXN6 &row, vec_ar9_CMatPP &mat, int Nu) {
    for (int alpha = 0; alpha <= Nu; alpha++) 
        for (int i = 0; i < 6; i++) 
            for (int j = 0; j < nPntEdge; j++)
                mat[alpha][i].block(j, 0, 1, nPntEdge) 
                    = row.block(alpha, nPE * i + nPntEdge * j, 1, nPntEdge);
}
