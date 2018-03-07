// eigenc.h
// created by Kuangdai on 27-Mar-2016 
// typedef of eigen matrices and vectors

#pragma once

#include "eigen.h"
#include <vector>
#include <array>
#include "global.h"

// pointwise fields
typedef Eigen::Matrix<Real, Eigen::Dynamic, 1> RColX;
typedef Eigen::Matrix<Real, Eigen::Dynamic, 3> RMatX3;
typedef Eigen::Matrix<Complex, Eigen::Dynamic, 1> CColX;
typedef Eigen::Matrix<Complex, Eigen::Dynamic, 3> CMatX3;
typedef std::array<CMatX3, nPntElem> arPP_CMatX3; // source 
typedef std::vector<arPP_CMatX3> vec_arPP_CMatX3; // off-axis source 
typedef Eigen::Matrix<Real, 1, 3> RRow3;         // receiver
typedef Eigen::Matrix<Real, 1, 6> RRow6;         // receiver
typedef Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic> CMatXX; // mpi buffer
typedef Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> RMatXX;
typedef Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RMatXX_RM;
typedef Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> CMatXX_RM;

// elemental fields - flat 
typedef Eigen::Matrix<Complex, Eigen::Dynamic, nPntElem * 1> CMatXN;
typedef Eigen::Matrix<Complex, Eigen::Dynamic, nPntElem * 3> CMatXN3;
typedef Eigen::Matrix<Complex, Eigen::Dynamic, nPntElem * 6> CMatXN6;
typedef Eigen::Matrix<Complex, Eigen::Dynamic, nPntElem * 9> CMatXN9;
typedef Eigen::Matrix<Real, Eigen::Dynamic, nPntElem * 1> RMatXN;
typedef Eigen::Matrix<Real, Eigen::Dynamic, nPntElem * 3> RMatXN3;
typedef Eigen::Matrix<Real, Eigen::Dynamic, nPntElem * 6> RMatXN6;
typedef Eigen::Matrix<Real, Eigen::Dynamic, nPntElem * 9> RMatXN9;
typedef Eigen::Matrix<Real, Eigen::Dynamic, nPntElem * 4> RMatXN4;  // particle relabelling
typedef Eigen::Matrix<Real, 1, nPntElem> RRowN;

// elemental fields - structured
typedef Eigen::Matrix<Real, nPntEdge, nPntEdge, Eigen::RowMajor> RMatPP;
typedef Eigen::Matrix<Complex, nPntEdge, nPntEdge, Eigen::RowMajor> CMatPP;
typedef std::array<CMatPP, 3> ar3_CMatPP;
typedef std::array<CMatPP, 6> ar6_CMatPP;
typedef std::array<CMatPP, 9> ar9_CMatPP;
typedef std::vector<CMatPP> vec_CMatPP;
typedef std::vector<ar3_CMatPP> vec_ar3_CMatPP;
typedef std::vector<ar6_CMatPP> vec_ar6_CMatPP;
typedef std::vector<ar9_CMatPP> vec_ar9_CMatPP;
const ar3_CMatPP zero_ar3_CMatPP = {CMatPP::Zero(), CMatPP::Zero(), CMatPP::Zero()};
const ar6_CMatPP zero_ar6_CMatPP = {CMatPP::Zero(), CMatPP::Zero(), CMatPP::Zero(), 
                                    CMatPP::Zero(), CMatPP::Zero(), CMatPP::Zero()};
const ar9_CMatPP zero_ar9_CMatPP = {CMatPP::Zero(), CMatPP::Zero(), CMatPP::Zero(), 
                                    CMatPP::Zero(), CMatPP::Zero(), CMatPP::Zero(), 
                                    CMatPP::Zero(), CMatPP::Zero(), CMatPP::Zero()};
                                    

