// eigenp.h
// created by Kuangdai on 3-May-2016 
// typedef of eigen matrices and vectors for preprocessor

#pragma once

#include "eigen.h"
#include <vector>
#include <array>
#include "global.h"

// 2D coordinates
typedef Eigen::Matrix<double, 2, 1> RDCol2;
typedef Eigen::Matrix<double, 2, 2> RDMat22;
typedef Eigen::Matrix<double, 2, 4> RDMat24;
typedef Eigen::Matrix<double, 1, 4> RDRow4;

// 3D coordinates
typedef Eigen::Matrix<double, 3, 3> RDMat33;
typedef Eigen::Matrix<double, 3, 1> RDCol3;
typedef Eigen::Matrix<double, Eigen::Dynamic, 3> RDMatX3;

// elemental fields 
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> RDColX;
typedef Eigen::Matrix<ComplexD, Eigen::Dynamic, 1> CDColX;
typedef Eigen::Matrix<double, Eigen::Dynamic, 4> RDMatX4;
typedef std::array<RDColX, nPntElem> arPP_RDColX;
// spectral 
typedef Eigen::Matrix<double, nPntEdge, nPntEdge, Eigen::RowMajor> RDMatPP;
typedef Eigen::Matrix<ComplexD, nPntEdge, nPntEdge, Eigen::RowMajor> CDMatPP;
typedef std::vector<CDMatPP> vec_CDMatPP;
typedef std::array<CDMatPP, 3> ar3_CDMatPP;
typedef std::vector<ar3_CDMatPP> vec_ar3_CDMatPP;
// cardinal
typedef Eigen::Matrix<double, 1, nPntElem> RDRowN;
typedef Eigen::Matrix<double, Eigen::Dynamic, nPntElem * 1> RDMatXN;

// local to global mapping
typedef Eigen::Matrix<int, nPntEdge, nPntEdge> IMatPP; 
typedef Eigen::Matrix<int, Eigen::Dynamic, 1> IColX;
typedef Eigen::Matrix<int, Eigen::Dynamic, 4> IMatX4;  

// others
typedef Eigen::Matrix<double, nPntEdge, 1> RDColP; // spectral constants
typedef Eigen::Matrix<double, Eigen::Dynamic, nPntElem * 4> RDMatXN4;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> RDMatXX;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RDMatXX_RM;
typedef Eigen::Matrix<ComplexD, Eigen::Dynamic, Eigen::Dynamic> CDMatXX;

// connectivity
typedef Eigen::Matrix<int, 1, 4> IRow4;
typedef Eigen::Matrix<int, 1, 2> IRow2; // ipol, jpol
typedef Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> IMatXX;


