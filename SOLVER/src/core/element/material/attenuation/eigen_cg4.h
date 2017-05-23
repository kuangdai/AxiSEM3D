// eigen_cg4.h
// created by Kuangdai on 27-Mar-2016 
// typedef of eigen matrices and vectors

#pragma once

#include "eigen.h"
#include <vector>
#include <array>
#include "global.h"

// elemental fields - flat 
static const int nCG = 4;
typedef Eigen::Matrix<Real, Eigen::Dynamic, nCG * 6> RMatX46;
typedef Eigen::Matrix<Real, Eigen::Dynamic, nCG> RMatX4;


// elemental fields - structured
typedef Eigen::Matrix<Real, 1, 4> RRow4; 
typedef Eigen::Matrix<Complex, 1, 4> CRow4;
typedef std::array<CRow4, 6> ar6_CRow4;
typedef std::vector<ar6_CRow4> vec_ar6_CRow4;
const ar6_CRow4 zero_ar6_CRow4 = {CRow4::Zero(), CRow4::Zero(), CRow4::Zero(), 
                                  CRow4::Zero(), CRow4::Zero(), CRow4::Zero()};

