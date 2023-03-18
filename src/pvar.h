#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <math.h>
#include <iostream>
#include <math.h>
#include <vector>
#include <string>
#include <tuple>
#include <cstring>
#include <fstream>

// #define EIGEN_USE_MKL_ALL
// #define EIGEN_INITIALIZE_MATRICES_BY_ZERO
// #define EIGEN_USE_BLAS
// #define EIGEN_USE_LAPACKE

#include <Eigen/Core>
#include <Eigen/StdVector>
#include <Eigen/QR>
#include <Eigen/LU>

using std::cout;
using std::endl;
using std::string;
using std::tuple;
using std::vector;

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::Ref;
using Eigen::VectorXcd;
using Eigen::VectorXi;

extern int DEP_lags;

using RowMatrixXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

// extern RowMatrixXd records;

#pragma omp declare reduction(+                             \
                              : RowMatrixXd                 \
                              : omp_out = omp_out + omp_in) \
    initializer(omp_priv = RowMatrixXd::Zero(omp_orig.rows(), omp_orig.cols()))
