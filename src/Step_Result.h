//
// Created by Tiger on 5/3/2022.
//

#ifndef UNTITLED_STEP_RESULT_H
#define UNTITLED_STEP_RESULT_H

#include <iostream>
#include <math.h>
#include <string>
#include <tuple>
#include <vector>

#include "Common_Functions.h"
#include "pvar.h"

#define EIGEN_INITIALIZE_MATRICES_BY_ZERO

using std::string;
using std::tuple;
using std::vector;

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::Ref;

using RowMatrixXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

class Step_Result {
  public:
    RowMatrixXd residual, _residual_t, XZ_W, W, W_inv, W_next, _M_XZ_W, zs, qs, ZuuZ, M, _zs_list, beta, vcov, std_err;

    Step_Result(RowMatrixXd &W_, RowMatrixXd &W_inv_, RowMatrixXd &W_next_, RowMatrixXd &XZ_W_, RowMatrixXd &M_,
                RowMatrixXd &_M_XZ_W_, RowMatrixXd &zs_, RowMatrixXd &qs_, RowMatrixXd &ZuuZ_, RowMatrixXd &_zs_list_,
                RowMatrixXd &beta_, RowMatrixXd &residual_, RowMatrixXd &_residual_t_);
    void add_vcov(RowMatrixXd &_vcov);
};

#endif   // UNTITLED_STEP_RESULT_H
