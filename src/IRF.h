#pragma once
#include "pvar.h"
#include <Eigen/MatrixFunctions>

RowMatrixXd PAR1_matrix(RowMatrixXd &beta, int lags);
vector<RowMatrixXd> IRF_matrix(RowMatrixXd &par1, int ahead, int num_dep);

RowMatrixXd vcov_residual(RowMatrixXd &residual, VectorXi &is_na);
RowMatrixXd irf(string method, RowMatrixXd &residual, VectorXi &is_na, RowMatrixXd &beta, int ahead, int lags);
RowMatrixXd oirf(RowMatrixXd &, VectorXi &, RowMatrixXd &beta, int ahead, int lags);
RowMatrixXd girf(RowMatrixXd &, VectorXi &, RowMatrixXd &beta, int ahead, int lags);

RowMatrixXd calculate_irf_matrix(int, int, vector<RowMatrixXd> &);