#pragma once
#include <Eigen/MatrixFunctions>
#include "pvar.h"

RowMatrixXd PAR1_matrix(RowMatrixXd &beta, int lags);
vector<RowMatrixXd> IRF_matrix(RowMatrixXd &par1, int ahead, int num_dep);

RowMatrixXd vcov_residual(RowMatrixXd &residual, VectorXi &is_na);
vector<RowMatrixXd>irf(string method, RowMatrixXd &residual, VectorXi &is_na, RowMatrixXd &beta, int ahead, int lags);
vector<RowMatrixXd>oirf(RowMatrixXd &, VectorXi &, RowMatrixXd &beta, int ahead, int lags);
vector<RowMatrixXd>girf(RowMatrixXd &, VectorXi &, RowMatrixXd &beta, int ahead, int lags);

