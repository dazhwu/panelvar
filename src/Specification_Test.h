#pragma once
#include <unordered_map>

#include "pvar.h"
#include "Common_Functions.h"
#include "IRF.h"
#include <boost/math/distributions/chi_squared.hpp>
#include <Eigen/Eigenvalues> 

// class AR_test_info {
// public:
//     int lag;
//     double AR;
//     double P_value;
//     AR_test_info();
//     AR_test_info(int _lag, double _AR, double _P);

// };

class Hansen_test_info {
public:
    double test_value;
    int df;
    double P_value;
    double critical_value;
    Hansen_test_info();
    Hansen_test_info(double, int, double, double);

};


VectorXcd stability_test(RowMatrixXd &beta, int lags);

Hansen_test_info hansen_overid(const Ref<const RowMatrixXd> &, const Ref<const RowMatrixXd> &, int, int,  int);

// vector<MatrixXd> AR_get_diff_XR(int N, 
//                                 std::unordered_map<string, RowMatrixXd> &final_xy_tables,
//                                 Eigen::Ref<RowMatrixXd> beta,
//                                 Eigen::Ref<RowMatrixXd> ori_residual,
//                                  int diff_width,
//                                 int r0_height, string transformation,
//                                 bool level); 

// vector<AR_test_info>
// AR_test(int , int , std::unordered_map<string, RowMatrixXd> &,Eigen::Ref<RowMatrixXd> ,
//         Eigen::Ref<RowMatrixXd> zs_list, Eigen::Ref<RowMatrixXd> ori_residual,
//         Eigen::Ref<RowMatrixXd> M_XZ_W, Eigen::Ref<RowMatrixXd> vcov,
//         Eigen::Ref<RowMatrixXd> beta,
//         int diff_width, string transformation,
//         bool level);