#pragma once
#include <unordered_map>

#include "Common_Functions.h"
#include "IRF.h"
#include "pvar.h"
#include <Eigen/Eigenvalues>
#include <boost/math/distributions/chi_squared.hpp>

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

Hansen_test_info hansen_overid(const Ref<const RowMatrixXd> &, const Ref<const RowMatrixXd> &, int, int, int);
