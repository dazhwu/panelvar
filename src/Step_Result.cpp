
#include "Step_Result.h"
// Step_Result(W, W_inv, W_next, _XZ_W, M,  _M_XZ_W, zs, qs, ZuuZ,
//                   _zs_list, beta,  residual, _residual_t);

Step_Result::Step_Result(RowMatrixXd &W_, RowMatrixXd &W_inv_, RowMatrixXd &W_next_, RowMatrixXd &XZ_W_,
                         RowMatrixXd &M_, RowMatrixXd &_M_XZ_W_, RowMatrixXd &zs_, RowMatrixXd &qs_, RowMatrixXd &ZuuZ_,
                         RowMatrixXd &_zs_list_, RowMatrixXd &beta_, RowMatrixXd &residual_,
                         RowMatrixXd &_residual_t_) {

    W = W_;

    W_inv = W_inv_;
    W_next = W_next_;
    XZ_W = XZ_W_;
    M = M_;
    _M_XZ_W = _M_XZ_W_;
    zs = zs_;
    qs = qs_;
    ZuuZ = ZuuZ_;

    _zs_list = _zs_list_;
    beta = beta_;

    residual = residual_;
    _residual_t = _residual_t_;
}

void
Step_Result::add_vcov(RowMatrixXd &_vcov) {
    vcov = _vcov;
    std_err = RowMatrixXd::Zero(beta.rows(), beta.cols());
    RowMatrixXd temp;

    temp = vcov.diagonal();

    for (int i = 0; i < beta.rows(); ++i)
        for (int j = 0; j < beta.cols(); ++j)
            std_err(i, j) = sqrt(temp(i * beta.cols() + j, 0));
}