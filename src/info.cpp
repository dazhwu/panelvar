#include "info.h"

model_options::model_options() {
    steps = 2;
    constant = false;
    level = true;
    beginner = false;
    timedumm = false;
    collapse = false;
    mmsc = "bic";
    transformation = "fd";
    irf = "girf";
}

basic_info::basic_info() {}

basic_info::basic_info(int _N, int _T, int _num_obs, int _num_instr, int _num_dep, int _num_dep_lags, int _num_indep,
                       int _diff_width, int _x_height, int _max_obs, int _min_obs, double _avg_obs,
                       vector<regular_variable> &dep_indep, vector<string> &_identifiers, model_options &options) {
    N = _N;
    T = _T;
    num_obs = _num_obs;
    num_instr = _num_instr;
    num_indep = _num_indep;
    num_dep = _num_dep;
    num_dep_lags = _num_dep_lags;
    diff_width = _diff_width;
    x_height = _x_height;
    max_obs = _max_obs;
    min_obs = _min_obs;
    avg_obs = _avg_obs;
    mmsc_lu = MMSC_LU();
    for (std::size_t i = 0, max = _identifiers.size(); i < max; ++i) {
        identifiers.push_back(_identifiers[i]);
    }

    for (int i = 0; i < num_dep; ++i)
        dep.push_back(dep_indep[i].name);

    for (std::size_t i = num_dep, max = dep_indep.size(); i < max; ++i) {
        int lag = dep_indep[i].lag;
        if (lag == 0)
            indep.push_back(dep_indep[i].name);
        else
            indep.push_back("L" + std::to_string(lag) + "." + dep_indep[i].name);
    }

    // for (std::size_t i = 0, max = dep_indep.size(); i < max; ++i)
    // {
    // 	if (i < num_dep)
    // 		dep.push_back(dep_indep[i].name);

    // 	else
    // 	{
    // 		int lag = dep_indep[i].lag;
    // 		if (lag == 0)
    // 			indep.push_back(dep_indep[i].name);
    // 		else
    // 			indep.push_back("L" + std::to_string(lag) + "." +
    // dep_indep[i].name);
    // 	}
    // }

    if (options.level && options.constant)
        indep.push_back("_con");
}

Regression::Regression() {}

Regression::Regression(RowMatrixXd &_beta, RowMatrixXd &_vcov, RowMatrixXd &_std_err, RowMatrixXd &w,
                       Ref<RowMatrixXd> _z, Ref<RowMatrixXd> _Cy, Ref<RowMatrixXd> _Cx, RowMatrixXd &_residual) {

    beta = _beta;
    vcov = _vcov;
    std_err = _std_err;
    Z_values = RowMatrixXd::Zero(beta.rows(), beta.cols());
    P_values = RowMatrixXd::Zero(beta.rows(), beta.cols());
    Residual = _residual;
    for (int i = 0; i < beta.rows(); ++i) {
        for (int j = 0; j < beta.cols(); ++j) {
            Z_values(i, j) = beta(i, j) / std_err(i, j);
            P_values(i, j) = 2 * (1 - standard_normalCDF(abs(Z_values(i, j))));
        }
    }

    weighting = w;
    Instruments = _z;
    Cy = _Cy;
    Cx = _Cx;
}

df_info::df_info() {}

df_info::df_info(int _N, int _T, int _f_d_i, int _f_lev_i, int _l_d_i, int _l_lev_i, int _max_lag) {
    N = _N;
    T = _T;
    first_diff_index = _f_d_i;
    first_level_index = _f_lev_i;
    last_diff_index = _l_d_i;
    last_level_index = _l_lev_i;
    max_lag = _max_lag;
}

MMSC_LU::MMSC_LU() {
    BIC = 0;
    HQIC = 0;
    AIC = 0;
}

MMSC_LU::MMSC_LU(double bic, double hqic, double aic) {
    BIC = bic;
    HQIC = hqic;
    AIC = aic;
}

z_info::z_info() {
    diff_height = 0;
    diff_width = 0;
    level_width = 0;
    level_height = 0;
    z_height = 0;
    z_width = 0;
    num_Dgmm_instr = 0;
    num_Lgmm_instr = 0;
    num_instr = 0;
}
z_info::z_info(int _diff_height, int _diff_width, int _level_width, int _level_height, int _z_height, int _z_width,
               int _num_Dgmm_instr, int _num_Lgmm_instr) {
    diff_height = _diff_height;
    diff_width = _diff_width;
    level_width = _level_width;
    level_height = _level_height;
    z_height = _z_height;
    z_width = _z_width;
    num_Dgmm_instr = _num_Dgmm_instr;
    num_Lgmm_instr = _num_Lgmm_instr;
    num_instr = _z_height;
}
