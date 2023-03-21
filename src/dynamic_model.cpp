#include "dynamic_model.h"

int N, T;

RowMatrixXd Dgmm_table, iv_table, Delta_iv_table, Lgmm_table;
RowMatrixXd ori_x_table, ori_y_table, ori_Diff_x_table, ori_Diff_y_table;

std::unordered_map<string, RowMatrixXd> xy_tables, final_xy_tables;

struct df_info
get_info(model_options options, model m) {
    int max_lag = 0;
    int max_Dgmm_minlag = 0;
    int max_Lgmm_minlag = 0;
    for (auto var : m.dep_indep) {
        if (var.lag > max_lag)
            max_lag = var.lag;
    }
    for (auto var : m.iv_vars) {
        if (var.lag > max_lag)
            max_lag = var.lag;
    }
    for (auto var : m.Dgmm_vars) {
        if (var.min_lag > max_Dgmm_minlag)
            max_Dgmm_minlag = var.min_lag;
    }
    for (auto var : m.Lgmm_vars) {
        if (var.min_lag > max_Lgmm_minlag)
            max_Lgmm_minlag = var.min_lag;
    }
    int last_level_index = T - 1;   // last period as 0 based
    int last_diff_index = T - 1;

    int first_level_index, first_diff_index;
    if (max_lag > max_Lgmm_minlag)
        first_level_index = max_lag;
    else
        first_level_index = max_Lgmm_minlag;
    if ((options.level) && (options.transformation == "fod"))
        first_diff_index = first_level_index;
    else
        first_diff_index = first_level_index + 1;   // # max(max_lag + 1, max_Dgmm_minlag)

    if (first_diff_index + 2 > last_diff_index)   // # to do: change 3 to something rated to AR(p)
        throw std::invalid_argument("Not enough periods to run the model");
    struct df_info tbr = df_info(N, T, first_diff_index, first_level_index, last_diff_index, last_level_index, max_lag);
    return tbr;
}

vector<returned_result>
prepare_data(vector<string> &identifiers, Ref<RowMatrixXd> pdata, vector<int> &dimensions, vector<string> &df_cols,
             vector<string> &col_timedumm, int ahead, int num_draws) {
    N = dimensions[0];
    T = dimensions[1];
    RowMatrixXd fd_data = get_first_diff_table(pdata, N);

    vector<model> the_models = generate_list_models(T);
    vector<returned_result> results;
    VectorXi na_records;
    int m_num = 1;
    for (std::size_t m = 0, max = the_models.size(); m < max; ++m) {
        xy_tables.clear();
        final_xy_tables.clear();
        struct df_info info;
        try {
            info = get_info(options, the_models[m]);
        } catch (const std::exception &e) {
            continue;
        }
        if (options.timedumm) {
            update_time_dummies(col_timedumm, info, the_models[m]);
        }
        handle_tables(pdata, fd_data, the_models[m], df_cols, info);
        RowMatrixXd z_table;
        struct z_info z_information;
        std::tie(z_table, z_information) = get_z_table(
            N, T, info, Dgmm_table, Lgmm_table, iv_table, Delta_iv_table, the_models[m].Dgmm_vars,
            the_models[m].Lgmm_vars, the_models[m].iv_vars, options.level, options.transformation, options.collapse);
        na_records = VectorXi::Zero(final_xy_tables["Cy"].rows());
        struct basic_info model_info =
            generate_model_info(N, T, z_table, na_records, z_information, the_models[m], identifiers);

        // //saveData("z.csv", z_table);
        returned_result model_result;
        temp_result tmp_res;
        try {
            std::tie(model_result, tmp_res) = regular_process(z_table, final_xy_tables, model_info, options);
            model_result.dep_indep = the_models[m].dep_indep;
            model_result.model_info.model_name = "m" + std::to_string(m_num);
            model_result.command_str = the_models[m].command_str;
            model_result.irf.reserve(3);   //[0] irf, [1] lower limit, [2] Upper limit
            process_residual(model_result, na_records);
            RowMatrixXd irf_result =
                irf(options.irf, model_result.regression_result.Residual, na_records,
                    model_result.regression_result.beta, ahead, model_result.model_info.num_dep_lags);

            RowMatrixXd upper, lower;
            std::tie(lower, upper) =
                Bootstrapping(options.irf, num_draws, ahead, z_table, na_records, tmp_res, model_info, options);
            model_result.irf.push_back(irf_result);
            model_result.irf.push_back(lower);
            model_result.irf.push_back(upper);
            results.push_back(model_result);
            m_num += 1;
        } catch (const std::exception &e) {
            continue;
        }
    }
    return results;
}
void
process_residual(returned_result &m, VectorXi &na_records) {

    int num_rows = m.regression_result.Residual.rows();
    int num_cols = m.regression_result.Residual.cols();

    RowMatrixXd NA_Row = RowMatrixXd::Constant(1, num_cols, NAN);
    for (int i = 0; i < num_rows; ++i) {
        if (na_records[i] == 1)
            m.regression_result.Residual.row(i) = NA_Row;
    }
}
void
update_time_dummies(vector<string> &col_timedumm, struct df_info &info, model &m) {
    for (int i = info.first_diff_index; i <= info.last_diff_index; ++i) {
        string var_name = col_timedumm[i];
        m.dep_indep.push_back(regular_variable(var_name, 0));
        m.iv_vars.push_back(regular_variable(var_name, 0));
    }
}

struct basic_info
generate_model_info(int N, int T, Ref<RowMatrixXd> z_table, VectorXi &na_records, struct z_info &z_information,
                    model &m, vector<string> &identifiers) {
    int num_obs, max_obs, min_obs;
    double avg_obs;
    std::tie(num_obs, max_obs, min_obs, avg_obs) =
        prepare_reg(z_table, final_xy_tables["Cx"], final_xy_tables["Cy"], na_records, z_information,
                    options.transformation, options.level);
    int num_indep = final_xy_tables["Cx"].cols();
    int num_dep = m.num_dep;   // final_xy_tables["Cy"].cols();

    int num_dep_lags = m.num_dep_lags;

    struct basic_info model_info(N, T, num_obs, z_information.num_instr, num_dep, num_dep_lags, num_indep,
                                 z_information.diff_width, z_information.z_width, max_obs, min_obs, avg_obs,
                                 m.dep_indep, identifiers, options);

    return model_info;
}

void
handle_tables(Ref<RowMatrixXd> pdata, Ref<RowMatrixXd> fd_data, model &m, vector<string> &df_cols, df_info info) {
    get_gmm_tables(pdata, fd_data, m.Dgmm_vars, m.iv_vars, m.Lgmm_vars, df_cols, options.level);
    get_xy_table_dict(pdata, m, df_cols, options.transformation);
    get_final_xy_tables(info, options.level, options.transformation);
}

void
get_gmm_tables(Ref<RowMatrixXd> pdata, Ref<RowMatrixXd> fd_data, vector<gmm_var> &Dgmm_vars,
               vector<regular_variable> &iv_vars, vector<gmm_var> &Lgmm_vars, vector<string> &df_cols, bool level) {
    Dgmm_table = gen_table(pdata, Dgmm_vars, df_cols);
    iv_table = gen_table(pdata, iv_vars, df_cols);
    Delta_iv_table = gen_table(fd_data, iv_vars, df_cols);
    if (level)
        Lgmm_table = gen_table(fd_data, Lgmm_vars, df_cols);
}

template <class ty>
RowMatrixXd
gen_table(Ref<RowMatrixXd> ori_data, vector<ty> variable_list, vector<string> df_cols) {
    int num_variables = variable_list.size();

    vector<int> which_col;
    for (auto var : variable_list) {
        int the_index = getIndex(df_cols, var.name);
        which_col.push_back(the_index);
    }

    int height = T;   // end_row - start_row + 1
    RowMatrixXd tbr = RowMatrixXd::Constant(height * N, num_variables, NAN);

#pragma omp parallel for
    for (int j = 0; j < num_variables; ++j) {
        int lag = variable_list[j].lag;
        if (lag == 0)
            tbr.col(j) = ori_data.col(which_col[j]);
        else {
#pragma omp parallel for
            for (int i = 0; i < N; ++i) {
                // int col = 0;
                Ref<RowMatrixXd> ori_i = ori_data.block(i * T, which_col[j], T, 1);
                Ref<RowMatrixXd> tbr_i = tbr.block(i * height, j, height, 1);
                tbr_i.block(0, 0, lag, 1) = RowMatrixXd::Constant(lag, 1, NAN);
                for (int k = lag; k < height; ++k) {
                    if (!isnan(ori_i(k, 0)))
                        tbr_i(k, 0) = ori_i(k - lag, 0);
                }
            }
        }
    }
    return tbr;
}

void
get_xy_table_dict(Ref<RowMatrixXd> pdata, model &m, vector<string> &df_cols, string transformation) {
    int num_var = m.dep_indep.size();

    vector<regular_variable> dep{&m.dep_indep[0], &m.dep_indep[m.num_dep]};
    vector<regular_variable> indeps{&m.dep_indep[m.num_dep], &m.dep_indep[num_var]};

    RowMatrixXd ori_y_table = gen_table(pdata, dep, df_cols);
    RowMatrixXd ori_x_table = gen_table(pdata, indeps, df_cols);

    xy_tables["x"] = ori_x_table;
    xy_tables["y"] = ori_y_table;
    RowMatrixXd ori_Diff_y_table = get_first_diff_table(ori_y_table, N);
    RowMatrixXd ori_Diff_x_table = get_first_diff_table(ori_x_table, N);
    if (transformation == "fd") {
        xy_tables["Dy"] = ori_Diff_y_table;
        xy_tables["Dx"] = ori_Diff_x_table;
    } else {   // fod

        xy_tables["Dy"] = get_fod_table(ori_y_table, ori_y_table, N, 0, m.num_dep, m.num_dep_lags);
        xy_tables["Dx"] = get_fod_table(ori_x_table, ori_y_table, N, 1, m.num_dep, m.num_dep_lags);

        xy_tables["Diff_y"] = ori_Diff_y_table;
        xy_tables["Diff_x"] = ori_Diff_x_table;
    }
}

void
get_final_xy_tables(df_info info, bool level, string transformation) {
    vector<int> Dcut = {info.first_diff_index, info.last_diff_index};
    vector<int> cut = {info.first_level_index, info.last_level_index};

    int Dcut_height = Dcut[1] - Dcut[0] + 1;
    int cut_height = cut[1] - cut[0] + 1;

    RowMatrixXd Cy, Cx, Diff_y, Diff_x;
    if (level) {
        std::tie(Cy, Cx) = get_final_xy_systemGMM(xy_tables["Dx"], xy_tables["Dy"], xy_tables["x"], xy_tables["y"], cut,
                                                  Dcut, cut_height, Dcut_height, transformation);
    } else {
        std::tie(Cy, Cx) = get_final_xy_diffGMM(xy_tables["Dy"], xy_tables["Dx"], Dcut, Dcut_height);
    }
    final_xy_tables["Cy"] = Cy;
    final_xy_tables["Cx"] = Cx;
    if (transformation == "fod") {
        std::tie(Diff_y, Diff_x) = get_final_xy_diffGMM(xy_tables["Diff_y"], xy_tables["Diff_x"], Dcut, Dcut_height);
        if (level) {
            int height = Diff_y.rows();
            RowMatrixXd zeros = RowMatrixXd::Zero(height, 1);
            Diff_x.conservativeResize(height, Diff_x.cols() + 1);
            Diff_x.col(Diff_x.cols() - 1) = zeros;
        }
        final_xy_tables["Diff_y"] = Diff_y;
        final_xy_tables["Diff_x"] = Diff_x;
    }
}

std::tuple<RowMatrixXd, RowMatrixXd>
get_final_xy_systemGMM(Ref<RowMatrixXd> Dx, Ref<RowMatrixXd> Dy, Ref<RowMatrixXd> x, Ref<RowMatrixXd> y,
                       vector<int> cut, vector<int> Dcut, int cut_height, int Dcut_height, string transformation) {
    int height_total = Dcut_height + cut_height;

    int width;   //= x.cols();
    int num_dep = y.cols();
    int Dx_height = Dx.rows() / N;
    int x_height = x.rows() / N;
    if (options.constant)
        width = x.cols() + 1;
    else
        width = x.cols();
    RowMatrixXd Cx(height_total * N, width);
    RowMatrixXd Cy(height_total * N, y.cols());
#pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        Ref<RowMatrixXd> Dy_i = Dy.block(i * Dx_height, 0, Dx_height, num_dep);
        Ref<RowMatrixXd> y_i = y.block(i * x_height, 0, x_height, num_dep);

        Ref<RowMatrixXd> temp_y = Cy.block(height_total * i, 0, height_total, num_dep);
        temp_y.block(0, 0, Dcut_height, num_dep) = Dy_i.block(Dcut[0], 0, Dcut_height, num_dep);
        temp_y.block(Dcut_height, 0, height_total - Dcut_height, num_dep) =
            y_i.block(cut[0], 0, height_total - Dcut_height, num_dep);
        Ref<RowMatrixXd> Dx_i = Dx.middleRows(i * Dx_height, Dx_height);
        Ref<RowMatrixXd> x_i = x.middleRows(i * x_height, x_height);
        Ref<RowMatrixXd> temp_x = Cx.middleRows(height_total * i, height_total);
        temp_x.block(0, 0, Dcut_height, x_i.cols()) = Dx_i.block(Dcut[0], 0, Dcut_height, x_i.cols());
        ;
        temp_x.block(Dcut_height, 0, height_total - Dcut_height, x_i.cols()) =
            x_i.block(cut[0], 0, height_total - Dcut_height, x_i.cols());
        if (options.constant) {
            temp_x.block(0, width - 1, Dcut_height, 1) = RowMatrixXd::Zero(Dcut_height, 1);
            temp_x.block(Dcut_height, width - 1, height_total - Dcut_height, 1) =
                RowMatrixXd::Ones(height_total - Dcut_height, 1);
        }
        if (transformation == "fod") {
            temp_y.row(0) = RowMatrixXd::Constant(1, num_dep, NAN);
        }
    }
    //   saveData("before_x.csv", Cx);
    //  saveData("before_y.csv", Cy);
    return std::make_tuple(Cy, Cx);
}

std::tuple<RowMatrixXd, RowMatrixXd>
get_final_xy_diffGMM(Ref<RowMatrixXd> Dy, Ref<RowMatrixXd> Dx, vector<int> Dcut, int Dcut_height) {
    int Dx_height = Dx.rows() / N;

    int width = Dx.cols();
    int num_dep = Dy.cols();

    RowMatrixXd Cx(Dcut_height * N, width);
    RowMatrixXd Cy(Dcut_height * N, num_dep);
#pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        Ref<RowMatrixXd> Dy_i = Dy.block(i * Dx_height, 0, Dx_height, num_dep);
        Ref<RowMatrixXd> temp_y = Cy.block(Dcut_height * i, 0, Dcut_height, num_dep);
        temp_y.block(0, 0, Dcut_height, num_dep) = Dy_i.block(Dcut[0], 0, Dcut_height, num_dep);
        Ref<RowMatrixXd> Dx_i = Dx.block(i * Dx_height, 0, Dx_height, width);

        Ref<RowMatrixXd> temp_x = Cx.block(Dcut_height * i, 0, Dcut_height, width);
        temp_x.block(0, 0, Dcut_height, width) = Dx_i.block(Dcut[0], 0, Dcut_height, width);
    }
    return std::make_tuple(Cy, Cx);
}

void
prepare_reg_fod(Ref<RowMatrixXd> Diff_x, Ref<RowMatrixXd> Diff_y) {

    vector<bool> row_if_nan_y = row_has_nan(Diff_y);
    vector<bool> row_if_nan_x = row_has_nan(Diff_x);

#pragma omp parallel for
    for (int i = 0; i < Diff_y.rows(); ++i) {
        if ((row_if_nan_y[i]) || (row_if_nan_x[i])) {
            Diff_x.row(i) = RowMatrixXd::Zero(1, Diff_x.cols());
            Diff_y.row(i) = RowMatrixXd::Zero(1, Diff_y.cols());
        }
    }
}

std::tuple<int, int, int, double>
prepare_reg(Ref<RowMatrixXd> z_list, Ref<RowMatrixXd> Cx, Ref<RowMatrixXd> Cy, VectorXi &na_records,
            z_info z_information, string transformation, bool level) {
    int xy_height = Cy.rows() / N;
    int z_height = z_list.rows() / N;
    int Cx_width = Cx.cols();
    int Cy_width = Cy.cols();
    int z_width = z_list.cols();

    int num_NA = 0;
    // int total = 0;
    int max = 0;
    int min = 0;
    // //saveData("before_x.csv", final_xy_tables["Diff_x"]);
    if (transformation == "fod")
        prepare_reg_fod(final_xy_tables["Diff_x"], final_xy_tables["Diff_y"]);

    // #pragma omp parallel for reduction(+ : num_NA)

    for (int i = 0; i < N; ++i) {
        Ref<RowMatrixXd> x = Cx.middleRows(i * xy_height, xy_height);
        Ref<RowMatrixXd> y = Cy.middleRows(i * xy_height, xy_height);

        Ref<RowMatrixXd> z = z_list.middleRows(i * z_height, z_height);

        vector<bool> row_if_nan_x = row_has_nan(x);
        vector<bool> row_if_nan_y = row_has_nan(y);

        int temp = 0;

        // #pragma omp parallel for reduction(+ : temp)
        for (int j = 0; j < z_width; ++j) {
            if (row_if_nan_x[j] || row_if_nan_y[j]) {
                na_records(z_width * i + j) = 1;
                x.row(j) = RowMatrixXd::Zero(1, Cx_width);
                y.row(j) = RowMatrixXd::Zero(1, Cy_width);
                z.col(j) = RowMatrixXd::Zero(z_height, 1);
                if (j < z_information.diff_width)
                    temp += 1;
            }
        }
        num_NA += temp;
        if (temp > max)
            max = temp;
        if (temp < min)
            min = temp;
    }
    int width = z_information.diff_width;
    // if (level)
    // 	width = z_information.level_width;

    // else
    // 	width = z_information.diff_width;
    int nobs = width * N - num_NA;
    int max_obs = width - min;
    int min_obs = width - max;
    double avg_obs = width * 1.0 - (1.0 / N) * num_NA;
    return std::make_tuple(nobs, max_obs, min_obs, avg_obs);
}

std::tuple<RowMatrixXd, RowMatrixXd>
Bootstrapping(string method, int num_draws, int ahead, Ref<RowMatrixXd> z_table, VectorXi &na_records,
              temp_result &tmp_res, struct basic_info &model_info, struct model_options &options) {
    int N = model_info.N;
    vector<int> ids = gen_random_draws(N * num_draws, 0, N - 1);
    vector<RowMatrixXd> All_Mats(num_draws);

    Ref<RowMatrixXd> Cx_list = final_xy_tables["Cx"];
    Ref<RowMatrixXd> Cy_list = final_xy_tables["Cy"];
    // #pragma omp parallel for
    for (int i = 0; i < num_draws; ++i) {
        cout << "bootstrapping " << i << endl;
        RowMatrixXd pseudo_Cx_list, pseudo_Cy_list, pseudo_z_list;
        RowMatrixXd pseudo_xz_list, pseudo_zy_list, pseudo_zHz_list;

        VectorXi pseudo_na_records;

        pseudo_Cx_list = RowMatrixXd::Zero(N * model_info.x_height, model_info.num_indep);
        pseudo_Cy_list = RowMatrixXd::Zero(N * model_info.x_height, model_info.num_dep);
        pseudo_z_list = RowMatrixXd::Zero(N * model_info.num_instr, model_info.x_height);
        pseudo_na_records = VectorXi::Zero(N * model_info.x_height);

        pseudo_xz_list = RowMatrixXd::Zero(tmp_res.xz_list.rows(), tmp_res.xz_list.cols());
        pseudo_zy_list = RowMatrixXd::Zero(tmp_res.zy_list.rows(), tmp_res.zy_list.cols());
        pseudo_zHz_list = RowMatrixXd::Zero(tmp_res.zHz_list.rows(), tmp_res.zHz_list.cols());
        // the height of the three above is z_height *N

        for (int j = 0; j < N; ++j) {

            int the_id = ids[i * N + j];
            ////cout<<"random id: " << the_id <<endl;
            pseudo_Cx_list.middleRows(j * model_info.x_height, model_info.x_height) =
                Cx_list.middleRows(the_id * model_info.x_height, model_info.x_height);
            pseudo_Cy_list.middleRows(j * model_info.x_height, model_info.x_height) =
                Cy_list.middleRows(the_id * model_info.x_height, model_info.x_height);
            pseudo_z_list.middleRows(j * model_info.num_instr, model_info.num_instr) =
                z_table.middleRows(the_id * model_info.num_instr, model_info.num_instr);
            pseudo_na_records.middleRows(j * model_info.x_height, model_info.x_height) =
                na_records.middleRows(the_id * model_info.x_height, model_info.x_height);

            pseudo_xz_list.middleRows(j * model_info.num_instr, model_info.num_instr) =
                tmp_res.xz_list.middleRows(the_id * model_info.num_instr, model_info.num_instr);
            pseudo_zy_list.middleRows(j * model_info.num_instr, model_info.num_instr) =
                tmp_res.zy_list.middleRows(the_id * model_info.num_instr, model_info.num_instr);
            pseudo_zHz_list.middleRows(j * model_info.num_instr, model_info.num_instr) =
                tmp_res.zHz_list.middleRows(the_id * model_info.num_instr, model_info.num_instr);
        }
        temp_result pseudo_t_res = temp_result(tmp_res.H, pseudo_xz_list, pseudo_zy_list, pseudo_zHz_list);
        RowMatrixXd residual, beta;
        std::tie(beta, residual) =
            special_process(pseudo_z_list, pseudo_Cx_list, pseudo_Cy_list, pseudo_t_res, model_info, options);
        All_Mats[i] = irf(method, residual, pseudo_na_records, beta, ahead, model_info.num_dep_lags);
    }

    
    // upper.reserve(model_info.num_dep);
    // lower.reserve(model_info.num_dep);
    // for (int i = 0; i < model_info.num_dep; ++i) {
    //     upper.push_back(RowMatrixXd::Zero(ahead, model_info.num_dep));
    //     lower.push_back(RowMatrixXd::Zero(ahead, model_info.num_dep));
    // }

    
    return choose_U_L(All_Mats, num_draws, model_info.num_dep, ahead);
}

std::tuple<RowMatrixXd, RowMatrixXd>
choose_U_L(vector<RowMatrixXd> &All_Mats,  int num_draws, int num_dep,
           int ahead) {

    int L = num_draws * 0.025 - 1;
    int U = num_draws * 0.975 - 1;

    RowMatrixXd upper(ahead, num_dep*num_dep);
    RowMatrixXd lower(ahead, num_dep*num_dep);

    for(int i=0; i<ahead; ++i){
        for (int j=0; j<num_dep*num_dep;++j){
            vector<double> temp(num_draws);
            for(int m=0; m<num_draws; ++m)
                temp[m]=(All_Mats[m])(i,j);
            std::nth_element(temp.begin(), temp.begin() + L, temp.end());
            
                lower(i, j) = temp[L];

                std::nth_element(temp.begin(), temp.begin() + U, temp.end());
                

                upper(i, j) = temp[U];
            
    }
    }
    

return std::make_tuple(lower, upper);
    
}
