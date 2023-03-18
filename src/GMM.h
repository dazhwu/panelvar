//
// Created by Tiger on 5/3/2022.
//

#ifndef UNTITLED_GMM_H
#define UNTITLED_GMM_H
#define EIGEN_INITIALIZE_MATRICES_BY_ZERO

#include <unordered_map>
#include "pvar.h"
#include "info.h"
#include "Step_Result.h"
#include "Common_Functions.h"
#include "Specification_Test.h"
#include <Eigen/KroneckerProduct>

class returned_result {
public:
	struct Regression regression_result;
	vector<vector<RowMatrixXd>> irf;
	Hansen_test_info hansen;
	//vector<AR_test_info> AR_list;
	VectorXcd stability;
	vector<regular_variable> dep_indep;
	struct basic_info model_info;
	struct model_options model_options;
	string command_str;
	returned_result();
	returned_result(struct Regression &, Hansen_test_info &, VectorXcd &, struct basic_info &, struct model_options &, string);

};

class temp_result {
public:
	RowMatrixXd xz_list;
	RowMatrixXd zy_list;
	RowMatrixXd zHz_list;
	RowMatrixXd H;
	temp_result();
	temp_result(RowMatrixXd &, RowMatrixXd &xz_list, RowMatrixXd &zy_list, RowMatrixXd &zHz_list);
};

std::tuple<RowMatrixXd, RowMatrixXd> special_process(Ref<RowMatrixXd> z_list, Ref<RowMatrixXd> Cx_list, Ref<RowMatrixXd> Cy_list,
                                                     temp_result &tmp_res,struct basic_info &model_info, struct model_options &options);

// std::tuple<vector<Step_Result> , Hansen_test_info , vector<AR_test_info> >
std::tuple<returned_result, temp_result> regular_process(Eigen::Ref <RowMatrixXd>,
                                std::unordered_map<string, RowMatrixXd> &, struct basic_info &,
                                struct model_options &);

void core_GMM(temp_result &tmp_res, Ref<RowMatrixXd> Cx_list, Ref<RowMatrixXd> Cy_list, Ref<RowMatrixXd>  z_list, MatrixXd &_z_t_list, int N, int T, int num_obs, int x_width,
              int y_width, int z_height, int diff_width, string transformation, bool level);

Step_Result GMM_step(int N, int num_obs, Ref <RowMatrixXd> z_list,
         MatrixXd &z_list_t, Ref <RowMatrixXd> Cx_list,
         Ref <RowMatrixXd> Cy_list, RowMatrixXd &_XZ, RowMatrixXd &_XZ_t,
         RowMatrixXd &_Zy, RowMatrixXd &_Zy_t, RowMatrixXd &W
);

std::tuple<RowMatrixXd, RowMatrixXd, RowMatrixXd, RowMatrixXd>
calculate_ZuuZ(int N, Ref <RowMatrixXd> z_list, RowMatrixXd &residual, RowMatrixXd &,
               int y_width);

std::tuple<RowMatrixXd, RowMatrixXd, RowMatrixXd, RowMatrixXd, RowMatrixXd>
calculate_XZ_ZY_W1(RowMatrixXd &xz_list, RowMatrixXd &zy_list, RowMatrixXd &zHz_list,
                   int N, int x_width, int y_width, int z_height);

std::tuple<RowMatrixXd, RowMatrixXd, RowMatrixXd> calculate_basic(Ref <RowMatrixXd> z_list,
                                                                  Ref <RowMatrixXd> Cx_list,
                                                                  Ref <RowMatrixXd> Cy_list, RowMatrixXd &H,
                                                                  int N, int z_height, int x_height, int x_width, int y_width);

//RowMatrixXd calculate_W(Ref <RowMatrixXd> H,
//                        Ref <RowMatrixXd> z_list,
//                        Ref <MatrixXd> _z_t_list, int step,
//                        int N);

RowMatrixXd calculate_residual(Ref <RowMatrixXd> y_list,
                               Ref <RowMatrixXd> x_list,
                               Ref <RowMatrixXd> beta);

RowMatrixXd vcov_step_1(Ref <RowMatrixXd> _M_XZ_W, Ref <RowMatrixXd> W2, int N);

RowMatrixXd vcov(Ref <RowMatrixXd> z_list,
                 Ref <RowMatrixXd> Cx, Ref <RowMatrixXd> M2, Ref <RowMatrixXd> _M2_XZ_W2, Ref <RowMatrixXd> _W2_inv,
                 Ref <RowMatrixXd> zs2,
                 int step, int N);

void update_MMSC_LU(struct basic_info &model_info, Hansen_test_info &hansen);

Hansen_test_info perform_hansen_test(int step, int num_instru, int num_indep, int num_dep, int N);

// vector <AR_test_info> perform_AR_test(std::unordered_map<string, RowMatrixXd> &,
//                                       Ref <RowMatrixXd> z_list, Ref <RowMatrixXd> _zs_list,
//                                       int step, int diff_width, int N, string transformation, bool level);

/*
std::tuple<Hansen_test_info, vector<AR_test_info>>  perform_test( Ref<RowMatrixXd> ,  Ref<RowMatrixXd> ,
                     Ref<RowMatrixXd> , Ref<RowMatrixXd> ,
                  int num_instru, int num_indep, int step, int didd_width,	int N, string transformation, bool level);
*/
RowMatrixXd get_H1(int width, int diff_width, int T,
                   string transformation, bool level);

RowMatrixXd get_H1_fod(int width, int diff_width, int T, bool level);

RowMatrixXd generate_D_matrix(int height, int T);

RowMatrixXd Windmeijer(Ref <RowMatrixXd> M2, Ref <RowMatrixXd> _M2_XZ_W2,
                       Ref <RowMatrixXd> W2_inv, Ref <RowMatrixXd> zs2,
                       Ref <RowMatrixXd> vcov_step1, Ref <RowMatrixXd> Cx_list,
                       Ref <RowMatrixXd> z_list, Ref <RowMatrixXd> residual1_t, int N);

#endif //UNTITLED_GMM_H