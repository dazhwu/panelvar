#pragma once
#include <unordered_map>
#include <math.h>
#include <random>
#include <iostream>
#include <bits/stdc++.h>

#include <algorithm>

#include "pvar.h"
#include "info.h"
#include "model_organizer.h"
#include "variable.h"
#include "instruments.h"
#include "Common_Functions.h"
#include "List_Variables.h"
#include "Command.h"
#include "GMM.h"

using Eigen::VectorXi;


struct df_info get_info(model_options options, model m);

vector<returned_result> prepare_data(vector<string> &,Ref <RowMatrixXd>,    vector<int> &,  vector<string> &, vector<string> &, int, int);

void update_time_dummies(vector<string> &col_timedumm, struct df_info &info,
                         model &m);
void process_residual(returned_result &, VectorXi &);

void handle_tables(Ref<RowMatrixXd> pdata, Ref<RowMatrixXd> fd_data, model &m, vector<string> &df_cols, df_info info);

struct basic_info generate_model_info(int N, int T, Ref<RowMatrixXd> z_table, VectorXi &na_records, struct z_info &z_information, model &m, vector<string> &);

void get_gmm_tables(Ref<RowMatrixXd> , Ref<RowMatrixXd> ,vector <gmm_var> &, vector <regular_variable> &,
                    vector <gmm_var> &, vector<string> &,bool level) ;
template<class ty>
RowMatrixXd gen_table(Ref <RowMatrixXd> ori_data, vector <ty> variable_list,
                      vector <string> df_cols) ;

void get_xy_table_dict(Ref <RowMatrixXd> , model &, vector<string> &, string );

void get_final_xy_tables(df_info info, bool level, string transformation) ;

std::tuple <RowMatrixXd, RowMatrixXd> get_final_xy_systemGMM(
		Ref <RowMatrixXd> Dx, Ref <RowMatrixXd> Dy, Ref <RowMatrixXd> x, Ref <RowMatrixXd> y,
		vector<int> cut, vector<int> Dcut, int cut_height, int Dcut_height, string transformation) ;



std::tuple <RowMatrixXd, RowMatrixXd>
get_final_xy_diffGMM(Ref <RowMatrixXd> Dy, Ref <RowMatrixXd> Dx, vector<int> Dcut, int Dcut_height);

void prepare_reg_fod(Ref <RowMatrixXd> Diff_x, Ref <RowMatrixXd> Diff_y);

std::tuple<int, int, int, double>
prepare_reg(Ref <RowMatrixXd> z_list, Ref <RowMatrixXd> Cx, Ref <RowMatrixXd> Cy, VectorXi &,z_info,
             string transformation, bool level);

std::tuple<vector<RowMatrixXd>,vector<RowMatrixXd> >  Bootstrapping(string, int , int ,  Ref<RowMatrixXd> z_table,
                   VectorXi &na_records, temp_result &tmp_res,struct basic_info &model_info, struct model_options &options);
void choose_U_L(vector<vector<RowMatrixXd>> &All_Mats, vector<RowMatrixXd> &upper, vector<RowMatrixXd> &lower, int num_draws, int num_dep, int ahead );
vector<int> gen_random_draws(int total_num_draws, int from, int to);