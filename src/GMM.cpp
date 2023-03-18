
#include "GMM.h"

#define EIGEN_INITIALIZE_MATRICES_BY_ZERO
std::vector<Step_Result> results;
int steps;
int num_dep_lags;
// std::tuple<vector<Step_Result> , Hansen_test_info , vector<AR_test_info>

returned_result::returned_result() {

}
returned_result::returned_result(struct Regression &reg,
                                 Hansen_test_info &_hansen,
                                 VectorXcd &_stability,
                                 struct basic_info &mo_info, struct model_options &_m_options,  string com_str) {
	regression_result = reg;
	hansen = _hansen;
	stability = _stability;
	model_info = mo_info;
	model_options = _m_options;
	//irf=_irf;
	command_str = com_str;
}

temp_result::temp_result() {

}

temp_result::temp_result(RowMatrixXd &_H, RowMatrixXd &_xz_list, RowMatrixXd &_zy_list, RowMatrixXd &_zHz_list) {
	H = _H;
	xz_list = _xz_list;
	zy_list = _zy_list;
	zHz_list = _zHz_list;
}

std::tuple<RowMatrixXd, RowMatrixXd> special_process(Ref<RowMatrixXd> z_list, Ref<RowMatrixXd> Cx_list, Ref<RowMatrixXd> Cy_list,
								temp_result &tmp_res,struct basic_info &model_info, struct model_options &options){
	string transformation = options.transformation;
	steps = options.steps;
	bool level = options.level;
	int N = model_info.N;
	int T = model_info.T;
	int num_steps = options.steps;
	num_dep_lags = model_info.num_dep_lags;
	// int num_instru = model_info.num_instr;
	// int num_indep = model_info.num_indep;
	//int num_obs = model_info.num_obs;
	//int diff_width = model_info.diff_width;

//	int z_height = model_info.num_instr;
//	int x_height = model_info.x_height;
//	int x_width = model_info.num_indep;
//	int y_width = model_info.num_dep;


	MatrixXd _z_t_list = z_list.transpose();
	results.clear();


	core_GMM(tmp_res, Cx_list, Cy_list, z_list, _z_t_list, N, T, model_info.num_obs, model_info.num_indep,
	         model_info.num_dep, model_info.num_instr, model_info.diff_width, transformation, level);

	return std::make_tuple(results[num_steps-1].beta, results[num_steps-1].residual);

}

std::tuple<returned_result, temp_result>
regular_process(Eigen::Ref <RowMatrixXd> z_list,
                std::unordered_map<string, RowMatrixXd> &final_xy_tables,
                struct basic_info &model_info, struct model_options &options) {
	Eigen::Ref <RowMatrixXd> Cx_list = final_xy_tables["Cx"];
	Eigen::Ref <RowMatrixXd> Cy_list = final_xy_tables["Cy"];

	string transformation = options.transformation;
	steps = options.steps;
	bool level = options.level;
	int N = model_info.N;
	int T = model_info.T;
	int num_steps = options.steps;
	num_dep_lags = model_info.num_dep_lags;
	// int num_instru = model_info.num_instr;
	// int num_indep = model_info.num_indep;
	//int num_obs = model_info.num_obs;
	int diff_width = model_info.diff_width;

	int z_height = model_info.num_instr;
	int x_height = model_info.x_height;
	int x_width = model_info.num_indep;
	int y_width = model_info.num_dep;


	// //saveData("cy.csv", Cy_list);

	// //saveData("z.csv", z_list);




	MatrixXd _z_t_list = z_list.transpose();
	results.clear();
	////cout <<"----calculate basic" << endl;
	RowMatrixXd xz_list, zy_list, H1, zHz_list;
	H1 = get_H1(z_list.cols(), diff_width, T, transformation, level);
	std::tie(xz_list, zy_list, zHz_list) = calculate_basic(z_list, Cx_list, Cy_list, H1, N, z_height, x_height, x_width, y_width);
	temp_result tmp_res=temp_result(H1, xz_list, zy_list, zHz_list);
	core_GMM(tmp_res, Cx_list, Cy_list, z_list, _z_t_list, N, T, model_info.num_obs, model_info.num_indep,
	         model_info.num_dep, model_info.num_instr, model_info.diff_width, transformation, level);



	RowMatrixXd _vcov;
	_vcov = vcov_step_1(results[0]._M_XZ_W, results[0].W_next, N);
	results[0].add_vcov(_vcov);

	if (num_steps == 2) {
		_vcov = vcov(z_list, Cx_list, results[1].M, results[1]._M_XZ_W, results[1].W_inv, results[1].qs, 2, N);
		results[1].add_vcov(_vcov);

	}
	Hansen_test_info hansen =
			hansen_overid(results[1].W_inv, results[1].zs, model_info.num_instr * model_info.num_dep, model_info.num_indep * model_info.num_dep, N);
// 考虑把W_inv计算加到W_next后边去   //！！！！查在onestep条件下，是否用result2.zs, 关系到是否把 ZuuZ计算单列


	//查在onestep条件下，是否用result2.beta
	VectorXcd stability =
			stability_test(results[1].beta, model_info.num_dep_lags);

	Regression tbr =
			Regression(results[num_steps - 1].beta, results[num_steps - 1].vcov,
			           results[num_steps - 1].std_err,
			           results[num_steps - 1].W, z_list, Cy_list, Cx_list, results[num_steps - 1].residual
			);  //_residual_t

	//model_info.SS = results[num_steps - 1].SS;
	//model_info.actual_steps = steps;
	model_info.actual_steps = num_steps;
	update_MMSC_LU(model_info, hansen);
	////cout << "MMSC_LU " <<endl;
	returned_result the_result =
			returned_result(tbr, hansen, stability, model_info, options, "");
	return std::make_tuple(the_result, tmp_res);

}

void core_GMM(temp_result &tmp_res, Ref<RowMatrixXd> Cx_list, Ref<RowMatrixXd> Cy_list, Ref<RowMatrixXd>  z_list, MatrixXd &_z_t_list, int N, int T, int num_obs, int x_width,
              int y_width, int z_height, int diff_width, string transformation, bool level) {

	RowMatrixXd _XZ, _Zy, _XZ_t, _Zy_t, W1;
	std::tie(_XZ, _Zy, _XZ_t, _Zy_t, W1) = calculate_XZ_ZY_W1(tmp_res.xz_list, tmp_res.zy_list, tmp_res.zHz_list, N, x_width, y_width, z_height);



////cout <<"----GMM" << endl;
	Step_Result result1 = GMM_step(N,  num_obs, z_list, _z_t_list, Cx_list, Cy_list, _XZ,
	                               _XZ_t, _Zy, _Zy_t, W1);

	Step_Result result2 = GMM_step(N,  num_obs,  z_list, _z_t_list, Cx_list, Cy_list, _XZ,
	                               _XZ_t, _Zy, _Zy_t, result1.W_next);

	results.push_back(result1);
	results.push_back(result2);

}

void update_MMSC_LU(struct basic_info &model_info, Hansen_test_info &hansen) {
	double log_n = log(model_info.num_obs);
	// double dif = (model_info.num_instr -
	// model_info.num_indep)*model_info.num_dep;
	double dif = model_info.num_instr * model_info.num_dep - model_info.num_indep;
	model_info.mmsc_lu.BIC = hansen.test_value - (dif) * log_n;
	model_info.mmsc_lu.HQIC = hansen.test_value - dif * log(log_n) * 2.1;
	model_info.mmsc_lu.AIC = hansen.test_value - (dif) * 2.0;

}

Step_Result GMM_step(int N, int num_obs, Ref <RowMatrixXd> z_list,
                     MatrixXd &z_list_t, Ref <RowMatrixXd> Cx_list,
                     Ref <RowMatrixXd> Cy_list, RowMatrixXd &_XZ, RowMatrixXd &_XZ_t,
                     RowMatrixXd &_Zy, RowMatrixXd &_Zy_t, RowMatrixXd &W
) {
	RowMatrixXd W_inv = common_inv(W);
	RowMatrixXd _XZ_W = _XZ * W_inv;

	RowMatrixXd _M_inv = _XZ_W * _XZ_t;
	RowMatrixXd M = common_inv(_M_inv);

	RowMatrixXd _M_XZ_W = M * _XZ_W;

	RowMatrixXd beta =
			(_M_XZ_W * _Zy_t)
					.reshaped<Eigen::RowMajor>(Cx_list.cols(), Cy_list.cols());

////cout << beta << endl;
	RowMatrixXd residual = calculate_residual(Cy_list, Cx_list, beta);

	RowMatrixXd _residual_t = residual.transpose();
	// RowMatrixXd SS = (_residual_t * residual) * (1.0 / 2 / num_obs);
	RowMatrixXd diag_dep = RowMatrixXd::Identity(Cy_list.cols(), Cy_list.cols());
	RowMatrixXd qs, _zs_list, zs, ZuuZ;
	std::tie(qs, _zs_list, zs, ZuuZ) = calculate_ZuuZ(N, z_list, residual, diag_dep, Cy_list.cols());
//  //cout << ZuuZ << endl;
	RowMatrixXd W_next = ZuuZ * (1.0 / N);
/*
  RowMatrixXd _vcov;

  if (step == 1)
    _vcov = vcov_step_1(_M_XZ_W, W_next, N);
  else
    _vcov = vcov(z_list, Cx_list, M, _M_XZ_W, W_inv, qs, step, N);

  ////cout << _vcov << endl;
  RowMatrixXd std_err(Cx_list.cols(), Cy_list.cols()), temp;

  temp = _vcov.diagonal();

  for (int i = 0; i < Cx_list.cols(); ++i)
    for (int j=0; j< Cy_list.cols(); ++j)
      std_err(i, j) = sqrt(temp(i*Cy_list.cols()+j, 0));
  //cout<<std_err<<endl;
*/
	Step_Result current_step =
			Step_Result(W, W_inv, W_next, _XZ_W, M, _M_XZ_W, zs, qs, ZuuZ,
			            _zs_list, beta, residual, _residual_t);
	return (current_step);
}

std::tuple<RowMatrixXd, RowMatrixXd, RowMatrixXd, RowMatrixXd>
calculate_ZuuZ(int N, Ref <RowMatrixXd> z_list, RowMatrixXd &residual, RowMatrixXd &diag_dep,
               int y_width) {
	int z_height = z_list.rows() / N;
	// int z_width = z_list.cols();
	int r_height = residual.rows() / N;

	int zs_height = z_height * y_width;
	RowMatrixXd _zs_list = RowMatrixXd::Zero(N * zs_height, 1);
	RowMatrixXd zs = RowMatrixXd::Zero(zs_height, 1);
	RowMatrixXd qs = RowMatrixXd::Zero(z_height, y_width);
	RowMatrixXd ZuuZ = RowMatrixXd::Zero(zs_height, zs_height);
#pragma omp parallel for reduction(+ : zs, qs, ZuuZ)
	for (int i = 0; i < N; ++i) {
		Map <RowMatrixXd> u(residual.middleRows(i * r_height, r_height).data(),
		                    r_height * y_width, 1);

//    if(i==0)
//      //cout <<z_list.middleRows(i * z_height, z_height) <<endl;

		qs.noalias() += z_list.middleRows(i * z_height, z_height) *
		                residual.middleRows(i * r_height, r_height);
		_zs_list.middleRows(i * zs_height, zs_height) =
				Eigen::kroneckerProduct(z_list.middleRows(i * z_height, z_height), diag_dep) * u;
		Ref <RowMatrixXd> temp_zs = _zs_list.middleRows(i * zs_height, zs_height);
		zs.noalias() += temp_zs;
		ZuuZ.noalias() += temp_zs * temp_zs.transpose();
	}
	return std::make_tuple(qs, _zs_list, zs, ZuuZ);
}

//
std::tuple<RowMatrixXd, RowMatrixXd, RowMatrixXd, RowMatrixXd, RowMatrixXd>
calculate_XZ_ZY_W1(RowMatrixXd &xz_list, RowMatrixXd &zy_list, RowMatrixXd &zHz_list,
                   int N, int x_width, int y_width, int z_height) {
	RowMatrixXd temp_xz = RowMatrixXd::Zero(z_height, x_width);
	RowMatrixXd temp_zy = RowMatrixXd::Zero(z_height, y_width);
	RowMatrixXd temp_W = RowMatrixXd::Zero(z_height, z_height);
	for (int i = 0; i < N; ++i) {
		Ref <RowMatrixXd> zx_block = xz_list.middleRows(i * z_height, z_height);
		Ref <RowMatrixXd> zy_block = zy_list.middleRows(i * z_height, z_height);
		Ref <RowMatrixXd> zHz_block = zHz_list.middleRows(i * z_height, z_height);
		temp_xz += zx_block;   //  .transpose();
		temp_zy += zy_block; //.transpose();
		temp_W += zHz_block;
	}
	RowMatrixXd diag_dep = RowMatrixXd::Identity(y_width, y_width);
	RowMatrixXd tbr_xz = Eigen::kroneckerProduct(temp_xz, diag_dep);

	RowMatrixXd W1 = Eigen::kroneckerProduct(temp_W, diag_dep);
	Map <RowMatrixXd> tbr_zy(temp_zy.data(), temp_zy.cols() * temp_zy.rows(), 1);
	RowMatrixXd t_xz = tbr_xz.transpose();
	RowMatrixXd t_zy = tbr_zy.transpose();
	return std::make_tuple(t_xz, t_zy, tbr_xz, tbr_zy, W1);

}

std::tuple<RowMatrixXd, RowMatrixXd, RowMatrixXd> calculate_basic(Ref <RowMatrixXd> z_list,
                                                                  Ref <RowMatrixXd> Cx_list,
                                                                  Ref <RowMatrixXd> Cy_list, RowMatrixXd &H,
                                                                  int N, int z_height, int x_height, int x_width, int y_width) {
	RowMatrixXd zx_list = RowMatrixXd::Zero(N * z_height, x_width);
	RowMatrixXd zy_list = RowMatrixXd::Zero(N * z_height, y_width);
	RowMatrixXd zHz_list = RowMatrixXd::Zero(N * z_height, z_height);

	//#pragma omp parallel for reduction(+:temp_xz, temp_zy)
	for (int i = 0; i < N; ++i) {
		Ref <RowMatrixXd> z = z_list.middleRows(z_height * i, z_height);
		Ref <RowMatrixXd> x = Cx_list.block(x_height * i, 0, x_height, x_width);
		Ref <RowMatrixXd> y = Cy_list.block(x_height * i, 0, x_height, y_width);

		Ref <RowMatrixXd> zx = zx_list.middleRows(i * z_height, z_height);
		Ref <RowMatrixXd> zy = zy_list.middleRows(i * z_height, z_height);
		Ref <RowMatrixXd> zHz = zHz_list.middleRows(i * z_height, z_height);
		zx = z * x;
		zy = z * y;
		zHz = z * H * z.transpose();
	}
	return std::make_tuple(zx_list, zy_list, zHz_list);

}

RowMatrixXd calculate_residual(Ref <RowMatrixXd> y_list, Ref <RowMatrixXd> x_list,
                               Ref <RowMatrixXd> beta) {
	RowMatrixXd tbr = y_list - x_list * beta;
	//saveData("residual.csv", tbr);
	return tbr;
}

RowMatrixXd vcov_step_1(Ref <RowMatrixXd> _M_XZ_W, Ref <RowMatrixXd> W2, int N) {
	// Ref<RowMatrixXd> W2 = step_1.W_next;
	RowMatrixXd tbr = N * (_M_XZ_W * W2 * _M_XZ_W.transpose());
	return tbr;
}

RowMatrixXd vcov(Ref <RowMatrixXd> z_list, Ref <RowMatrixXd> Cx,
                 Ref <RowMatrixXd> M2, Ref <RowMatrixXd> _M2_XZ_W2,
                 Ref <RowMatrixXd> _W2_inv, Ref <RowMatrixXd> zs2, int step,
                 int N) {
	RowMatrixXd tbr;

	Step_Result previous_step = results[step - 2];
	Ref <RowMatrixXd> vcov_step_previous = previous_step.vcov;
	Ref <RowMatrixXd> residual1 = previous_step.residual;
	tbr = Windmeijer(M2, _M2_XZ_W2, _W2_inv, zs2, vcov_step_previous, Cx, z_list,
	                 residual1, N);
	return tbr;
}

/*
Hansen_test_info perform_hansen_test(int step, int num_instru, int num_indep,
                                     int num_dep, int N) {
	Step_Result step1 = results[0];
	Step_Result step2 = results[1];
	Step_Result current_step = step2; // to be changed
	//	if (step <= 2)
	//		current_step = step2;
	//	else
	//		current_step = results[step - 1];

	Ref <RowMatrixXd> W2_inv = current_step.W_inv;
	Ref <RowMatrixXd> zs = current_step.zs;

	Hansen_test_info hansen =
			hansen_overid(W2_inv, zs, num_instru * num_dep, num_indep * num_dep, N);
	return hansen;
}
*/
RowMatrixXd get_H1(int width, int diff_width, int T,
                   string transformation, bool level) {
	//int width = z_list.cols();
	RowMatrixXd tbr;
	if (transformation == "fd") {
		tbr = RowMatrixXd::Zero(width, width);
		for (int i = 0; i < diff_width; ++i) {
			tbr(i, i) = 2;
			if (i >= 1)
				tbr(i - 1, i) = -1;
			if (i < diff_width - 1)
				tbr(i + 1, i) = -1;
		}
		if (width > diff_width) {
			for (int i = diff_width; i < width; ++i)
				tbr(i, i) = 1;
			Eigen::Ref <RowMatrixXd> low_left =
					tbr.block(diff_width, 0, width - diff_width, diff_width);
			for (int i = 0; i < diff_width; ++i) {
				low_left(i, i) = -1;
				low_left(i + 1, i) = 1;
			}
			Eigen::Ref <RowMatrixXd> up_right =
					tbr.block(0, diff_width, diff_width, width - diff_width);
			for (int i = 0; i < diff_width; ++i) {
				up_right(i, i) = -1;
				up_right(i, i + 1) = 1;
			}
		}
	} else {
		tbr = get_H1_fod(width, diff_width, T, level);
	}
	return tbr;
}

RowMatrixXd get_H1_fod(int width, int diff_width, int T, bool level) {
	RowMatrixXd D_up = generate_D_matrix(diff_width, T);
	if (level) {
		RowMatrixXd D = RowMatrixXd::Zero(width, T);
		D.block(0, 0, diff_width, T) = D_up;
		Eigen::Ref <RowMatrixXd> low_right =
				D.block(diff_width, T - (width - diff_width), width - diff_width,
				        width - diff_width);
		for (int i = 0; i < width - diff_width; ++i) {
			low_right(i, i) = 1.0;
		}
		////cout<<D<<endl;
		return D * D.transpose();
	} else {
		return D_up * D_up.transpose();
	}
}

RowMatrixXd generate_D_matrix(int height, int T) {
	//	//cout<<height<<endl;
	RowMatrixXd temp = RowMatrixXd::Zero(T, T);
	RowMatrixXd D = RowMatrixXd::Zero(height, T);
	int start;
	for (int i = 0; i < T; ++i) {
		for (int j = i; j < T; ++j) {
			if (i == j)
				temp(i, j) = sqrt((T - i - 1) * 1.0 / (T - i));
			else
				temp(i, j) = (-1) * sqrt(1.0 / ((T - i) * 1.0 * (T - i - 1)));
		}
	}
	start = num_dep_lags - 1;
	////cout<<temp<<endl;
	D = temp.block(start, 0, height, T);
	return D;
}

RowMatrixXd Windmeijer(Ref <RowMatrixXd> M2, Ref <RowMatrixXd> _M2_XZ_W2,
                       Ref <RowMatrixXd> W2_inv, Ref <RowMatrixXd> zs2,
                       Ref <RowMatrixXd> vcov_step1, Ref <RowMatrixXd> Cx_list,
                       Ref <RowMatrixXd> z_list, Ref <RowMatrixXd> residual1,
                       int N) {
	int x_height = Cx_list.rows() / N;
	int z_height = z_list.rows() / N;
	int x_width = Cx_list.cols();
	// int zs_height= zs2.rows();
	int y_width = residual1.cols();

	RowMatrixXd diag_dep = RowMatrixXd::Identity(y_width, y_width);
	RowMatrixXd D = RowMatrixXd::Zero(z_height * y_width * y_width * x_width,
	                                  z_height * y_width);
	int z_width = z_list.cols();

	RowMatrixXd zx;
	for (int i = 0; i < N; ++i) {
		Ref <RowMatrixXd> x = Cx_list.middleRows(i * x_height, x_height);
		Ref <RowMatrixXd> z = z_list.block(i * z_height, 0, z_height, z_width);

		Ref <RowMatrixXd> u = residual1.middleRows(i * x_height, x_height);

		RowMatrixXd xu; // = x * u_t;

		RowMatrixXd temp_xz = z * x;
		zx = Eigen::kroneckerProduct(
				temp_xz, diag_dep); // Eigen::kroneckerProduct(z*x, diag_dep);



		RowMatrixXd temp_uz = z * u;

		Map <RowMatrixXd> u_z(temp_uz.data(), 1, zx.rows());
		for (int j = 0; j < zx.cols(); ++j) {
			Ref <RowMatrixXd> D_j = D.middleRows(j * zx.rows(), zx.rows());
			D_j.noalias() += zx.col(j) * u_z + (zx.col(j) * u_z).transpose();
		}
	}
	D = (-1.0 / N) * D;
	RowMatrixXd D_W = RowMatrixXd::Zero(M2.rows(), M2.cols());
	Map <RowMatrixXd> zs_vec(zs2.data(), zs2.rows() * zs2.cols(), 1);
	for (int k = 0; k < D_W.cols(); ++k) {
		Ref <RowMatrixXd> partial_dir = D.middleRows(k * zx.rows(), zx.rows());

		RowMatrixXd tt = _M2_XZ_W2 * partial_dir * W2_inv * zs_vec;
		D_W.col(k) = (-1) * _M2_XZ_W2 * partial_dir * W2_inv * zs_vec;

	}
	RowMatrixXd temp = N * M2 + D_W * N * M2 + N * M2 * D_W.transpose();
	temp += D_W * vcov_step1 * D_W.transpose();
	return temp;
}
