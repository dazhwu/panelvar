#pragma once



#include "pvar.h"
#include "List_Variables.h"
#include "variable.h"
//#include <boost/math/distributions/chi_squared.hpp>



struct model_options {
  int steps;             // = 2;
  bool constant;
  bool level;            // = true;
  bool beginner;         // = false;
  bool timedumm;         // = false;
  bool collapse;         // = false;
  string mmsc;           // = "bic";
  string transformation; // = "fd";
  string irf;
  model_options();
};

extern struct model_options options;

struct df_info {
	int N;
	int T;
	//vector<string> ids;
	int first_diff_index;
	int last_diff_index;
	int first_level_index;
	int last_level_index;
	int max_lag;
	df_info();
	df_info(int _N, int _T, int _f_d_i, int _f_lev_i, int _l_d_i, int _l_lev_i, int _max_lag);
	// last_fod_index;
	//# first_fod_index; int
};

struct z_info {
  int diff_width;
  int diff_height;
  int level_width;
  int level_height;
  int z_width;
  int z_height;
  int num_Dgmm_instr;
  int num_Lgmm_instr;
  int num_instr;
  z_info();
  z_info(int, int, int, int, int, int, int, int);
  //# int num_vars
  //# int num_gmm_instr
};



struct Regression {
   RowMatrixXd beta;
   RowMatrixXd vcov;
   RowMatrixXd std_err;
   RowMatrixXd Z_values;
   RowMatrixXd P_values;
   RowMatrixXd weighting;
   RowMatrixXd Instruments;
   RowMatrixXd Cy;
   RowMatrixXd Cx;
   RowMatrixXd Residual;
   
  Regression();
  Regression(RowMatrixXd &_beta, 
RowMatrixXd &_vcov, 
RowMatrixXd &_std_err,
RowMatrixXd &w,
Ref<RowMatrixXd> _z, Ref<RowMatrixXd> _Cy, Ref<RowMatrixXd> _Cx, RowMatrixXd & _residual);

};

struct MMSC_LU {
  double BIC;
  double HQIC;
  double AIC;
  MMSC_LU();
  MMSC_LU(double bic, double hqic, double aic);
};

struct basic_info {
    string model_name;
    vector<string> identifiers;
    vector<string> dep;
    vector<string> indep;
    int N;
    int T;
    int actual_steps;
    int num_obs;
    int num_instr;
    int num_dep;
    int num_dep_lags;
    int num_indep;
    int diff_width;
	int x_height;
    int max_obs;
    int min_obs;
    double avg_obs;
    RowMatrixXd SS;
    struct MMSC_LU mmsc_lu;
    basic_info();

	basic_info(int _N, int _T, int _num_obs, int _num_instr, int _num_dep, int _num_dep_lags, int _num_indep, int _diff_width, int _x_height, int _max_obs, int _min_obs, double _avg_obs, vector<regular_variable> &dep_indep, vector<string> &identifiers, model_options &options);
};


