#pragma once
#include "List_Variables.h"
#include "pvar.h"
#include "variable.h"

class model {
  public:
    string name;
    int num_dep;
    int num_dep_lags;
    vector<regular_variable> dep_indep, iv_vars;
    vector<gmm_var> Dgmm_vars, Lgmm_vars;
    string command_str;
    // model();
};

// extern vector<model> list_models;
vector<model> generate_list_models(int);
// void list_to_gmm_m(int , vector<model> &, List_Variables &);

void update_commandStr(model &m);
void list_to_Dep(int T, vector<model> &list_models);
void add_Dep(model &m, int first_lag, int last_lag);

void list_to_exog_iv(int T, vector<model> &list_models, List_Variables &the_list, int dest);
void list_to_gmm(int T, vector<model> &list_models, List_Variables &the_list);

void add_gmm(model &m, string var_name, int min_lag, int max_lag);

void explode_gmm(int min_lag, int T, vector<model> &list_models, string var_name);

void explode_exog_iv(int max_lag, int T, vector<model> &list_models, string var_name, int dest);
