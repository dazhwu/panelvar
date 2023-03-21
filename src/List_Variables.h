#pragma once
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include "Common_Functions.h"
#include "pvar.h"

// using std::string;
// using std::vector;

class List_Variables {
  public:
    vector<string> names;
    vector<vector<int>> lags;
    vector<int> min_lags;
    vector<int> max_lags;
    vector<bool> adjustable_min_lags;
    vector<bool> adjustable_max_lags;
    List_Variables();
    // List_Variables(const vector<string> &);
    void append(string, vector<int>, bool, bool, vector<string> &);
    string purge();
};

extern List_Variables Endo_list;
extern List_Variables LaggedDep_list;
extern List_Variables Exog_list;

extern List_Variables IV_list;
extern List_Variables LGMM_list;
extern List_Variables DGMM_list;