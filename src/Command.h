#pragma once

#include "Common_Functions.h"
#include "List_Variables.h"
#include "pvar.h"
#include "info.h"
#include <algorithm>
#include <climits>
#include <numeric>
#include <regex>
#include <string>
#include <vector>

using std::regex;
using std::smatch;
using std::string;
using std::vector;

class Command {
public:
  int largest_T;

  // string command_str;
  vector<string> cols;

  Command(int, string, int, string, string, string, vector<string> &);

  void parse_dep(string, int);

  void parse_exog(string);

  void parse_gmm_iv(string);

  // bool variable_exists(string);
  bool parse_spaced_vars_range(string, List_Variables &);

  bool parse_spaced_vars_single(string, List_Variables &);

  bool parse_spaced_vars_auto(string, List_Variables &);

  void parse_spaced_vars(std::vector<std::string> list_vars,
                         List_Variables &dest_list);

  void parse_gmmStyle(vector<string> &, string);

  void process_GMM(vector<string> &, int, int, string, bool);

  void parse_IV(vector<string> &, string);

  void parse_options(string);

  void check_dep();

  void check_exog();

  void check_GMM();

  void check_iv();

  void check_three_lists();
};

std::tuple<vector<string>, struct model_options>
process_command(int, string, int, string, string, string, vector<string> &);
