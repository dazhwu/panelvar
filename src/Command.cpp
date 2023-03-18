#include "Command.h"

using namespace std;
List_Variables Endo_list;
List_Variables LaggedDep_list;
List_Variables Exog_list;
List_Variables IV_list;
List_Variables LGMM_list;
List_Variables DGMM_list;

struct model_options options;
int DEP_lags;

std::tuple<vector<string>, struct model_options>
process_command(int T_MAX, string dep, int lags, string exog, string gmm,
                string options_str, vector<string> &df_col_names) {
  IV_list = List_Variables();
  Endo_list = List_Variables();
  LaggedDep_list=List_Variables();
  Exog_list = List_Variables();
  LGMM_list = List_Variables();
  DGMM_list = List_Variables();
  options = model_options();

  DEP_lags=lags;

  Command new_command = Command(T_MAX, dep, lags, exog, gmm, options_str, df_col_names);

  vector<string> names = combined_vector(Endo_list.names, Exog_list.names);

  names = combined_vector(names, IV_list.names);
  names = combined_vector(names, DGMM_list.names);
  names = combined_vector(names, LGMM_list.names);

  return std::make_tuple(names, options);
}
// pvar_module.process_command(pdata.T, dep, lags, exog, gmm, options_str,
// df.columns)
Command::Command(int T_MAX, string dep, int lags, string exog, string gmm,
                 string options_str, vector<string> &df_col_names) {
  // command_str = commandstr;
  largest_T = T_MAX;
  cols = df_col_names;

  parse_dep(dep, lags);
  parse_exog(exog);
  parse_gmm_iv(gmm);
  parse_options(options_str);
  check_dep();
  check_exog();
  check_GMM();
  check_iv();
  check_three_lists();
}

void Command::parse_dep(string dep, int lags) {

  vector<string> list_vars = splitString(dep, ' ');
  parse_spaced_vars(list_vars, Endo_list);

  parse_spaced_vars(list_vars, LaggedDep_list);

  for (std::size_t i = 0; i != LaggedDep_list.names.size(); ++i){
    vector<int> Livec(lags);
    std::iota(Livec.begin(), Livec.end(), 1);
    LaggedDep_list.lags[i]=Livec;
  }




}

void Command::parse_exog(string exog) {
  vector<string> list_vars = splitString(exog, ' ');
  parse_spaced_vars(list_vars, Exog_list);
}

void Command::parse_gmm_iv(string gmm) {
  vector<string> matching_parts;
  parse_gmmStyle(matching_parts, gmm);
  //parse_endo_pred(matching_parts);
  parse_IV(matching_parts, gmm);
}

bool Command::parse_spaced_vars_range(string var, List_Variables &dest_list) {
  string lag_range =
      "^L[(]([0-9]{1,})[:]([0-9]{1,})[)][.]([a-zA-Z_]{1,}[a-zA-Z_0-9]{0,})$";
  regex pt_range(lag_range);
  smatch match_groups_multiple;
  bool isMatch = regex_match(var, match_groups_multiple, pt_range);
  if (isMatch) {
    int LB = std::stoi(match_groups_multiple[1]);
    int UB = std::stoi(match_groups_multiple[2]);
    if (LB > UB) {
      int temp = LB;
      LB = UB;
      UB = temp;
    }
    string name = match_groups_multiple[3];

    vector<int> ivec(UB - LB + 1);
    std::iota(ivec.begin(), ivec.end(), LB);
    dest_list.append(name, ivec, false, false, cols);
    return true;
  } else {
    return false;
  }
}

bool Command::parse_spaced_vars_single(string var, List_Variables &dest_list) {
  string lag_single = "^L([0-9]{1,})[.]([a-zA-Z_]{1,}[a-zA-Z_0-9]{0,})$";
  regex pt_single(lag_single);
  smatch match_groups_single;
  bool isMatch = regex_match(var, match_groups_single, pt_single);
  if (isMatch) {
    int lag = stoi(match_groups_single[1]);
    string name = match_groups_single[2];

    vector<int> ivec(1);
    ivec[0] = lag;
    dest_list.append(name, ivec, false, false, cols);
    return true;
  } else
    return false;
}

bool Command::parse_spaced_vars_auto(string var, List_Variables &dest_list) {
  string lag_auto =
      "^L[(]([0-9]{1,})[:]([?])[)][.]([a-zA-Z_]{1,}[a-zA-Z_0-9]{0,})$";
  regex pt_auto(lag_auto);
  smatch match_groups_auto;
  bool isMatch = regex_match(var, match_groups_auto, pt_auto);
  if (isMatch) {
    int LB = stoi(match_groups_auto[1]);
    string name = match_groups_auto[3];
    options.beginner = true;
    vector<int> ivec(1);
    ivec[0] = LB;  
    dest_list.append(name, ivec, false, true, cols);
    return true;
  } else
    return false;
}

// bool Command::variable_exists(string name) {
//     if (std::find(cols.begin(), cols.end(), name) == cols.end()) {
//         return false;
//     }
//     return true;
// }
void Command::parse_spaced_vars(vector<string> list_vars,
                                List_Variables &dest_list) {
  bool isMatch;
  for (string var : list_vars) {
    isMatch = parse_spaced_vars_range(var, dest_list); // e.g., L(1:3).variable
    if (!isMatch) {
      isMatch = parse_spaced_vars_single(var, dest_list); // e.g., L1.variable
      if (!isMatch)
        isMatch =
            parse_spaced_vars_auto(var, dest_list); // e.g., L(1:?).variable
      if (!isMatch) { // if variable name is not in the forms above, then just
                      // assume no lag
        vector<int> ivec(1);
        ivec[0] = 0;
        dest_list.append(var, ivec, false, false, cols);
      }
    }
  }
}

void Command::parse_gmmStyle(vector<string> &matching_parts, string gmm) {
  smatch match;
  string temp = gmm;
  regex r("gmm[(][a-zA-Z_0-9 ]{1,}[,][ ]{0,}[0-9]{1,}[ ]{0,}[:][ "
          "]{0,}(?:(?:[.])|(?:[?])|(?:[0-9]{1,}))[ ]{0,}[)]");
  int min_lag;//, max_lag;
  // int i = 1;
  while (regex_search(temp, match, r)) {
    string part = match.str(0);
    matching_parts.push_back(part);
    smatch match_groups_multiple;
    regex prog_1("^gmm[(]([a-zA-Z_0-9 ]{1,})[,][ ]{0,}([0-9]{1,})[ ]{0,}[:][ "
                 "]{0,}((?:[.])|(?:[?])|(?:[0-9]{1,}))[ ]{0,}[)]$");
    regex_match(part, match_groups_multiple, prog_1);
    vector<string> vars = splitString(match_groups_multiple[1], ' ');
    min_lag = stoi(match_groups_multiple[2]);
    if (match_groups_multiple[3] == '.')
      process_GMM(vars, min_lag, largest_T, part, false);      
    else if(match_groups_multiple[3] == '?')
      process_GMM(vars, min_lag, min_lag, part, true);      
    else
      process_GMM(vars, min_lag, stoi(match_groups_multiple[3]), part, false);            
    
    temp = match.suffix().str();
  }
}

void Command::process_GMM(vector<string> &vars, int min_lag, int max_lag,
                          string part, bool auto_mode) {
  if (min_lag > max_lag)
    throw std::invalid_argument(
        part + ": minimum lag cannot be greater than maximum lag");
  if (min_lag < 0)
    throw std::invalid_argument(part + ": lags must be non-negative");
  if (vars.size() == 0)
    throw std::invalid_argument(part + ": no variable is included");
  for (string var : vars) {
    if (!variable_exists(var, cols))
      throw std::invalid_argument(part + ": " + var + " does not exist");
    if (variable_exists(var, DGMM_list.names))
      throw std::invalid_argument(
          part + ": " + var +
          " cannot be declared as GMM variable for twice or more");
    vector<int> ivec(max_lag - min_lag + 1);
    std::iota(ivec.begin(), ivec.end(), min_lag);
    DGMM_list.append(var, ivec, false, auto_mode, cols);
    int Lmin_lag;
    if (min_lag - 1 < 0)
      Lmin_lag = 0;
    else
      Lmin_lag = min_lag - 1;
    vector<int> Livec(max_lag - Lmin_lag + 1);
    std::iota(Livec.begin(), Livec.end(), Lmin_lag);
    LGMM_list.append(var, Livec, false, auto_mode, cols);
  }
}

void Command::parse_IV(vector<string> &matching_parts, string gmm) {
  smatch match;
  string temp = gmm;
  regex r("iv[(].{1,}[)]");
  while (regex_search(temp, match, r)) {
    string part = match.str(0);
    matching_parts.push_back(part);
    smatch match_groups_multiple;
    regex prog_1("^iv[(](.{1,})[)]$");
    regex_match(part, match_groups_multiple, prog_1);
    vector<string> vars = splitString(match_groups_multiple[1], ' ');
    parse_spaced_vars(vars, IV_list);
    temp = match.suffix().str();
  }
}

void Command::parse_options(string options_str) {
  vector<string> list_options = splitString(options_str, ' ');
  for (string option : list_options) {
    
    if (option == "onestep")
      options.steps = 1;
	else if (option =="constant")
		options.constant= true;
    //else if (option == "iterated")
    //  options.steps = 1000;
    else if (option == "nolevel")
      options.level = false;
    else if (option == "hqic")
      options.mmsc = "hqic";
    else if (option == "fod")
      options.transformation = "fod";
    else if (option == "timedumm")
      options.timedumm = true;
    else if (option == "collapse")
      options.collapse = true;
	else if (option == "oirf")
		options.irf="oirf";
    else
      throw std::invalid_argument(option + ": is not a valid option");
  }
  if (options.constant && !options.level) {
  	throw std::invalid_argument("Options constant and nolevel are mutually exclusive");
  }
}

void Command::check_dep() {
  string ret = Endo_list.purge();
  if (ret != "")
    throw std::invalid_argument(ret);

  for (std::size_t i = 0; i != Endo_list.names.size(); ++i) {
    string dep_name = Endo_list.names[i];
    vector<int> list_dep_lags = Endo_list.lags[i];
    if (list_dep_lags[0] != 0)
      throw std::invalid_argument("dependent variable should not be lagged on "
                                  "the left hand side of the model");
    if (list_dep_lags.size() == 0)
      throw std::invalid_argument(
          "lagged dependent variable should be included");
  }
}

void Command::check_exog() {
  string ret;
  if (Exog_list.names.size() >= 1) {
    ret = Exog_list.purge();
    if (ret != "")
      throw std::invalid_argument(ret);
  }
}

void Command::check_GMM() {
  vector<string> dep_names = Endo_list.names;

  for (std::size_t i = 0; i != DGMM_list.names.size(); ++i) {

    if (getIndex(dep_names,DGMM_list.names[i] )>=0) {
      if (DGMM_list.lags[i][0] < 2)
        throw std::invalid_argument("must use lag 2 or earlier of the "
                                    "dependent variable as instruments");
      // break;
    }
  }
}

void Command::check_iv() {
  string ret = IV_list.purge();
  if (ret != "")
    throw std::invalid_argument(ret);
}

void Command::check_three_lists() {
  for (auto iv_name : IV_list.names) {
    if (variable_exists(iv_name, Endo_list.names))
      throw std::invalid_argument("variable " + iv_name +
                                  ": a dependent variable cannot be instrument variable");

    if (variable_exists(iv_name, DGMM_list.names))
      throw std::invalid_argument("variable " + iv_name +
                                  ": a variable can be either in GMM style or "
                                  "in IV style, but not both");
                  
  }
  int i = 0;

  
  for (auto var_name : Exog_list.names) {
    vector<int> var_lags = Exog_list.lags[i];
    bool bool_GMM = variable_exists(var_name, DGMM_list.names);
    if (!bool_GMM) {
      bool bool_IV = variable_exists(var_name, IV_list.names);
      if (!bool_IV) {
        IV_list.append(var_name, var_lags, false, false, cols);
      }
    }
    i += 1;
  }
}
