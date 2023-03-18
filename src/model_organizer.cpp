#include "model_organizer.h"
#include <iostream>

vector<model> generate_list_models(int T) {
	vector<model> list_models;
	model first;
	list_models.push_back(first);
	list_to_Dep(T, list_models);
	list_to_exog_iv(T, list_models, Exog_list, 2);
	list_to_exog_iv(T, list_models, IV_list, 3);
	list_to_gmm(T, list_models, DGMM_list);
	for (std::size_t i = 0; i < list_models.size(); ++i) {
		list_models[i].num_dep = Endo_list.names.size();
		////cout<<"----update"<<endl;
		update_commandStr(list_models[i]);
		////cout << list_models[i].command_str<< endl;
	}

	// for (auto m : list_models){
	//   //cout<<"-------------------------"<<endl;
	//   //cout<<list_models.size()<<endl;
	//   //cout << m.dep_indep.size() << endl;
	//   //cout << m.Dgmm_vars.size() << endl;
	//   //cout << m.Lgmm_vars.size() << endl;
	// }

	////cout << " models: " << list_models.size() << endl;
	return list_models;
}

void update_commandStr(model &m) {
	m.command_str = "";
	for (gmm_var var: m.Dgmm_vars) {
		m.command_str += (" gmm(" + var.name + "," + std::to_string(var.min_lag) + ":" + std::to_string(var.max_lag) + ") ");
	}
	int num_iv = m.iv_vars.size();
	string prefix;
	for (int i = 0; i < num_iv; ++i) {
		prefix = " ";
		if (i == 0) {
			prefix = "";
			m.command_str += " iv(";
		}
		string tmp_str;
		for (regular_variable var: m.iv_vars) {
			if (var.lag > 0)
				tmp_str = "L" + std::to_string(var.lag) + "." + var.name + " ";
			else
				tmp_str = var.name;
			m.command_str += tmp_str;
		}
		m.command_str += ")";
	}
	////cout<<m.command_str<<endl;
	//m.command_str[2]="iv(";

}





// void list_to_gmm_m(int T, vector<model> &list_models,
//                    List_Variables &temp_list) {

//   for (std::size_t i = 0, max = temp_list.names.size(); i < max; ++i) {
//     string var_name = temp_list.names[i];
//     vector<int> lags = temp_list.lags[i];
//     int num_lags = lags.size();

//     int min_lag = lags[0];

//     int max_lag = min_lag + num_lags - 1;

//     for (std::size_t m = 0, max = list_models.size(); m < max; ++m) {
//       gmm_var new_var = gmm_var(var_name, min_lag, max_lag, 0);
//       list_models[m].Dgmm_vars.push_back(new_var);

//       int Lmin_lag = min_lag - 1;
//       if (Lmin_lag < 0)
//         Lmin_lag = 0;
//       gmm_var temp_var = gmm_var(var_name, Lmin_lag, min_lag,
//                                  0); // # the 3rd argument (min_lag) not used
//       list_models[m].Lgmm_vars.push_back(temp_var);
//     }
//   }
// }

void list_to_exog_iv(int T, vector<model> &list_models, List_Variables &the_list,
                     int dest) {
	int num_var = the_list.names.size();
	for (int i = 0; i < num_var; ++i) {
		string name = the_list.names[i];
		vector<int> lags = the_list.lags[i];
		int num_lags = lags.size();
		int min_lag = lags[0];
		int max_lag = min_lag + num_lags - 1;

		// navigate through models
		for (std::size_t m = 0, max = list_models.size(); m < max; ++m) {
			switch (dest) {
				case 2:
					for (int j = min_lag; j <= max_lag; ++j) {
						list_models[m].dep_indep.push_back(regular_variable(name, j));
					}
					break;
				case 3:
					for (int j = min_lag; j <= max_lag; ++j) {
						list_models[m].iv_vars.push_back(regular_variable(name, j));
					}
					break;
					// case 4:
					//   gmm_var new_var = gmm_var(var_name, min_lag, max_lag, 0);
					//   list_models[m].Dgmm_vars.push_back(new_var);

					//   int Lmin_lag = min_lag - 1;
					//   Lmin_lag = 0;
					//   gmm_var temp_var = gmm_var(name, Lmin_lag, min_lag,
					//                              0); // # the 3rd argument (min_lag) not used
					//   list_models[m].Lgmm_vars.push_back(temp_var);
					//   break;
			}
		}
		if (the_list.adjustable_max_lags[i])
			explode_exog_iv(max_lag, T, list_models, name, dest);
	}
}

// void list_to_exog_iv(int T, vector<model> &list_models,
//                      List_Variables &the_list, int dest) {
//   int num_var = the_list.names.size();

//   for (int i = 0; i < num_var; ++i) {
//     string name = the_list.names[i];
//     vector<int> lags = the_list.lags[i];
//     int num_lags = lags.size();
//     int min_lag = lags[0];
//     int max_lag = min_lag + num_lags - 1;
//     for (int j = min_lag; j <= max_lag; ++j) {
//       for (std::size_t m = 0, max = list_models.size(); m < max; ++m) {
//         if (dest != 3) {
//           list_models[m].dep_indep.push_back(regular_variable(name, j));
//           // if ((j==max_lag) && (dest==1))
//           //   list_models[m].num_dep_lags=max_lag;
//           if (j == 0)
//             list_models[m].command_str += (name + " ");
//           else
//             list_models[m].command_str +=
//                 ("L" + std::to_string(j) + "." + name + " ");
//         } else
//           list_models[m].iv_vars.push_back(regular_variable(name, j));
//       }
//     }
//   }
// }

void list_to_Dep(int T, vector<model> &list_models) {
	//int num_var = Endo_list.names.size();

	if (DEP_lags >= 1) {
		add_Dep(list_models[0], 0, DEP_lags);
		list_models[0].num_dep_lags = DEP_lags;
	} else {
		add_Dep(list_models[0], 0, 1);
		list_models[0].num_dep_lags = 1;
		vector<model> new_models;
		for (int last_lag = 2; last_lag <= T - 2; ++last_lag) {
			model new_m = list_models[0];
			add_Dep(new_m, 2, last_lag);
			new_m.num_dep_lags = last_lag;
			new_models.push_back(new_m);
		}
		list_models.insert(list_models.end(), new_models.begin(), new_models.end());

	}
}

void add_Dep(model &m, int first_lag, int last_lag) {
	int num_var = Endo_list.names.size();
	for (int lag = first_lag; lag <= last_lag; ++lag) {
		for (int i = 0; i < num_var; ++i) {
			string name = Endo_list.names[i];
			m.dep_indep.push_back(regular_variable(name, lag));
		}
	}

}

// void list_to_dep_indep_iv_m(int T, vector<model>& list_models, List_Variables
// &the_list, int dest) {
//   int num_var = the_list.names.size();

//   for (int i = 0; i < num_var; ++i) {
//     string name = the_list.names[i];
//     vector<int> lags = the_list.lags[i];
//     int num_lags = lags.size();
//     int min_lag = lags[0];
//     int max_lag = min_lag + num_lags - 1;
//     for (int j = min_lag; j <= max_lag; ++j) {
//       for (std::size_t m=0, max=list_models.size(); m<max; ++m) {
//         if (dest!=3){
//           list_models[m].dep_indep.push_back(regular_variable(name, j));
//           if ((j==max_lag) && (dest==1))
//             list_models[m].num_dep_lags=max_lag;
//           if (j==0)
//             list_models[m].command_str += (name + " ");
//           else
//             list_models[m].command_str += ("L" + std::to_string(j)+ "."+ name
//             + " ");
//         }
//         else
//           list_models[m].iv_vars.push_back(regular_variable(name, j));

//       }
//     }

//     if (the_list.adjustable_max_lags[i]) {
//       vector<model> new_models;
//       for (std::size_t m=0, max=list_models.size(); m<max; ++m) {
//         for (int j = max_lag + 1; j <= T - 1; ++j) {
//           model new_m = list_models[m];
//           for (int k = max_lag + 1; k < j + 1; ++k) {
//             if (dest != 3){
//               new_m.dep_indep.push_back(regular_variable(name, k));
//               if((k==j) && (dest==1))
//                 new_m.num_dep_lags=j;
//               if (k==0)
//                 new_m.command_str += (name + " ");
//               else
//                 new_m.command_str += ("L" + std::to_string(k)+ "."+ name + "
//                 ");
//             }
//             else
//               new_m.iv_vars.push_back(regular_variable(name, k));
//           }
//           new_models.push_back(new_m);
//         }
//       }

//       list_models.insert(list_models.end(), new_models.begin(),
//                          new_models.end());
//     }
//   }

// }

//

void list_to_gmm(int T, vector<model> &list_models, List_Variables &the_list) {
	int num_var = the_list.names.size();
	for (int i = 0; i < num_var; ++i) {
		string var_name = the_list.names[i];
		vector<int> lags = the_list.lags[i];
		int min_lag = lags[0];
		if (the_list.adjustable_max_lags[i]) {
			////cout << "explo" << endl;
			explode_gmm(min_lag, T, list_models, var_name);

		} else {
			int num_lags = lags.size();
			int min_lag = lags[0];
			int max_lag = min_lag + num_lags - 1;
			for (std::size_t m = 0, max = list_models.size(); m < max; ++m)
				add_gmm(list_models[m], var_name, min_lag, max_lag);
		}
	}
}

void add_gmm(model &m, string var_name, int min_lag, int max_lag) {
	gmm_var new_var = gmm_var(var_name, min_lag, max_lag, 0);
	m.Dgmm_vars.push_back(new_var);
	int Lmin_lag = min_lag - 1;
	if (Lmin_lag < 0)
		Lmin_lag = 0;
	gmm_var temp_var = gmm_var(var_name, Lmin_lag, min_lag,
	                           0); // # the 3rd argument (min_lag) not used
	m.Lgmm_vars.push_back(temp_var);
}

void explode_gmm(int min_lag, int T, vector<model> &list_models,
                 string var_name) {
	vector<model> new_models;
	for (std::size_t m = 0, max = list_models.size(); m < max; ++m) {
		for (int j = min_lag + 1; j <= T; ++j) {
			model new_m = list_models[m];
			add_gmm(new_m, var_name, min_lag, j);
			new_models.push_back(new_m);
		}
		add_gmm(list_models[m], var_name, min_lag, min_lag);
	}
	list_models.insert(list_models.end(), new_models.begin(), new_models.end());
}

void explode_exog_iv(int max_lag, int T, vector<model> &list_models,
                     string var_name, int dest) {
	vector<model> new_models;
	for (std::size_t m = 0, max = list_models.size(); m < max; ++m) {
		for (int j = max_lag + 1; j <= T - 1; ++j) {
			model new_m = list_models[m];
			for (int k = max_lag + 1; k < j + 1; ++k) {
				if (dest == 2) // exog
					new_m.dep_indep.push_back(regular_variable(var_name, k));
				else if (dest == 3) // iv
					new_m.iv_vars.push_back(regular_variable(var_name, k));
			}
			new_models.push_back(new_m);
		}
	}
	list_models.insert(list_models.end(), new_models.begin(), new_models.end());
}
