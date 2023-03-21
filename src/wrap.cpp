#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "Command.h"
#include "Common_Functions.h"
#include "GMM.h"
#include "List_Variables.h"
#include "dynamic_model.h"
#include "info.h"
#include "instruments.h"
#include "pvar.h"
#include "variable.h"

namespace py = pybind11;

PYBIND11_MODULE(pvar_module, m) {

    using RowMatrixXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

    py::class_<model_options>(m, "model_options")
        .def(py::init<>())
        .def_readwrite("steps", &model_options::steps)
        .def_readwrite("constant", &model_options::constant)
        .def_readwrite("level", &model_options::level)
        .def_readwrite("beginner", &model_options::beginner)
        .def_readwrite("timedumm", &model_options::timedumm)
        .def_readwrite("collapse", &model_options::collapse)
        .def_readwrite("mmsc", &model_options::mmsc)
        .def_readwrite("irf", &model_options::irf)
        .def_readwrite("transformation", &model_options::transformation);

    py::class_<Regression>(m, "Regression")
        .def(py::init<>())
        .def(py::init<RowMatrixXd &, RowMatrixXd &, RowMatrixXd &, RowMatrixXd &, Ref<RowMatrixXd>, Ref<RowMatrixXd>,
                      Ref<RowMatrixXd>, RowMatrixXd &>())
        .def_readwrite("beta", &Regression::beta)
        .def_readwrite("vcov", &Regression::vcov)
        .def_readwrite("std_err", &Regression::std_err)
        .def_readwrite("Z_values", &Regression::Z_values)
        .def_readwrite("P_values", &Regression::P_values)
        .def_readwrite("weighting", &Regression::weighting)
        .def_readwrite("Instruments", &Regression::Instruments)
        .def_readwrite("Y", &Regression::Cy)
        .def_readwrite("X", &Regression::Cx)
        .def_readwrite("Residual", &Regression::Residual);

    py::class_<basic_info>(m, "basic_info")
        .def(py::init<>())
        .def(py::init<int, int, int, int, int, int, int, int, int, int, int, double, vector<regular_variable> &,
                      vector<string> &, model_options &>())
        .def_readwrite("model_name", &basic_info::model_name)
        .def_readwrite("identifiers", &basic_info::identifiers)
        .def_readwrite("dep", &basic_info::dep)
        .def_readwrite("indep", &basic_info::indep)
        .def_readwrite("N", &basic_info::N)
        .def_readwrite("T", &basic_info::T)
        .def_readwrite("actual_steps", &basic_info::actual_steps)
        .def_readwrite("num_obs", &basic_info::num_obs)
        .def_readwrite("num_instr", &basic_info::num_instr)
        .def_readwrite("num_dep", &basic_info::num_dep)
        .def_readwrite("num_dep_lags", &basic_info::num_dep_lags)
        .def_readwrite("num_indep", &basic_info::num_indep)
        .def_readwrite("diff_width", &basic_info::diff_width)
        .def_readwrite("x_height", &basic_info::x_height)
        .def_readwrite("max_obs", &basic_info::max_obs)
        .def_readwrite("min_obs", &basic_info::min_obs)
        .def_readwrite("avg_obs", &basic_info::avg_obs)
        .def_readwrite("SS", &basic_info::SS)
        .def_readwrite("mmsc_lu", &basic_info::mmsc_lu);

    py::class_<z_info>(m, "z_info")
        .def(py::init<>())
        .def(py::init<int, int, int, int, int, int, int, int>())
        .def_readwrite("diff_width", &z_info::diff_width)
        .def_readwrite("diff_height", &z_info::diff_height)
        .def_readwrite("level_width", &z_info::level_width)
        .def_readwrite("level_height", &z_info::level_height)
        .def_readwrite("z_width", &z_info::z_width)
        .def_readwrite("z_height", &z_info::z_height)

        .def_readwrite("num_Dgmm_instr", &z_info::num_Dgmm_instr)
        .def_readwrite("num_Lgmm_instr", &z_info::num_Lgmm_instr)
        .def_readwrite("num_instr", &z_info::num_instr);

    py::class_<List_Variables>(m, "List_Variables")
        .def(py::init<>())
        .def_readwrite("names", &List_Variables::names)
        .def_readwrite("lags", &List_Variables::lags)
        .def_readwrite("min_lags", &List_Variables::min_lags)
        .def_readwrite("max_lags", &List_Variables::max_lags)
        .def_readwrite("adjustable_min_lags", &List_Variables::adjustable_min_lags)
        .def_readwrite("adjustable_max_lags", &List_Variables::adjustable_max_lags);

    py::class_<Step_Result>(m, "Step_Result")
        .def(py::init<RowMatrixXd &, RowMatrixXd &, RowMatrixXd &, RowMatrixXd &, RowMatrixXd &, RowMatrixXd &,
                      RowMatrixXd &, RowMatrixXd &, RowMatrixXd &, RowMatrixXd &, RowMatrixXd &, RowMatrixXd &,
                      RowMatrixXd &>())
        .def_readwrite("residual", &Step_Result::residual)
        .def_readwrite("_residual_t", &Step_Result::_residual_t)
        .def_readwrite("XZ_W", &Step_Result::XZ_W)
        .def_readwrite("W", &Step_Result::W)
        .def_readwrite("W_inv", &Step_Result::W_inv)
        .def_readwrite("W_next", &Step_Result::W_next)
        .def_readwrite("_M_XZ_W", &Step_Result::_M_XZ_W)
        .def_readwrite("zs", &Step_Result::zs)
        .def_readwrite("ZuuZ", &Step_Result::ZuuZ)
        .def_readwrite("vcov", &Step_Result::vcov)
        .def_readwrite("M", &Step_Result::M)
        .def_readwrite("_zs_list", &Step_Result::_zs_list)
        .def_readwrite("beta", &Step_Result::beta)
        //.def_readwrite("SS", &Step_Result::SS)
        .def_readwrite("std_err", &Step_Result::std_err);

    // py::class_<AR_test_info>(m, "AR_test_info")
    //     .def(py::init<>())
    //     .def(py::init<int, double, double>())
    //     .def_readwrite("lag", &AR_test_info::lag)
    //     .def_readwrite("AR", &AR_test_info::AR)
    //     .def_readwrite("P_value", &AR_test_info::P_value);

    py::class_<Hansen_test_info>(m, "Hansen_test_info")
        .def(py::init<>())
        .def(py::init<double, int, double, double>())
        .def_readwrite("test_value", &Hansen_test_info::test_value)
        .def_readwrite("df", &Hansen_test_info::df)
        .def_readwrite("critical_value", &Hansen_test_info::critical_value)
        .def_readwrite("P_value", &Hansen_test_info::P_value);

    py::class_<df_info>(m, "df_info")
        .def(py::init<>())
        .def(py::init<int, int, int, int, int, int, int>())
        .def_readwrite("N", &df_info::N)
        .def_readwrite("T", &df_info::T)
        //.def_readwrite("ids", &df_info::ids)
        .def_readwrite("first_diff_index", &df_info::first_diff_index)
        .def_readwrite("last_diff_index", &df_info::last_diff_index)
        .def_readwrite("first_level_index", &df_info::first_level_index)
        .def_readwrite("last_level_index", &df_info::last_level_index)
        .def_readwrite("max_lag", &df_info::max_lag);

    py::class_<regular_variable>(m, "regular_variable")
        //.def(py::init<string, int>())
        .def_readwrite("name", &regular_variable::name)
        .def_readwrite("lag", &regular_variable::lag);

    py::class_<gmm_var>(m, "gmm_var")
        //.def(py::init<string, int, int, int>())
        .def_readwrite("name", &gmm_var::name)
        .def_readwrite("min_lag", &gmm_var::min_lag)
        .def_readwrite("max_lag", &gmm_var::max_lag)
        .def_readwrite("lag", &gmm_var::lag);

    py::class_<MMSC_LU>(m, "MMSC_LU")
        .def(py::init<>())
        .def(py::init<double, double, double>())
        .def_readwrite("BIC", &MMSC_LU::BIC)
        .def_readwrite("HQIC", &MMSC_LU::HQIC)
        .def_readwrite("AIC", &MMSC_LU::AIC);

    py::class_<returned_result>(m, "returned_result")
        .def(py::init<struct Regression &, Hansen_test_info &, VectorXcd &, struct basic_info &, struct model_options &,
                      string>())
        .def_readwrite("dep_indep", &returned_result::dep_indep)
        .def_readwrite("regression_result", &returned_result::regression_result)
        .def_readwrite("hansen", &returned_result::hansen)
        .def_readwrite("stability", &returned_result::stability)
        //.def_readwrite("AR_list", &returned_result::AR_list)
        .def_readwrite("model_info", &returned_result::model_info)
        .def_readwrite("model_options", &returned_result::model_options)
        .def_readwrite("irf", &returned_result::irf)
        .def_readwrite("command_str", &returned_result::command_str);

    // m.def("AR_test", &AR_test);
    m.def("process_command", &process_command);
    m.def("regular_process", &regular_process);
    m.def("get_z_table", &get_z_table);
    m.def("prepare_data", &prepare_data);
}

/*
g++  -O3 -Wall -shared -std=c++17 -fopenmp  -mavx2 `python3 -m pybind11 --includes` -I/usr/local/include
-I/usr/include/eigen3 `python3-config --ldflags` -o gmm_module.so GMM.cpp info.cpp Specification_Test.cpp
Step_Result.cpp Common_Functions.cpp wrap.cpp Command.cpp List_Variables.cpp instruments.cpp dynamic_model.cpp
variable.cpp



*/
