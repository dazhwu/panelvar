import panelvar.pvar_module as pvar_module

import warnings
import numpy as np
from pandas import DataFrame

from panelvar.model_summary import model_summary
from panelvar.panel_data import panel_data
from matplotlib import pyplot as plt
from sys import exit

import sys

warnings.filterwarnings("ignore", category=RuntimeWarning)


class pvar:

    def __init__(self, dep, df: DataFrame, identifiers: list, exog="", lags=1, gmm="", options="", ahead=10, draw=100):

        if len(identifiers) != 2:
            print('two variables needed')
            exit()

        if len(dep.strip()) == 0:
            print('no dependent variable specified')
            exit()

        if lags < 1:
            print('argument lags needs to be at least 1')
            exit()

        pdata = panel_data(df, identifiers)
        try:
            (names, model_options) = pvar_module.process_command(pdata.T, dep, lags, exog, gmm, options, df.columns)
        except Exception as e:
            print(e)
            sys.exit(1)

        pdata.export_data(df, names, model_options.timedumm)

        self.list_models = pvar_module.prepare_data(identifiers, pdata.data, [pdata.N, pdata.T], pdata.cols,
                                                    pdata.col_timedumm,
                                                    ahead, draw)

        self._good_models = []
        self._bad_models = []

        for m in self.list_models:

            # new_model = dynamic_panel_model(identifiers, m.dep_indep, m.regression_result, m.model_info, m.hansen,
            #                                 m.stability, m.irf, m.model_options, m.command_str)

            if (self.check_model(m) == True):
                self._good_models.append(m)
            else:
                self._bad_models.append(m)

        for m in self._good_models:
            self.form_results(m)
            print(m.command_str)
            m.plot_irf()

    def form_results(self, model):

        ms = model_summary()
        reg_tables = self.form_regression_table(model)
        ms.print_summary(model, reg_tables)

    def check_model(self, model):
        tbr = False

        if (max(model.stability) <= 1):
            if model.hansen.P_value > 0.05 and model.hansen.P_value < 0.99999:
                tbr = True

        return (tbr)


    def form_regression_table(self, model):
        var_names = model.model_info.indep
        regression_tables = {}
        i = 0
        for dep in model.model_info.dep:
            coeff = model.regression_result.beta[:, i]
            std_err = model.regression_result.std_err[:, i]
            z_value = model.regression_result.Z_values[:, i]
            p_value = model.regression_result.P_values[:, i]
            sig = ['***' if p <= 0.001 else ('**' if p <= 0.01 else ('*' if p <= 0.05 else ' ')) for p in p_value]
            regression_tables[dep] = DataFrame(list(zip(var_names, coeff, std_err, z_value, p_value, sig)),
                                               columns=['variable', 'coefficient', 'std_err', 'z_value', 'p_value',
                                                        'sig'])

            i = i + 1
        return (regression_tables)

    def plot_irf(self, model):
        print("plot_irf")

        ahead=model.irf[0].shape[0]
        num_dep=model.model_info.num_dep;
        plt.rcParams.update({'font.size': 22})
        x = (np.arange(0, ahead).reshape(ahead, 1))[:, 0]
        for i in range(0, num_dep):  #from
            for j in range(0, num_dep):  #to
                # fig = plt.figure()
                # ax = fig.add_subplot(1, 1, 1)
                y = model.irf[0][:, i*num_dep+j]
                l = model.irf[1][:, i*num_dep+j]
                u = model.irf[2][:, i*num_dep+j]
                print("image")
                fig, ax = plt.subplots()
                fig.set_figwidth(20)
                fig.set_figheight(10)
                ax.plot(x, y)
                ax.fill_between(x, l, u, color='b', alpha=.1)
                plt.grid(color = 'green', linestyle = '--', linewidth = 0.5)
                plt.ylabel(model.model_info.dep[i] + " on " +model.model_info.dep[j])

                # ax.scatter(x, y,color='b')
                # plt.show()
                fig.savefig(model.model_info.dep[i] + " on " +model.model_info.dep[j] + '.png',bbox_inches='tight')
