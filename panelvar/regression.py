
import panelvar.pvar_module as pvar_module

import warnings
from pandas import DataFrame
from panelvar.dynamic_panel_model import dynamic_panel_model

from panelvar.model_summary import model_summary
from panelvar.panel_data import panel_data

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

        self.list_models = pvar_module.prepare_data(identifiers, pdata.data, [pdata.N, pdata.T], pdata.cols, pdata.col_timedumm,
                                                    ahead, draw)

        self._good_models = []
        self._bad_models = []

        for m in self.list_models:

            new_model = dynamic_panel_model(identifiers, m.dep_indep, m.regression_result, m.model_info, m.hansen,
                                            m.stability, m.irf, m.model_options, m.command_str)
            print(new_model.identifiers)
            if (self.check_model(new_model) == True):
                self._good_models.append(new_model)
            else:
                self._bad_models.append(new_model)

        for m in self._good_models:
            self.form_results(m)
            print(m.command_str)
            m.plot_irf()

    def form_results(self, model):

        ms = model_summary()
        ms.print_summary(model)

    def check_model(self, model):
        tbr = False

        if (max(model.stability) <= 1):
            if model.hansen.P_value > 0.05 and model.hansen.P_value < 0.99999:
                tbr = True

        return (tbr)
