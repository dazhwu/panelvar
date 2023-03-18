import numpy as np

import panelvar.pvar_module as pvar_module

import time
import warnings

import panelvar.dynamic_panel_model as dynm
from pandas import DataFrame
from panelvar.dynamic_panel_model import dynamic_panel_model

from panelvar.model_summary import model_summary
from panelvar.panel_data import panel_data
from panelvar.irf import oirf
from sys import exit


import sys

warnings.filterwarnings("ignore", category=RuntimeWarning)





class pvar:

    def __init__(self, dep, df:DataFrame, identifiers:list, exog="", lags=1, gmm="", options="", ahead=10,draw=100):

        if len(identifiers) != 2:
            print('two variables needed')
            exit()

        if len(dep.strip()) == 0:
            print('no dependent variable specified')
            exit()

        # if lags < 1:
        #     print('argument lags needs to be at least 1')
        #     exit()
        pdata = panel_data(df, identifiers)
        try:
            (names, model_options) = pvar_module.process_command(pdata.T, dep, lags, exog, gmm, options, df.columns)
        except Exception as e:
            print(e)
            sys.exit(1)

        pdata.export_data(df, names, model_options.timedumm)

        self.list_models = pvar_module.prepare_data(pdata.data, [pdata.N, pdata.T], pdata.cols, pdata.col_timedumm, ahead, draw)



        self._good_models = []
        self._bad_models = []
        for m in self.list_models:
        #     #print(m.hansen.test_value)

            new_model=dynamic_panel_model(identifiers, m.dep_indep, m.regression_result,   m.model_info, m.hansen, m.stability, m.irf, m.model_options,  m.command_str)
            #print(oirf(new_model,3))
            #
            #self.form_results(new_model)

        #
        #     if (options.beginner==True) & (len(list_models)>=2):
            if (self.check_model(new_model)==True):

                self._good_models.append(new_model)
            else:
                self._bad_models.append(new_model)


        for m in self._good_models:
            self.form_results(m)
            print(m.command_str)
            m.plot_irf()
        # for m in self._bad_models:
        #     self.form_results(m)
        #     print(m.command_str)
        #     m.plot_irf()

        #
        # ms = model_summary()
        # if len(self._good_models) >= 2:
        #     ms.print_good_list(self._good_models, options.level, options.mmsc)
        #
        # if len(self._bad_models) >= 1:
        #     print('\nThe following model(s) did not pass specification tests:')
        #     ms.print_bad_list(self._bad_models)

    

 

    def form_results(self, model):
  
        # if model.name != '':
        #     print(' ' + model.name)

 
        ms = model_summary()
        ms.print_summary(model)

        #self.models.append(model)  # results = {}

    def check_model(self, model):
        tbr = False

        if (max(model.stability) <= 1):
            if model.hansen.P_value > 0.05 and model.hansen.P_value < 0.99999:
                tbr=True

        return (tbr)