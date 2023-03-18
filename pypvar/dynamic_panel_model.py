#import math

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
#import scipy

#import pydynpd.gmm_module as gmm_module
#from pydynpd.info import df_info, options_info
#from pydynpd.panel_data import panel_data


class dynamic_panel_model(object):
    def __init__(self, identifiers, dep_indep, regression, basic_info, hansen, stability, irf, options, command_str):
        self.name = basic_info.model_name
        self.dep_indep=dep_indep

        self.indeps=basic_info.indep
        self.identifiers=identifiers
        self.T = basic_info.T
        self.N = basic_info.N
        self.actual_steps=basic_info.actual_steps
        self.num_obs=basic_info.num_obs
        self.num_indep=basic_info.num_indep
        self.num_dep=basic_info.num_dep
        self.num_dep_lags=basic_info.num_dep_lags
        self.num_instr=basic_info.num_instr
        self.dep=basic_info.dep
        #self.indep=basic_info.indep
        self.beta=regression.beta
        self.std_err=regression.std_err
        self.residual=regression.Residual
        self.vcov=regression.vcov
        self.form_regression_table(regression)
        self.model_options = options
        self.max_obs=basic_info.max_obs
        self.min_obs=basic_info.min_obs
        self.avg_obs=basic_info.avg_obs
        self.irf=irf[0]
        self.irf_L=irf[1]
        self.irf_U=irf[2]
        self.MMSC_LU = {}
        self.MMSC_LU["bic"] = basic_info.mmsc_lu.BIC
        self.MMSC_LU["hqic"] = basic_info.mmsc_lu.HQIC
        self.MMSC_LU["aic"] = basic_info.mmsc_lu.AIC

        self.hansen=hansen
        self.stability=abs(stability)
        self.command_str=command_str
        # self.command_str = list_parts[0] + '|' + list_parts[1]
        # if list_parts[2] != '':
        #     self.command_str += '|' + list_parts[2]
        
   

    def form_regression_table(self, regression):

        var_names=self.indeps
        self.regression_tables={}
        i=0
        for dep in self.dep:
            coeff = regression.beta[:, i]
            std_err = regression.std_err[:, i]
            z_value = regression.Z_values[:, i]
            p_value = regression.P_values[:, i]
            sig = ['***' if p <= 0.001 else ('**' if p <= 0.01 else ('*' if p <= 0.05 else ' ')) for p in p_value]
            self.regression_tables[dep] = pd.DataFrame(list(zip(var_names, coeff, std_err, z_value, p_value, sig)),
                                                 columns=['variable', 'coefficient', 'std_err', 'z_value', 'p_value',
                                                          'sig'])

            i=i+1


    def plot_irf(self):
        print("plot_irf")

        ahead=self.irf[0].shape[0]
        plt.rcParams.update({'font.size': 22})
        x = (np.arange(0, ahead).reshape(ahead, 1))[:, 0]
        for i in range(0, self.num_dep):  # matrix
            for j in range(0, self.num_dep):  # column
                # fig = plt.figure()
                # ax = fig.add_subplot(1, 1, 1)
                y = self.irf[i][:, j]
                l = self.irf_L[i][:, j]
                u = self.irf_U[i][:, j]
                print("image")
                fig, ax = plt.subplots()
                fig.set_figwidth(20)
                fig.set_figheight(10)
                ax.plot(x, y)
                ax.fill_between(x, l, u, color='b', alpha=.1)
                plt.grid(color = 'green', linestyle = '--', linewidth = 0.5)
                plt.ylabel(self.dep[i] + " on " +self.dep[j])

                # ax.scatter(x, y,color='b')
                # plt.show()
                fig.savefig(self.dep[i] + "on" +self.dep[j] + '.png',bbox_inches='tight')



   