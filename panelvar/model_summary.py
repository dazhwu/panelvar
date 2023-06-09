import os

from prettytable import PrettyTable


class model_list(object):
    def __init__(self, m_list: list):
        self.names = []
        self.min_lags = []
        self.max_lags = []
        for m in m_list:
            for i in range(1, len(m.dep_indep)):
                indep = m.dep_indep[i]
                if not indep.name in self.names:
                    self.names.append(indep.name)
                    self.min_lags.append(indep.lag)
                    self.max_lags.append(indep.lag)
                else:
                    the_index = self.names.index(indep.name)
                    if indep.lag > self.max_lags[the_index]:
                        self.max_lags[the_index] = indep.lag

                    if indep.lag < self.min_lags[the_index]:
                        self.min_lags[the_index] = indep.lag


class model_summary(object):

    def print_summary(self, model, reg_tables):

        if model.model_options.steps == 2:
            str_steps = 'two-step '
        elif model.model_options.steps == 1:
            str_steps = 'one-step '
        else:
            
            str_steps = str(model.model_info.actual_steps) + '-step '

        if model.model_options.level:
            str_gmm = 'system GMM'
        else:
            str_gmm = 'difference GMM'

        to_print = []
        # to_print.append(model.command_str)
        to_print.append(' Dynamic Panel VAR Estimation, ' + str_steps + str_gmm)
        to_print.append(self.basic_information(model))
        to_print.append(self.regression_table(model, reg_tables))
        to_print.append(self.test_results(model))
        for line in to_print:
            print(line)

    def basic_information(self, model):

        basic_table = PrettyTable()
        middle_space = '                   '
        basic_table.field_names = ["    ", "   ", "  "]
        basic_table.border = False
        basic_table.header = False
        basic_table.align = 'l'
        basic_table.add_row(
            ['Group variable: ' + model.model_info.identifiers[0], middle_space, 'Number of obs = ' + str(model.model_info.num_obs)])
        basic_table.add_row(
            ['Time variable: ' + model.model_info.identifiers[1], middle_space, 'Min obs per group: ' + str(model.model_info.min_obs)])
        basic_table.add_row(['Number of instruments = ' + str(model.model_info.num_instr), middle_space,
                             'Max obs per group: ' + str(model.model_info.max_obs)])
        basic_table.add_row(
            ['Number of groups = ' + str(model.model_info.N), middle_space,
             'Avg obs per group: ' + '{0:.2f}'.format(model.model_info.avg_obs)])

        return (basic_table.get_string())

    def test_results(self, model):
        
        str_toprint = 'Hansen test of overid. restrictions: chi(' + str(model.hansen.df) + ') = ' + '{:.3f}'.format(
            model.hansen.test_value)
        str_toprint = str_toprint + ' Prob > Chi2 = ' + '{:.3f}'.format(model.hansen.P_value) + '\n'
        if(max(model.stability)<=1):
            str_toprint = str_toprint + 'All the eigenvalues lie inside the unit circle.' + '\n'
            str_toprint = str_toprint + 'PVAR satisfies stability condition.' +'\n'
        else:
            str_toprint = str_toprint + 'Not all the eigenvalues lie inside the unit circle.' + '\n'
            str_toprint = str_toprint + 'PVAR does not satisfy stability condition.' +'\n'


        return (str_toprint)

    def regression_table(self, model,regression_tables):
        r_table = PrettyTable()

        var_names=model.model_info.indep

        i = 0
        r_table.field_names = ["Equation","   ", "coef.", "Corrected Std. Err.", "z", "P>|z|", " "]
        r_table.float_format = '.7'
        for dep in model.model_info.dep:
            r_table.add_row([dep, "","","","","",""])
            reg_table=regression_tables[dep]
            for j in range(len(model.model_info.indep)):
                var_name = reg_table['variable'][j]
                coeff = reg_table['coefficient'][j]
                stderr = reg_table['std_err'][j]
    
                z = reg_table['z_value'][j]
                p = reg_table['p_value'][j]
                sig = reg_table['sig'][j]
                r_table.add_row([" ", var_name, coeff, stderr, z, p, sig])

        return r_table.get_string()

        
    def print_good_list(self, the_list: list, level: bool, mmsc: str):
        the_list.sort(key=lambda x: x.MMSC_LU[mmsc])
        m_list = model_list(the_list)

        r_table = PrettyTable()
        # col_names=['variables'] + [m.name for m in the_list]
        variable_names = []
        for i in range(len(m_list.names)):
            v = m_list.names[i]
            for j in range(m_list.min_lags[i], m_list.max_lags[i] + 1):
                if j == 0:
                    variable_names.append(v)
                else:
                    variable_names.append('L' + str(j) + '.' + v)
        if level:
            variable_names.append('_con')
        r_table.add_column('variables', variable_names)
        for m in the_list:
            new_col = []
            for i in range(len(variable_names)):
                new_col.append('')
            rt = m.regression_table
            ind = [variable_names.index(v) for v in rt['variable']]
            j = 0
            for i in ind:
                new_col[i] = '{:.3f}'.format(rt['coefficient'][j]) + rt['sig'][j] + '\n(' + '{:.3f}'.format(
                    rt['std_err'][j]) + ')'
                j += 1

            r_table.add_column(m.name, new_col)

        print('models are sorted by ' + mmsc)
        print(r_table)

        try:
            with open('output.html', 'w') as f:
                # print(f.__dir__())
                print('HTML output named "output.html" is located in folder ' + os.getcwd())
                f.write(r_table.get_html_string(
                    attributes={'border': 1, 'style': 'border-width: 1px; border-collapse: collapse;'}))
        except Exception as e:
            print(e)

        print('\n')
        print('MMSC_LU scores:')
        mmsc_table = PrettyTable()
        mmsc_table.field_names = ['model', 'aic', 'bic', 'hqic', 'command str']
        mmsc_table.align['command str'] = "l"

        for m in the_list:
            mmsc_table.add_row([m.name, '{:.3f}'.format(m.MMSC_LU['aic']), '{:.3f}'.format(m.MMSC_LU['bic']),
                                '{:.3f}'.format(m.MMSC_LU['hqic']), m.command_str])

        print(mmsc_table)

    def print_bad_list(self, the_list: list):

        bad_table = PrettyTable()
        bad_table.field_names = ['model', 'command str']
        bad_table.align['command str'] = "l"
        for m in the_list:
            bad_table.add_row([m.name, m.command_str])

        print(bad_table)
