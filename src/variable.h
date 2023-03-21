

#ifndef PYDYNPD_VARIABLE_H
#define PYDYNPD_VARIABLE_H

#include "List_Variables.h"
#include "pvar.h"

class regular_variable {
  public:
    string name;
    int lag;
    regular_variable(string, int);
};

class gmm_var {
  public:
    string name;
    int min_lag;
    int max_lag;
    int lag;
    gmm_var(string, int, int, int);
};

#endif   // PYDYNPD_VARIABLE_H
