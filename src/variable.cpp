
#include "variable.h"

regular_variable::regular_variable(string _name, int _lag) {
    name = _name;
    lag = _lag;
}

gmm_var::gmm_var(string _name, int _min, int _max, int _lag) {
    name = _name;
    min_lag = _min;
    max_lag = _max;
    lag = _lag;
}
