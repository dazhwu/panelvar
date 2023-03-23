//
// Created by Tiger on 5/4/2022.
//

#include "Common_Functions.h"

#define EIGEN_INITIALIZE_MATRICES_BY_NAN
using namespace std;
using namespace Eigen;
using std::string;
using std::vector;

vector<int>
gen_random_draws(int total_num_draws, int from, int to) {
    // generate a vector of integers that are uniformly distributed between "from" and "to"

    std::default_random_engine gen(123);
    std::uniform_int_distribution<int> distrib(from, to);

    vector<int> tbr(total_num_draws);
//#pragma omp parallel for
    for (int n = 0; n < total_num_draws; ++n)
        tbr[n] = (distrib(gen));

    return (tbr);
}

vector<bool>
row_has_nan(const RowMatrixXd &x) {
    int num_rows = x.rows();
    vector<bool> tbr(num_rows);
    for (int i = 0; i < num_rows; ++i) {
        tbr[i] = !((x.row(i).array() == x.row(i).array()).all());
    }
    return tbr;
}

RowMatrixXd
common_inv(Ref<RowMatrixXd> ori) {
    Eigen::FullPivLU<RowMatrixXd> Lu(ori);
    return ori.completeOrthogonalDecomposition().pseudoInverse();
    // if (Lu.isInvertible())
    //      return ori.inverse();
    //   else
    //        return ori.completeOrthogonalDecomposition().pseudoInverse();
}

double
standard_normalCDF(double x)   // Phi(-∞, x) aka N(x)
{
    return std::erfc(-x / std::sqrt(2)) / 2;
}

/* The function checks if the array elements are consecutive
If elements are consecutive, then returns true, else returns
false */
bool
areConsecutive(int *arr, int n)   // int arr[]
{
    // int n = sizeof(arr) / sizeof(arr[0]);
    if (n < 1)
        return false;
    int min = *min_element(arr, arr + n);

    int max = *max_element(arr, arr + n);
    if (max - min + 1 == n) {
        /* Create a temp array to hold visited flag of all elements.
        Note that, calloc is used here so that all values are initialized
        as false */
        bool *visited = (bool *) calloc(n, sizeof(bool));
        int i;
        for (i = 0; i < n; i++) {
            /* If we see an element again, then return false */
            if (visited[arr[i] - min] != false)
                return false;

            /* If visited first time, then mark the element as visited */
            visited[arr[i] - min] = true;
        }

        /* If all elements occur once, then return true */
        return true;
    }
    return false;   // if (max - min  + 1 != n)
}

// void tokenize(std::string const &str, const char* delim,
//             std::vector<std::string> &out)
//{
//     char *token = strtok_s(const_cast<char*>(str.c_str()), delim);
//     while (token != nullptr)
//     {
//         out.push_back(std::string(token));
//         token = strtok_s(nullptr, delim);
//     }
// }

vector<string>
splitString(string str, char splitter) {
    vector<string> result;
    string current = "";
    for (std::size_t i = 0, max = str.size(); i < max; i++) {
        if (str[i] == splitter) {
            if (current != "") {
                result.push_back(current);
                current = "";
            }
            continue;
        }
        current += str[i];
    }
    if (current.size() != 0)
        result.push_back(current);
    return result;
}

bool
variable_exists(string var_name, vector<string> &cols) {
    // for (string var_name : names)
    //	if (var_name == var_name)
    //		return true;

    // return false;

    // if (std::find(cols.begin(), cols.end(), var_name) == cols.end())
    // 	return false;
    // return true;
    if (getIndex(cols, var_name) >= 0)
        return true;
    else
        return false;
}

int
getIndex(vector<string> &v, string K) {
    auto it = std::find(v.begin(), v.end(), K);

    // If element was found
    if (it != v.end()) {
        int index = it - v.begin();
        return index;
    } else {
        // If the element is not
        // present in the vector
        return -1;
    }
}

std::vector<string>
combined_vector(std::vector<string> v1, std::vector<string> v2) {
    std::vector<string> combined;   //(v1.size()+v2.size());
    // std:://cout<<"1"<<std::endl;
    combined.reserve(v1.size() + v2.size());
    combined.insert(combined.end(), v1.begin(), v1.end());
    combined.insert(combined.end(), v2.begin(), v2.end());
    sort(combined.begin(), combined.end());
    auto last = unique(combined.begin(), combined.end());
    // v now holds {1 2 3 4 5 6 7 x x x x x x}, where 'x' is indeterminate
    combined.erase(last, combined.end());
    return combined;
}

void
lag(Ref<RowMatrixXd> mat, Ref<RowMatrixXd> lagged, int N, int lag_number, double fill) {
    int height = mat.rows() / N;
    int width = mat.cols();
#pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        // start_row = i * height;
        // end_row = start_row + height;
        Ref<RowMatrixXd> mat_i = mat.block(i * height, 0, height, width);
        Ref<RowMatrixXd> lagged_i = lagged.block(i * height, 0, height, width);

        // if (!isnan(fill))
        //   lagged_i.block(0, 0, lag_number, width) =
        //       RowMatrixXd::Zero(lag_number, width);
        // else
        lagged_i.block(0, 0, lag_number, width) = RowMatrixXd::Constant(lag_number, width, fill);
        lagged_i.block(lag_number, 0, height - lag_number, width) = mat_i.block(0, 0, height - lag_number, width);
    }
}

RowMatrixXd
get_first_diff_table(Ref<RowMatrixXd> ori_arr, int N) {
    int num_cols = ori_arr.cols();
    int num_rows = ori_arr.rows();
    // int height = num_rows / N;

    RowMatrixXd lag_arr(num_rows, num_cols);
    RowMatrixXd tbr_arr(num_rows, num_cols);
    lag(ori_arr, lag_arr, N, 1, NAN);
    // //saveData("ori.csv", ori_arr);
    //  //saveData("lag.csv", lag_arr);
    tbr_arr = ori_arr - lag_arr;
    //  //saveData("tbr.csv", tbr_arr);
    return tbr_arr;
}

MatrixXi
is_dep_NAs(Ref<RowMatrixXd> deps) {
    int num_cols = deps.cols();
    int num_rows = deps.rows();
    MatrixXi tbr = MatrixXi::Zero(num_rows, num_cols);
    for (int i = 0; i < num_rows; ++i) {
        for (int j = 0; j < num_cols; ++j) {
            if (isnan(deps(i, j)))
                tbr(i, j) = 1;
        }
    }
    return tbr;
}

RowMatrixXd
get_fod_table(Ref<RowMatrixXd> ori_arr, Ref<RowMatrixXd> ori_dep, int N, int type, int num_dep, int num_dep_lags) {
    int num_cols = ori_arr.cols();
    int num_rows = ori_arr.rows();
    int height = num_rows / N;
    MatrixXi dep_NAs = MatrixXi::Zero(ori_dep.rows(), ori_dep.cols());
    if (type == 1) {   // 0 if dep   1 if indep
        dep_NAs = is_dep_NAs(ori_dep);
    }
    RowMatrixXd tbr = RowMatrixXd::Constant(num_rows, num_cols, NAN);
    int this_count;
    for (int col_index = 0; col_index < num_cols; ++col_index) {
        vector<int> is_na(num_rows, 0);
        if ((type == 1) && (col_index < num_dep * num_dep_lags)) {
            int which_dep = col_index % num_dep;
            for (int i = 0; i < num_rows; ++i)
                is_na[i] = dep_NAs(i, which_dep);
        }
        for (int i = 0; i < N; ++i) {
            Ref<RowMatrixXd> ori_i = ori_arr.block(i * height, col_index, height, 1);
            Ref<RowMatrixXd> tbr_i = tbr.block(i * height, col_index, height, 1);

            // RowMatrixXd temp = RowMatrixXd::Constant(height, num_cols, NAN);
            double next_sum = NAN;
            int next_count = 0;
            double this_sum;
            // for j in range(height - 2, -1, -1):

            for (int j = height - 2; j >= 0; --j) {
                if (isnan(ori_i(j + 1, 0)) || (is_na[j + 1 + (i * height)] == 1)) {
                    this_count = next_count;
                    this_sum = next_sum;
                    if (j < height - 2)
                        tbr_i(j + 1, 0) = tbr_i(j + 2, 0);
                } else {
                    this_count = next_count + 1;
                    double temp_next_sum, temp_ori;
                    if (isnan(next_sum))
                        temp_next_sum = 0;
                    else
                        temp_next_sum = next_sum;
                    if (isnan(ori_i(j + 1, 0)))
                        temp_ori = 0;
                    else
                        temp_ori = ori_i(j + 1, 0);
                    this_sum = temp_next_sum + temp_ori;

                    // std:://cout << this_sum << std::endl;
                    // this_avg = this_sum * (1.0 / this_count);
                    tbr_i(j + 1, 0) =
                        (ori_i(j, 0) - this_sum * (1.0 / this_count)) * sqrt(this_count * 1.0 / (this_count + 1));
                }
                next_sum = this_sum;
                next_count = this_count;
            }
            tbr_i(0, 0) = NAN;
        }
    }
    return tbr;
}

void
saveRef2CSV(string fileName, Ref<RowMatrixXd> matrix) {
    const static IOFormat CSVFormat(FullPrecision, DontAlignCols, ", ", "\n");

    ofstream file(fileName);
    if (file.is_open()) {
        file << matrix.format(CSVFormat);
        file.close();
    }
}

void
saveData(string fileName, RowMatrixXd &matrix) {
    // https://eigen.tuxfamily.org/dox/structEigen_1_1IOFormat.html
    const static IOFormat CSVFormat(FullPrecision, DontAlignCols, ", ", "\n");

    ofstream file(fileName);
    if (file.is_open()) {
        file << matrix.format(CSVFormat);
        file.close();
    }
}
