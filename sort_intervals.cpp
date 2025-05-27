#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <vector>
#include <eigen3/Eigen/Dense>
#include <casadi/casadi.hpp>
#include "include/poc_function.hpp"

using namespace std;
using namespace Eigen;
using namespace casadi;

void POCFunction::_sort_intervals() {

    MatrixXd all_int_sorted = _all_int;

    //sort intervals
    MatrixXd int_sizes = MatrixXd::Zero(_n_rows, _n_cols * 0.5);
    MatrixXd TwoPiRow  = MatrixXd::Zero(_n_rows, _n_cols);
    TwoPiRow.col(1)    = VectorXd::Constant(_n_rows, 2.0 * M_PI);
    int p = 0;
    int pp;
    MatrixXd int_1;
    MatrixXd int_2;
    MatrixXd int1_new;
    MatrixXd int2_new;
    VectorXd size_int1;
    VectorXd size_int2;
    for (int k = 0; k < _n_cols; k += 2) {
        p++;
        pp = 0;
        for (int i = 0; i < _n_cols; i += 2) {
            pp++;
            int_1 = all_int_sorted.block(0, k, _n_rows, 2); //grab one interval
            size_int1 = int_1.col(1) - int_1.col(0);
            if (i != k) { //skip if I1 and I2 = I1 would be selected
                int_2 = all_int_sorted.block(0, i, _n_rows, 2); //interval at i and i+1
                size_int2 = int_2.col(1) - int_2.col(0);
                _int_unions(&int_1, &size_int1, &int_2, &size_int2, &int1_new, &int2_new);

                //update the sorted intervals back into the matrix
                all_int_sorted.block(0, k, _n_rows, 2) = int1_new;
                all_int_sorted.block(0, i, _n_rows, 2) = int2_new;
                if (i + 2 > _n_cols) {
                    break;
                }
            }
        }
    }
    
    //ensure the total measure of all intervals does not exceed 2pi
    int n_count = 0;
    for (int p = 0; p < _n_cols; p += 2) {
        VectorXd lb = all_int_sorted.col(p);
        VectorXd ub = all_int_sorted.col(p + 1);
        int_sizes.col(n_count) = ub - lb;
        n_count++;

        if (p + 2 > _n_cols) {
            break;
        }
    }

    //total interval sizes for each row
    VectorXd total_int_sizes = int_sizes.rowwise().sum();
    for (int i = 0; i < _n_rows; i++){
        //set rows where total interval sizes exceed 2pi
        if (total_int_sizes(i) > 2.0 * M_PI) {
            all_int_sorted.row(i) = TwoPiRow.row(i);
        }
    }

    //output
    _thetas             = all_int_sorted;
}