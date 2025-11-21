/*
Written by Leon Tolksdorf, 2025. Please refer to: 

@article{CollisionProbabilityEstimation,
  title   = {Collision Probability Estimation for Optimization-based Vehicular Motion Planning},
  author  = {Tolksdorf, Leon and Tejada, Arturo and Birkner, Christian and van de Wouw, Nathan},
  journal = {arXiv preprint arXiv:2505.21161} ,
  year    = {2025}
}
*/

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

#define DEBUG 0

void POCFunction::_sort_intervals() {

    MatrixXd all_int_sorted = _all_int;
    MatrixXd int_1, int_2, int1_new, int2_new;
    VectorXd size_int1, size_int2;
    int pp, n_count = 0, p = 0;

    #if DEBUG
        cout << "all intervals before sorting: \n" << all_int_sorted << endl;
    #endif

    for (int k = 0; k < _n_cols; k += 2) {
        p++;
        pp = 0;
        for (int i = 0; i < _n_cols; i += 2) {
            pp++;
            int_1 = all_int_sorted.block(0, k, _n_rows, 2); //grab one interval
            if (i != k) { //skip if I1 and I2 = I1 would be selected
                int_2 = all_int_sorted.block(0, i, _n_rows, 2); //grab another interval
                _int_unions(&int_1, &int_2, &int1_new, &int2_new);
                //update the intervals back into the matrix
                all_int_sorted.block(0, k, _n_rows, 2) = int1_new;
                all_int_sorted.block(0, i, _n_rows, 2) = int2_new;
                if (i + 2 > _n_cols) {
                    break;
                }
            }
        }
    }
    #if DEBUG
        cout << "all intervals after sorting: \n" << all_int_sorted << endl;
    #endif
    //output
    _thetas = all_int_sorted;
}