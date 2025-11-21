/*
Written by Leon Tolksdorf, 2025. Please refer to: 

@article{CollisionProbabilityEstimation,
  title   = {Collision Probability Estimation for Optimization-based Vehicular Motion Planning},
  author  = {Tolksdorf, Leon and Tejada, Arturo and Birkner, Christian and van de Wouw, Nathan},
  journal = {arXiv preprint arXiv:2505.21161} ,
  year    = {2025}
}
*/

#include <iostream>
#include <vector>
#include <eigen3/Eigen/Dense>
#include "include/poc_function.hpp"

#define DEBUG 0

using namespace Eigen;
using namespace std;

void POCFunction::_collision_angle_intervals(VectorXd* r, VectorXd* phi){

    vector<double> dists_e, dists_o;
    VectorXd r_d, theta_lb1, theta_ub1, theta_lb2, theta_ub2, theta_coll, theta_1, theta_2;
    double L_cir_e, L_cir_o, r_abriss;
    int n_it, cir_o_r, cir_o_f, counter = 0;
    _all_int = MatrixXd::Zero(_n_rows, _n_cols);
    _cir_nbr = MatrixXd::Zero(1, _n_cols);
    VectorXd x_o = (*r).array() * (*phi).array().cos();
    VectorXd y_o = (*r).array() * (*phi).array().sin();
    VectorXd phi_d(x_o.size());

    if (_paras.N_co % 2 == 0) {
        n_it = _paras.N_co/2;
    } else {
        n_it = (_paras.N_co + 1)/2;
    }
    _get_distances(_paras.d_ce, _paras.d_co, _paras.N_ce, _paras.N_co, &dists_e, &dists_o);
    #if DEBUG
        cout << "dists_e: " << dists_e << endl;
        cout << "dists_o: " << dists_o << endl;
        cout<< "dists_e.size() *  n_it/(4*_paras.N_ce): " << dists_e.size() *  n_it/(4*_paras.N_ce) << endl;
    #endif
    for (int i_e = 0; i_e < dists_e.size(); ++i_e) {
        L_cir_e = dists_e[_paras.N_ce - i_e - 1];
        r_d = ((x_o.array() - L_cir_e).pow(2) + y_o.array().pow(2)).sqrt();

        for (int i = 0; i < x_o.size(); ++i) {
            phi_d(i) = atan2(y_o(i), x_o(i) - L_cir_e);             // atan2(y_o, x_o - L_cir_e) Eigen does not support atan2
            phi_d(i) = fmod(phi_d(i) + 2 * M_PI, 2 * M_PI);         // Modulo 2 * M_PI, + 2 * M_PI to avoid negative numbers
        }

        cir_o_r = _paras.N_co;
        cir_o_f = 1;
        #if DEBUG
            cout << "L_cir_e: " << L_cir_e << endl;
        #endif
        for (int i_o = 0; i_o < n_it; ++i_o) {
            L_cir_o = dists_o[_paras.N_co - i_o - 1];
            #if DEBUG
                cout << "L_cir_o: " << L_cir_o << endl;
            #endif
            if (L_cir_o == 0) {
                theta_lb1 = VectorXd::Zero(_n_rows);
                theta_ub1 = VectorXd::Constant(_n_rows, 2 * M_PI);
                theta_ub1 = (r_d.array() > _paras.R).select(0, theta_ub1.array());

                theta_lb2 = theta_lb1;
                theta_ub2 = theta_ub1;
            }
            else {
                r_abriss = _paras.R + L_cir_o;
                theta_coll = ((L_cir_o * L_cir_o + r_d.array().square() - _paras.R * _paras.R) / (2 * L_cir_o * r_d.array())).acos();
                theta_coll = (r_d.array() >= r_abriss).select(0, theta_coll.array());
                
                theta_1 = phi_d.array() + M_PI - theta_coll.array();
                theta_2 = phi_d.array() + M_PI + theta_coll.array();

                if (_paras.R >= L_cir_o) {
                    theta_1 = (r_d.array() <= _paras.R - L_cir_o).select(0, theta_1.array());
                    theta_2 = (r_d.array() <= _paras.R - L_cir_o).select(2 * M_PI, theta_2.array());
                }
                else {
                    theta_1 = (r_d.array() <= L_cir_o - _paras.R).select(0, theta_1.array());
                    theta_2 = (r_d.array() <= L_cir_o - _paras.R).select(0, theta_2.array());
                }
                
                
                theta_lb1 = theta_1;
                theta_ub1 = theta_2;
                
                theta_lb1 = ((theta_ub1.array() - theta_lb1.array()) >= M_PI).select(0, theta_lb1.array());
                theta_ub1 = ((theta_ub1.array() - theta_lb1.array()) >= M_PI).select(M_PI, theta_ub1.array());
                
                theta_lb2 = theta_lb1.array() + M_PI;
                theta_ub2 = theta_ub1.array() + M_PI;
                
            }
            if (cir_o_f == cir_o_r){
                _all_int.block(0, counter, _n_rows, 2) << theta_lb1, theta_ub1;
                _cir_nbr.block(0, counter, 1, 2) << i_e + 1, cir_o_f;
                #if DEBUG
                    cout << "counter: " << counter << endl;
                    cout << "ego_cir: " << i_e + 1 << ", obj_cir_f: " << cir_o_f << endl;
                #endif
                counter += 2;
            } else{
                _all_int.block(0, counter, _n_rows, 4) << theta_lb1, theta_ub1, theta_lb2, theta_ub2;
                _cir_nbr.block(0, counter, 1, 4) << i_e + 1, cir_o_f, i_e + 1, cir_o_r;
                #if DEBUG
                    cout << "counter: " << counter << endl;
                    cout << "ego_cir: " << i_e + 1 << ", obj_cir_f: " << cir_o_f << ", obj_cir_r: " << cir_o_r << endl;
                #endif
                counter += 4;
            }
            #if DEBUG
                cout << "_all_int: \n" << _all_int << endl;
                cout << "_cir_nbr: " << _cir_nbr << endl;
            #endif
            cir_o_r--;
            cir_o_f++;
        }
        
    }

}