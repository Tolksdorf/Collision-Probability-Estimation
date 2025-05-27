#include <iostream>
#include <vector>
#include <eigen3/Eigen/Dense>
#include "include/poc_function.hpp"

using namespace Eigen;
using namespace std;

void POCFunction::_collision_angle_intervals(VectorXd* r, VectorXd* phi){

    vector<double> dists_e, dists_o;
    _get_distances(paras.L_e, paras.L_o, paras.n_cir_e, paras.n_cir_o, &dists_e, &dists_o);

    _all_int = MatrixXd::Zero(_n_rows, _n_cols);
    _cir_nbr = MatrixXd::Zero(1, _n_cols);
    int counter = 1;
    VectorXd x_o = (*r).array() * (*phi).array().cos();
    VectorXd y_o = (*r).array() * (*phi).array().sin();
    
    for (int i_e = 0; i_e < dists_e.size(); ++i_e) {
        double L_cir_e = dists_e[i_e];
        VectorXd r_d = ((x_o.array() - L_cir_e).pow(2) + y_o.array().pow(2)).sqrt();
        VectorXd phi_d(x_o.size());
        for (int i = 0; i < x_o.size(); ++i) {
            phi_d(i) = atan2(y_o(i), x_o(i) - L_cir_e);             // atan2(y_o, x_o - L_cir_e) Eigen does not support atan2
            phi_d(i) = fmod(phi_d(i) + 2 * M_PI, 2 * M_PI);         // Modulo 2 * M_PI, + 2 * M_PI to avoid negative numbers
        }

        int cir_o_r = paras.n_cir_o;
        int cir_o_f = 1;
        
        
        for (int i_o = 0; i_o < _n_cols/(4*paras.n_cir_e); ++i_o) {
            double L_cir_o = dists_o[paras.n_cir_o - i_o - 1];

            VectorXd theta_lb1, theta_ub1, theta_lb2, theta_ub2;

            if (L_cir_o == 0) {
                theta_lb1 = VectorXd::Zero(_n_rows);
                theta_ub1 = VectorXd::Constant(_n_rows, 2 * M_PI);
                theta_ub1 = (r_d.array() > paras.R).select(0, theta_ub1.array());

                theta_lb2 = theta_lb1;
                theta_ub2 = theta_ub1;
            }
            else {
                double r_abriss = paras.R + L_cir_o;
                VectorXd theta_coll = ((L_cir_o * L_cir_o + r_d.array().square() - paras.R * paras.R) / (2 * L_cir_o * r_d.array())).acos();
                theta_coll = (r_d.array() >= r_abriss).select(0, theta_coll.array());
                
                VectorXd theta_1 = phi_d.array() + M_PI - theta_coll.array();
                VectorXd theta_2 = phi_d.array() + M_PI + theta_coll.array();

                if (paras.R >= L_cir_o) {
                    theta_1 = (r_d.array() <= paras.R - L_cir_o).select(0, theta_1.array());
                    theta_2 = (r_d.array() <= paras.R - L_cir_o).select(2 * M_PI, theta_2.array());
                }
                else {
                    theta_1 = (r_d.array() <= L_cir_o - paras.R).select(0, theta_1.array());
                    theta_2 = (r_d.array() <= L_cir_o - paras.R).select(0, theta_2.array());
                }
                
                
                theta_lb1 = theta_1;
                theta_ub1 = theta_2;
                
                theta_lb1 = ((theta_ub1.array() - theta_lb1.array()) >= M_PI).select(0, theta_lb1.array());
                theta_ub1 = ((theta_ub1.array() - theta_lb1.array()) >= M_PI).select(M_PI, theta_ub1.array());
                
                theta_lb2 = theta_lb1.array() + M_PI;
                theta_ub2 = theta_ub1.array() + M_PI;
                
            }
            _all_int.block(0, 4 * (counter - 1), _n_rows, 4) << theta_lb1, theta_ub1, theta_lb2, theta_ub2;
            _cir_nbr.block(0, 4 * (counter - 1), 1, 4) << i_e + 1, cir_o_f, i_e + 1, cir_o_r;
            cir_o_r--;
            cir_o_f++;
            counter++;
            
        }
        
    }

}