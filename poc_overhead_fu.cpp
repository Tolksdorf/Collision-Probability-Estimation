#include <iostream>
#include <eigen3/Eigen/Dense>
#include <casadi/casadi.hpp>
#include "include/poc_function.hpp"

using namespace std;
using namespace Eigen;
using namespace casadi;

POCFunction::POCFunction(){
    //constructor
    poc_parameter_values(&paras);
    _n_rows = paras.Nbr_samples_dim * paras.Nbr_samples_dim;

    //determine twice the number of collision angle intervals
    if (paras.n_cir_o % 2 == 0) {_n_cols = paras.n_cir_o/2 * 4 * paras.n_cir_e;
    } else {_n_cols = (paras.n_cir_o + 1)/2 * 4 * paras.n_cir_e;}

    //initialize the integration ranges
    double r_lb = 0;
    double r_ub = paras.R + paras.L_o/2 + (paras.n_cir_o/2 - 1) * paras.L_o + (paras.n_cir_e - 1)/2 * paras.L_e; 
    double phi_lb = 0;
    double phi_ub = 2 * M_PI;
    VectorXd segment = VectorXd::Ones(paras.Nbr_samples_dim);
    VectorXd r_single_range   = VectorXd::LinSpaced(paras.Nbr_samples_dim, r_lb, r_ub);
    VectorXd phi_single_range = VectorXd::LinSpaced(paras.Nbr_samples_dim, phi_lb, phi_ub);
    _r_range   = VectorXd::Zero(_n_rows);
    _phi_range = VectorXd::Zero(_n_rows);
    _h_r   = r_single_range(1) - r_single_range(0);
    _h_phi = phi_single_range(1) - phi_single_range(0);

    for (int k = 0; k < paras.Nbr_samples_dim; k++){
        _r_range.segment(k*paras.Nbr_samples_dim, paras.Nbr_samples_dim) << segment * r_single_range(k);
        _phi_range.segment(k*paras.Nbr_samples_dim, paras.Nbr_samples_dim) <<  phi_single_range;
    }

    //get collision angle intervals
    _collision_angle_intervals(&_r_range, &_phi_range);

    //sort intervals
    _sort_intervals();	
    
    //numerical integration
    f_poc  = _twoD_integral();

}
