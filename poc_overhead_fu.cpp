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
#include <eigen3/Eigen/Dense>
#include <casadi/casadi.hpp>
#include "include/poc_function.hpp"

using namespace std;
using namespace Eigen;
using namespace casadi;

POCFunction::POCFunction(){
    //constructor
    poc_parameter_values(&_paras);
    _n_rows = _paras.Nbr_samples_dim * _paras.Nbr_samples_dim;

    //determine twice the number of collision angle intervals
    _n_cols = 2 * _paras.N_ce * _paras.N_co;

    //initialize the integration ranges
    double r_lb = 0;
    double r_ub = _paras.R + _paras.d_ce/2 * (_paras.N_co - 1) + _paras.d_ce/2 * (_paras.N_ce - 1); 
    double phi_lb = 0;
    double phi_ub = 2 * M_PI;
    VectorXd segment = VectorXd::Ones(_paras.Nbr_samples_dim);
    VectorXd r_single_range   = VectorXd::LinSpaced(_paras.Nbr_samples_dim, r_lb, r_ub);
    VectorXd phi_single_range = VectorXd::LinSpaced(_paras.Nbr_samples_dim, phi_lb, phi_ub);
    _r_range   = VectorXd::Zero(_n_rows);
    _phi_range = VectorXd::Zero(_n_rows);
    _h_r   = r_single_range(1) - r_single_range(0);
    _h_phi = phi_single_range(1) - phi_single_range(0);

    for (int k = 0; k < _paras.Nbr_samples_dim; k++){
        _r_range.segment(k*_paras.Nbr_samples_dim, _paras.Nbr_samples_dim) << segment * r_single_range(k);
        _phi_range.segment(k*_paras.Nbr_samples_dim, _paras.Nbr_samples_dim) <<  phi_single_range;
    }

    //get collision angle intervals
    _collision_angle_intervals(&_r_range, &_phi_range);

    //sort intervals
    _sort_intervals();	
    
    //numerical integration
    f_poc  = _twoD_integral();

}
