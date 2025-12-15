/*
Written by Leon Tolksdorf, 2025. Please refer to: 

@article{CollisionProbabilityEstimation,
  title   = {Collision Probability Estimation for Optimization-based Vehicular Motion Planning},
  author  = {Tolksdorf, Leon and Tejada, Arturo and Birkner, Christian and van de Wouw, Nathan},
  journal = {arXiv preprint arXiv:2505.21161} ,
  year    = {2025}
}
*/

#include <casadi/casadi.hpp>
#include <cmath>
#include <vector>
#include <iostream>
#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/SpecialFunctions>
#include "include/poc_function.hpp"

using namespace std;
using namespace Eigen;
using namespace casadi;


casadi::Function POCFunction::_twoD_integral(){
    
    SX poc_value, term1, term2;
    SX p_vector = SX::zeros(_n_rows, 1);
    SX sigma_x = _sigma_conf(0);
    SX sigma_y = _sigma_conf(1);
    SX sigma_t = _sigma_conf(2);
    SX sqrt2_sigma_t = sqrt(2) * sigma_t;
    int n_probs = _n_cols / 2;
    double two_pi = 2.0 * M_PI;
    SX d_mu = sqrt(pow(_mu_o(0), 2) + pow(_mu_o(1), 2));
    SX phi_mu = fmod(atan2(_mu_o(1), _mu_o(0)), two_pi);
    SX theta_mu = fmod(_mu_o(2), two_pi);
    SX p_theta = SX::zeros(_n_rows, n_probs);
    int count = 0;
    SX thetas = _eigenToCasadiSX(_thetas);

    for (int n_pr = 0; n_pr < n_probs; n_pr++) {
        for (int k = -static_cast<int>(_paras.N_beta); k <= static_cast<int>(_paras.N_beta); k++) {
            term1 = (thetas(Slice(), count + 1) - theta_mu + k * two_pi) / sqrt2_sigma_t;
            term2 = (thetas(Slice(), count    ) - theta_mu + k * two_pi) / sqrt2_sigma_t;
            p_theta(Slice(), n_pr) += erf(term1) - erf(term2);
        }
        count += 2;
    }

    VectorXd cos_phi = _phi_range.array().cos();
    VectorXd sin_phi = _phi_range.array().sin();
    
    SX cos_phi_cas = _eigenToCasadiSX(cos_phi);
    SX sin_phi_cas = _eigenToCasadiSX(sin_phi);
    SX r_cas       = _eigenToCasadiSX(_r_range);
    SX cos_phi_mu = cos(phi_mu);
    SX sin_phi_mu = sin(phi_mu);

    SX p_c1_c2 = r_cas * exp(-0.5 * ((pow(r_cas * cos_phi_cas - d_mu * cos_phi_mu, 2) / pow(sigma_x, 2)) +
                                     (pow(r_cas * sin_phi_cas - d_mu * sin_phi_mu, 2) / pow(sigma_y, 2))));
    
    p_theta = sum2(p_theta);
    p_vector = p_c1_c2 * p_theta;

    if (_paras.integral_method == 1){
        poc_value = _trapezoidal_2Dintegral(p_vector);
    } else if (_paras.integral_method == 2) {
        poc_value = _Simpson_2Dintegral(p_vector);
    } else{throw invalid_argument("incorrect integration method selected. Set either integral_method = 1 for trapezoidal or integration_method = 2 for Simpson");}

    poc_value *= 1/(4 * M_PI * sigma_x * sigma_y);
    poc_value = cse(poc_value); //simplyfies symbolic exspressions
    return Function("poc_function", {_mu_o, _sigma_conf}, {poc_value}, 
                                    {"_mu_o", "_sigma_conf"}, {"poc_value"});

}