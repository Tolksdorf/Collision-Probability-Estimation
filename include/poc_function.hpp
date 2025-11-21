/*
Written by Leon Tolksdorf, 2025. Please refer to: 

@article{CollisionProbabilityEstimation,
  title   = {Collision Probability Estimation for Optimization-based Vehicular Motion Planning},
  author  = {Tolksdorf, Leon and Tejada, Arturo and Birkner, Christian and van de Wouw, Nathan},
  journal = {arXiv preprint arXiv:2505.21161} ,
  year    = {2025}
}
*/

#pragma once
#include <cmath>
#include <iostream>
#include <vector>
#include <eigen3/Eigen/Dense>
#include <casadi/casadi.hpp>

class POCFunction 
{
    public:
        //constructor
        POCFunction(); 

        //structs for handling parameters and dynamic data
        struct parameters {
            double d_ce;
            double d_co;
            double R;
            int N_co;
            int N_ce;
            int N_beta;
            int Nbr_samples_dim;
            int integral_method; 
        };

        //define functions
        void poc_parameter_values(parameters* values);

        //global variables
        casadi::DMDict poc_fu_args;
        casadi::Function f_poc;

    private:

        //define functions        
        void _collision_angle_intervals(Eigen::VectorXd* r, Eigen::VectorXd* phi);
        void _get_distances(double L_e, double L_o, int n_cir_e, int n_cir_o, std::vector<double>* dists_e, std::vector<double>* dists_o);
        void _sort_intervals();
        void _int_unions(const Eigen::MatrixXd* Int1, const Eigen::MatrixXd* Int2, Eigen::MatrixXd *Int1_new, Eigen::MatrixXd *Int2_new);

        Eigen::VectorXd _modified_mod(Eigen::VectorXd bound);
        Eigen::VectorXd _mod_2pi(Eigen::VectorXd vec);

        casadi::SX _eigenToCasadiSX(const Eigen::MatrixXd& eigen_mat);
        casadi::SX _trapezoidal_2Dintegral(casadi::SX vector);
        casadi::SX _Simpson_2Dintegral(casadi::SX vector);
        casadi::Function _twoD_integral();

        //create struct instance
        parameters _paras;

        //private variables
        int _perspective;
        int _n_rows;
        int _n_cols;

        double _h_r;
        double _h_phi;

        Eigen::MatrixXd _thetas;
        Eigen::MatrixXd _all_int;
        Eigen::MatrixXd _cir_nbr;

        Eigen::VectorXd _r_range;
        Eigen::VectorXd _phi_range;

        casadi::SX _mu_o       = casadi::SX::sym("mu_o", 3);
        casadi::SX _sigma_conf = casadi::SX::sym("sigma_conf", 3);

        Eigen::VectorXd _lb_int1;
        Eigen::VectorXd _ub_int1;
        Eigen::VectorXd _lb_int2;
        Eigen::VectorXd _ub_int2;
        Eigen::VectorXd _m_lb_int1;
        Eigen::VectorXd _m_ub_int1;
        Eigen::VectorXd _m_lb_int2;
        Eigen::VectorXd _m_ub_int2;
        Eigen::VectorXd _lb_int1_new;
        Eigen::VectorXd _ub_int1_new;
        Eigen::VectorXd _lb_int2_new;
        Eigen::VectorXd _ub_int2_new;
        Eigen::VectorXd _logic_I1wrap;
        Eigen::VectorXd _logic_I2wrap;
        Eigen::VectorXd _logic_I1contI2;
        Eigen::VectorXd _logic_I2contI1;
        Eigen::VectorXd _logic_I1overlapI2;
        Eigen::VectorXd _logic_I2overlapI1;
        Eigen::VectorXd _logic_I1doubleoverlapI2;
        Eigen::VectorXd _logic_I2doubleoverlapI1;
        Eigen::VectorXd _logic_I1disjointI2;
};