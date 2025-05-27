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

SX POCFunction::_eigenToCasadiSX(const MatrixXd& eigen_mat) {
    //note that this function is mainly used to convert Eigen::VectorXd to SX. As this type is also allowed 
    int rows = eigen_mat.rows();
    int cols = eigen_mat.cols();
    vector<double> data(eigen_mat.data(), eigen_mat.data() + eigen_mat.size());
    return SX::reshape(SX(data), rows, cols);
}

VectorXd POCFunction::_mod_2pi(VectorXd vec){
    return vec.array().unaryExpr([](double x) {return fmod(x, 2.0 * M_PI);});
}

VectorXd POCFunction::_modified_mod(VectorXd bound){
    //usually mod(N*pi, 2*pi) = 0, this function returns 2*pi instead otherwise collision angle intervals of [0, 2*pi] would collaps.
    bound = (_mod_2pi(bound).array() == 0 && bound.array() > 2*M_PI).select(2*M_PI, bound); //for all positive multiples of pi
    return (bound.array() > 2*M_PI).select(_mod_2pi(bound).array(), bound);
}

SX POCFunction::_trapezoidal_2Dintegral(SX vector){
    SX val = 0;
    SX integral_value = 0;
    int r_phi_count = 0;

    for(int i_phi = 0; i_phi < paras.Nbr_samples_dim; ++i_phi)
    {
        for(int j_r = 0; j_r < paras.Nbr_samples_dim; ++j_r)
        {
            val = vector(r_phi_count);
            // note: for a 2D gird, all four corners are weighted by 1/4 (see first for if statements). The sides of the grid, excluding the corners, 
            // are weighted by 1/2 (see next four if statements). All remaining points are weighted by 1 (no if statements required).
            if(i_phi == 0 || i_phi == (paras.Nbr_samples_dim - 1)){
                val *= 0.5;
            }
            if(j_r == 0 || j_r == (paras.Nbr_samples_dim - 1)){
                val *= 0.5;
            }
            integral_value += val;
            r_phi_count++;
        }
    }
    return integral_value * _h_r * _h_phi;
}


SX POCFunction::_Simpson_2Dintegral(SX vector){
    SX val = SX::zeros(paras.Nbr_samples_dim,1);
    SX integral_value = 0;
    int r_phi_count = 0;
    for (int i_phi = 0; i_phi < paras.Nbr_samples_dim; i_phi++){

        for(int j_r = 0; j_r < paras.Nbr_samples_dim; j_r++){

            if (j_r == 0 || j_r == paras.Nbr_samples_dim - 1){
                val(i_phi) +=  vector(r_phi_count);
            } else if (j_r % 2 == 0){
                val(i_phi) += 2 * vector(r_phi_count);
            } else{
                val(i_phi) += 4 * vector(r_phi_count);
            }
            r_phi_count++;
        }
        val(i_phi) = val(i_phi) * _h_r / 3.0;
        if(i_phi == 0 || i_phi == paras.Nbr_samples_dim - 1){
            integral_value += val(i_phi);

        } else if (i_phi % 2 == 0){
            integral_value += 2 * val(i_phi);

        } else {
            integral_value += 4 * val(i_phi);
        } 
    }
    return integral_value * _h_phi / 3.0;
}

void POCFunction::_get_distances(double L_e, double L_o, int n_cir_e, int n_cir_o, std::vector<double>* dists_e, std::vector<double>* dists_o) {
	
    (*dists_o).resize(n_cir_o, 0.0);  //object distances of each circle of the center in local coordinates
    (*dists_e).resize(n_cir_e, 0.0);  //ego distances of each circle of the center in local coordinates

    //calculate ego distances
    if (n_cir_e % 2 == 0) {  //if even
        for (int i_ee = 0; i_ee < n_cir_e / 2; ++i_ee) {
            (*dists_e)[i_ee] = -(L_e / 2.0 + (n_cir_e / 2.0 - i_ee - 1.0) * L_e);
            (*dists_e)[n_cir_e - 1 - i_ee] = L_e / 2.0 + (n_cir_e / 2.0 - i_ee - 1.0) * L_e;
        }
    }
    else {  //if odd
        for (int i_ee = 0; i_ee < (n_cir_e - 1) / 2; ++i_ee) {
            (*dists_e)[i_ee] = -(L_e + ((n_cir_e - 1.0) / 2.0 - i_ee - 1.0) * L_e);
            (*dists_e)[n_cir_e - 1 - i_ee] = L_e + ((n_cir_e - 1.0) / 2.0 - i_ee - 1.0) * L_e;
        }
    }

    //calculate object distances
    if (n_cir_o % 2 == 0) {  //if even
        for (int i_oo = 0; i_oo < n_cir_o / 2; ++i_oo) {
            (*dists_o)[i_oo] = -(L_o / 2.0 + (n_cir_o / 2.0 - i_oo - 1.0) * L_o);
            (*dists_o)[n_cir_o - 1 - i_oo] = L_o / 2.0 + (n_cir_o / 2.0 - i_oo - 1.0) * L_o;
        }
    }
    else {  //if odd
        for (int i_oo = 0; i_oo < (n_cir_o - 1) / 2; ++i_oo) {
            (*dists_o)[i_oo] = -(L_o + ((n_cir_o - 1.0) / 2.0 - i_oo - 1.0) * L_o);
            (*dists_o)[n_cir_o - 1 - i_oo] = L_o + ((n_cir_o - 1.0) / 2.0 - i_oo - 1.0) * L_o;
        }
    }
}
