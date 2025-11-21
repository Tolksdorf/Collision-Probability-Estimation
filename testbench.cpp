/*
Written by Leon Tolksdorf, 2025. Please refer to: 

@article{CollisionProbabilityEstimation,
  title   = {Collision Probability Estimation for Optimization-based Vehicular Motion Planning},
  author  = {Tolksdorf, Leon and Tejada, Arturo and Birkner, Christian and van de Wouw, Nathan},
  journal = {arXiv preprint arXiv:2505.21161} ,
  year    = {2025}
}
*/

#include <eigen3/Eigen/Dense>
#include <iostream>
#include <random>
#include <chrono>
#include <vector>
#include <casadi/casadi.hpp>
#include <random>
#include <fstream>
#include "include/poc_function.hpp"
#include "include/mcs_rectangle.hpp"

using namespace Eigen;
using namespace std;
using namespace casadi;


void POCFunction::poc_parameter_values(parameters* values) {
        double length_e = 4.5;        //length ego rectangle [m]
        double width_e = 2;           //width ego rectangle [m]
        double length_o = 4.5;        //length object rectangle [m]
        double width_o = 2;           //width object rectangle [m]
        values->N_co = 3;             //number of circles approximating the object rectangle
        values->N_ce = 3;             //number of circles approximating the ego rectangle
        values->N_beta = 3;           //number of windings approximating the wrapped Gaussian distribution
        values->Nbr_samples_dim = 20; //number of samples per dimensions for the numerical integration in \rho and \phi
        values->integral_method = 1;  //1: trapezoidal method, 2: Simpson method, everything else throws an error
        double radius_e = sqrt(pow(length_e / 2, 2) / pow(values->N_ce, 2) + pow(width_e, 2)/ 4);
        double radius_o = sqrt(pow(length_o / 2, 2) / pow(values->N_co, 2) + pow(width_o, 2)/ 4);
        values->R = radius_e + radius_o;
        values->d_ce = 2 * sqrt(pow(radius_e, 2) - pow(width_e, 2) / 4);
        values->d_co = 2 * sqrt(pow(radius_o, 2) - pow(width_o, 2) / 4);;
    }

double random_double(double min_val, double max_val) {
        static random_device rd;
        static mt19937 gen(rd());  //Mersenne twister
        std::uniform_real_distribution<double> dist(min_val, max_val);
        return dist(gen);
    }

void save_vectors_CSV(const vector<double>& col1, const std::string& filename) {

    std::ofstream file(filename);  //open the file for writing
    if (!file) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }
    //column headers
    file << "two-circle" << "\n";
    size_t numRows = col1.size();
    //write data row-wise
    for (size_t i = 0; i < numRows; ++i) {
        file << col1[i] << "\n";
    }
    file.close(); 
    std::cout << "Data saved to " << filename << std::endl;
}

int main (void){

    Vector3d mu_o = {0, -2, 0.25*M_PI}; //the object's configuration is measured with respect to the ego vehicle, i.e., in the ego's coordiante system. \theta_o is measured with respect to the ego's x-axis.
    Vector3d sigma_o = {1, 1, 1};       //objects standard deviations in x, y, \theta. NOTE: if set too low, the Nbr_samples_dim must be increased, otherwise the intergration grid becomes too coarse.
    Vector2d dim_o = {4.5, 2};          //[length, width] of the object in [m].
    Vector2d dim_e = {4.5, 2};          //[length, width] of the ego in [m].
    int N_MCS_samples = 1000;           //number of Monte Carlo samples for comparative Monte Carlo method on rectangles.
    int N_simulations = 1000;           //number of simulations to calculate the average computional time over.

    vector<double> cir_calc_t;
    vector<double> init_calc_t;
    double yo_x, yo_y, yo_theta, sigma_x, sigma_y, sigma_theta;
    MX POC;
    auto start = chrono::high_resolution_clock::now();
    POCFunction poc_function;
    auto end = std::chrono::high_resolution_clock::now();
    chrono::duration<double, std::milli> duration = end - start;
    init_calc_t.push_back(duration.count());
    cout << "Initalizing POC function, calculation time: " << duration.count() << " ms." << endl;
    //calculate the average POC estimation time with random inputs
    for (int i = 0; i < N_simulations; i++){
        yo_x = random_double(0, 20);
        yo_y = random_double(0, 20);
        yo_theta = random_double(0, 6.28);
        sigma_x = random_double(0.1, 10);
        sigma_y = random_double(0.1, 10);
        sigma_theta = random_double(0.1, 5);

        start = std::chrono::high_resolution_clock::now();
        poc_function.poc_fu_args["_mu_o"]       = {yo_x, yo_y, yo_theta};
        poc_function.poc_fu_args["_sigma_conf"] = {sigma_x, sigma_y, sigma_theta};
        POC = poc_function.f_poc(poc_function.poc_fu_args).at("poc_value"); 
        end = chrono::high_resolution_clock::now();
        duration = end - start;
        cir_calc_t.push_back(duration.count());
    }
    save_vectors_CSV(cir_calc_t, "results.csv");
    double total_time = accumulate(cir_calc_t.begin(), cir_calc_t.end(), 0.0);

    //calculate the POC value for the given inputs
    poc_function.poc_fu_args["_mu_o"]       = {mu_o[0], mu_o[1], mu_o[2]};
    poc_function.poc_fu_args["_sigma_conf"] = {sigma_o[0], sigma_o[1], sigma_o[2]};
    POC = poc_function.f_poc(poc_function.poc_fu_args).at("poc_value"); 
    cout << "Symbolic POC calculation: " << POC << ", avg. calculation time: " << total_time/N_simulations << " ms." << endl;

    start = chrono::high_resolution_clock::now();
    double poc_MCS = MCS_rectangle(sigma_o, mu_o, dim_e, dim_o, N_MCS_samples);
    end = std::chrono::high_resolution_clock::now();
    duration = end - start;
    cout << "MCS POC calculation: " << poc_MCS << ", calculation time: " << duration.count() << " ms." << endl;
	return 0;
	
}