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
        double length_e = 4.5;
        double width_e = 2;
        double length_o = 4.5;
        double width_o = 2;
        values->n_cir_o = 2;
        values->n_cir_e = 2;
        values->n_w = 3;
        values->Nbr_samples_dim = 20;
        values->integral_method = 1; 
        double radius_e = sqrt(pow(length_e / 2, 2) / pow(values->n_cir_e, 2) + pow(width_e, 2)/ 4);
        double radius_o = sqrt(pow(length_o / 2, 2) / pow(values->n_cir_o, 2) + pow(width_o, 2)/ 4);
        values->R = radius_e + radius_o;
        values->L_e = 2 * sqrt(pow(radius_e, 2) - pow(width_e, 2) / 4);
        values->L_o = 2 * sqrt(pow(radius_o, 2) - pow(width_o, 2) / 4);;
    }

double random_double(double min_val, double max_val) {
        static random_device rd;
        static mt19937 gen(rd());  //Mersenne Twister RNG
        std::uniform_real_distribution<double> dist(min_val, max_val);
        return dist(gen);
    }

void saveVectorsToCSV(const vector<double>& col1, const std::string& filename) {

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

    file.close();  // Close the file
    std::cout << "Data saved to " << filename << std::endl;
}

int main (void){

    Vector3d y_e = {0, 0, 0};
    Vector3d mu_o = {0, 7, 0.785};
    Vector3d sigma_o = {3, 3, 2.121};
    Vector2d dim_o = {4.5, 2};
    Vector2d dim_e = {4.5, 2};
    int N_MCS_samples = 1000;
    int N_simulations = 10000;
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

    for (int i = 0; i < N_simulations; i++){
        yo_x = random_double(0, 20);
        yo_y = random_double(0, 20);
        yo_theta = random_double(0, 6.28);
        sigma_x = random_double(0.1, 10);
        sigma_y = random_double(0.1, 10);
        sigma_theta = random_double(0.1, 5);

        start = std::chrono::high_resolution_clock::now();
        poc_function.poc_fu_args["_y_e"]        = {y_e[0], y_e[1], y_e[2]};
        poc_function.poc_fu_args["_mu_o"]       = {yo_x, yo_y, yo_theta};
        poc_function.poc_fu_args["_sigma_conf"] = {sigma_x, sigma_y, sigma_theta};
        POC = poc_function.f_poc(poc_function.poc_fu_args).at("poc_value"); 
        end = chrono::high_resolution_clock::now();
        duration = end - start;
        cir_calc_t.push_back(duration.count());
    }
    saveVectorsToCSV(cir_calc_t, "results.csv");
    double total_time = accumulate(cir_calc_t.begin(), cir_calc_t.end(), 0.0);
    cout << "Symbolic POC calculation: " << POC << ", avg. calculation time: " << total_time/N_simulations << " ms." << endl;

    start = chrono::high_resolution_clock::now();
    double poc_MCS = MCS_rectangle(sigma_o, y_e, mu_o, dim_e, dim_o, N_MCS_samples);
    end = std::chrono::high_resolution_clock::now();
    duration = end - start;
    cout << "MCS POC calculation: " << poc_MCS << ", calculation time: " << duration.count() << " ms." << endl;
	return 0;
	
}