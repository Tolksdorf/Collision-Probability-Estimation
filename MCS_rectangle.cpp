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
#include <cmath>
#include <random>
#include <bits/stdc++.h>
#include <eigen3/Eigen/Dense>
#include "include/mcs_rectangle.hpp"

using namespace std;
using namespace Eigen;

double normrnd(double mean, double sigma) {
    random_device rd;
    mt19937 gen(rd());
    normal_distribution<double> dist(mean, sigma);
    return dist(gen);
}

Matrix2d rotation_matrix(double angle) {
    Matrix2d R;
    R << cos(angle), -sin(angle),
         sin(angle),  cos(angle);
    return R;
}

MatrixXd get_rectangle_corners(Vector3d states, double length, double width) {
    MatrixXd corners(4, 2);
    Matrix2d R = rotation_matrix(states(2));

    corners.row(2) = R * Vector2d(0.5 * length, 0.5 * width);
    corners.row(3) = R * Vector2d(0.5 * length, -0.5 * width);
    corners.row(1) = R * Vector2d(-0.5 * length, 0.5 * width);
    corners.row(0) = R * Vector2d(-0.5 * length, -0.5 * width);
    corners.rowwise() += states.head<2>().transpose(); 
    MatrixXd closed_corners(5, 2);
    closed_corners << corners, corners.row(0);
    return closed_corners;
}

bool check_side(const MatrixXd& rect1, const MatrixXd& rect2, Vector2d edgeNormal, Vector2d edgePoint) {
    VectorXd projections1 = (rect1.rowwise() - edgePoint.transpose()) * edgeNormal;
    VectorXd projections2 = (rect2.rowwise() - edgePoint.transpose()) * edgeNormal;

    return !(projections1.minCoeff() > projections2.maxCoeff() || projections1.maxCoeff() < projections2.minCoeff());
}

bool RectIntersect(const MatrixXd& rect1, const MatrixXd& rect2) {
    for (int j = 0; j < 4; ++j) {
        Vector2d edge = rect1.row((j + 1) % 4) - rect1.row(j);
        Vector2d edgeNormal(-edge.y(), edge.x());  //perpendicular vector

        if (!check_side(rect1, rect2, edgeNormal, rect1.row(j).transpose())) return false;
    }
    for (int j = 0; j < 4; ++j) {
        Vector2d edge = rect2.row((j + 1) % 4) - rect2.row(j);
        Vector2d edgeNormal(-edge.y(), edge.x());

        if (!check_side(rect1, rect2, edgeNormal, rect2.row(j).transpose())) return false;
    }
    return true;
}

bool rect_intersect(Vector3d states_o, double length_e, double width_e, double length_o, double width_o) {
    Vector3d states_e = {0, 0, 0};
    MatrixXd rect_e = get_rectangle_corners(states_e, length_e, width_e);
    MatrixXd rect_o = get_rectangle_corners(states_o, length_o, width_o);
    return RectIntersect(rect_e, rect_o);
}
double MCS_rectangle(Vector3d sigma_o, Vector3d y_o, Vector2d dim_e, Vector2d dim_o, int N) {
    int P_rect = 0;

    for (int j = 0; j < N; ++j) {
        Vector3d r_o(normrnd(y_o(0), sigma_o(0)), normrnd(y_o(1), sigma_o(1)), normrnd(y_o(2), sigma_o(2)));
        if (rect_intersect(r_o, dim_e(0), dim_e(1), dim_o(0), dim_o(1))) {
            P_rect++;
        }
    }
    return static_cast<double>(P_rect) / N;
}