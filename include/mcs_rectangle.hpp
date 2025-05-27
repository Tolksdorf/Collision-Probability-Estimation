#include <iostream>
#include <vector>
#include <eigen3/Eigen/Dense>

using namespace std;
using namespace Eigen;

double normrnd(double mean, double sigma);
Matrix2d rotation_matrix(double angle);
MatrixXd get_rectangle_corners(Vector3d states, double length, double width);
bool CheckSide(const MatrixXd& rect1, const MatrixXd& rect2, Vector2d edgeNormal, Vector2d edgePoint);
bool RectIntersect(const MatrixXd& rect1, const MatrixXd& rect2);
bool rect_intersect(Vector3d states_e, Vector3d states_o, double length_e, double width_e, double length_o, double width_o);
double MCS_rectangle(Vector3d sigma_o, Vector3d y_e, Vector3d y_o, Vector2d dim_e, Vector2d dim_o, int N);