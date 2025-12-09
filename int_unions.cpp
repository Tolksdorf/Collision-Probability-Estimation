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
#include <utility>
#include <iostream>
#include "include/poc_function.hpp"

using namespace Eigen;
using namespace casadi;
using namespace std;

#define DEBUG 0

void POCFunction::_int_unions(const MatrixXd* Int1, const MatrixXd* Int2, MatrixXd *Int1_new, MatrixXd *Int2_new) {
    
    //extract lower and upper bounds of int1 and int2
    _lb_int1 = (*Int1).col(0);
    _ub_int1 = (*Int1).col(1);
    _lb_int2 = (*Int2).col(0);
    _ub_int2 = (*Int2).col(1);
    _m_lb_int1 = _modified_mod(_lb_int1);
    _m_ub_int1 = _modified_mod(_ub_int1);
    _m_lb_int2 = _modified_mod(_lb_int2);
    _m_ub_int2 = _modified_mod(_ub_int2);

    //initialize new variables for int1 and int2
    _lb_int1_new = _lb_int1;
    _ub_int1_new = _ub_int1;
    _lb_int2_new = _lb_int2;
    _ub_int2_new = _ub_int2;

    //check if interval wraps around
    _logic_I1wrap = VectorXd::Zero(_n_rows);
    _logic_I2wrap = VectorXd::Zero(_n_rows);
    _logic_I1wrap = (_m_lb_int1.array() > _m_ub_int1.array()).select(1, _logic_I1wrap);
    _logic_I2wrap = (_m_lb_int2.array() > _m_ub_int2.array()).select(1, _logic_I2wrap);

    //check wheather one interval contains the other, if both intervals are equal they both contain each other
    _logic_I1contI2 = VectorXd::Zero(_n_rows);
    _logic_I2contI1 = VectorXd::Zero(_n_rows);
    //-->check if interval size is >= 2.0:
    _logic_I1contI2 = (abs(_ub_int1.array() - _lb_int1.array()) >= 2.0 * M_PI).select(1, _logic_I1contI2);
    _logic_I2contI1 = (abs(_ub_int2.array() - _lb_int2.array()) >= 2.0 * M_PI).select(1, _logic_I2contI1);
    //-->check non-wrapping case:
    _logic_I1contI2 = (_m_lb_int1.array() <= _m_lb_int2.array() && _m_ub_int1.array() >= _m_ub_int2.array() && _logic_I1wrap.array() == 0 && _logic_I2wrap.array() == 0 && _logic_I1contI2.array() == 0).select(1, _logic_I1contI2);
    _logic_I2contI1 = (_m_lb_int2.array() <= _m_lb_int1.array() && _m_ub_int2.array() >= _m_ub_int1.array() && _logic_I1wrap.array() == 0 && _logic_I2wrap.array() == 0 && _logic_I2contI1.array() == 0).select(1, _logic_I2contI1);
    //-->check wrapping case, int1 wrapping:
    _logic_I1contI2 = (_m_lb_int1.array() <= _m_ub_int2.array() && _m_lb_int1.array() <= _m_lb_int2.array() && _logic_I1wrap.array() == 1 && _logic_I2wrap.array() == 0 && _logic_I1contI2.array() == 0).select(1, _logic_I1contI2);
    _logic_I1contI2 = (_m_ub_int1.array() >= _m_ub_int2.array() && _m_lb_int1.array() >= _m_lb_int2.array() && _logic_I1wrap.array() == 1 && _logic_I2wrap.array() == 0 && _logic_I1contI2.array() == 0).select(1, _logic_I1contI2);
    //-->check wrapping case, int2 wrapping:
    _logic_I2contI1 = (_m_lb_int1.array() >= _m_lb_int2.array() && _m_ub_int1.array() >= _m_ub_int2.array() && _logic_I1wrap.array() == 0 && _logic_I2wrap.array() == 1 && _logic_I2contI1.array() == 0 ).select(1, _logic_I2contI1);
    _logic_I2contI1 = (_m_lb_int1.array() <= _m_lb_int2.array() && _m_ub_int1.array() <= _m_ub_int2.array() && _logic_I1wrap.array() == 0 && _logic_I2wrap.array() == 1 && _logic_I2contI1.array() == 0 ).select(1, _logic_I2contI1);
    //-->check wrapping case, int1 and int2 wrapping:
    _logic_I1contI2 = (_m_lb_int1.array() <= _m_lb_int2.array() && _m_ub_int1.array() >= _m_ub_int2.array() && _logic_I1wrap.array() == 1 && _logic_I2wrap.array() == 1 && _logic_I1contI2.array() == 0).select(1, _logic_I1contI2);
    _logic_I2contI1 = (_m_lb_int2.array() <= _m_lb_int1.array() && _m_ub_int2.array() >= _m_ub_int1.array() && _logic_I1wrap.array() == 1 && _logic_I2wrap.array() == 1 && _logic_I2contI1.array() == 0).select(1, _logic_I2contI1);
    
    //check wheather there is an overlap among both intervals
    _logic_I1overlapI2 = VectorXd::Zero(_n_rows);
    _logic_I2overlapI1 = VectorXd::Zero(_n_rows);
    _logic_I1doubleoverlapI2 = VectorXd::Zero(_n_rows);
    _logic_I2doubleoverlapI1 = VectorXd::Zero(_n_rows);

    //-->check non-wrapping case, without touching:
    _logic_I1overlapI2       = (_m_lb_int2.array() < _m_ub_int1.array() && _m_ub_int2.array() > _m_ub_int1.array() && _m_lb_int2.array() > _m_lb_int1.array() && _logic_I1wrap.array() == 0 && _logic_I2wrap.array() == 0 && _logic_I1contI2.array() == 0 && _logic_I2contI1.array() == 0).select(1, _logic_I1overlapI2);
    _logic_I2overlapI1       = (_m_lb_int2.array() < _m_lb_int1.array() && _m_ub_int2.array() > _m_lb_int1.array() && _m_ub_int1.array() > _m_ub_int2.array() && _logic_I1wrap.array() == 0 && _logic_I2wrap.array() == 0 && _logic_I1contI2.array() == 0 && _logic_I2contI1.array() == 0).select(1, _logic_I2overlapI1);
    //-->check wrapping case, int1 wrapping, without touching:
    _logic_I1overlapI2       = (_m_lb_int2.array() < _m_ub_int1.array() && _m_ub_int2.array() > _m_ub_int1.array() && _m_ub_int2.array() < _m_lb_int1.array() && _logic_I1wrap.array() == 1 && _logic_I2wrap.array() == 0 && _logic_I1contI2.array() == 0 && _logic_I2contI1.array() == 0 && _logic_I1overlapI2.array() == 0).select(1, _logic_I1overlapI2);
    _logic_I2overlapI1       = (_m_lb_int2.array() < _m_lb_int1.array() && _m_ub_int2.array() > _m_lb_int1.array() && _m_lb_int2.array() > _m_ub_int1.array() && _logic_I1wrap.array() == 1 && _logic_I2wrap.array() == 0 && _logic_I1contI2.array() == 0 && _logic_I2contI1.array() == 0 && _logic_I2overlapI1.array() == 0).select(1, _logic_I2overlapI1);
    _logic_I1doubleoverlapI2 = (_m_lb_int2.array() < _m_ub_int1.array() && _m_lb_int1.array() < _m_ub_int2.array() && _logic_I1wrap.array() == 1 && _logic_I2wrap.array() == 0 && _logic_I1contI2.array() == 0 && _logic_I2contI1.array() == 0 && _logic_I1overlapI2.array() == 0).select(1, _logic_I1doubleoverlapI2);
    //-->check wrapping case, int2 wrapping, without touching:
    _logic_I2overlapI1       = (_m_lb_int1.array() < _m_ub_int2.array() && _m_ub_int1.array() > _m_ub_int2.array() && _m_ub_int1.array() < _m_lb_int2.array() && _logic_I1wrap.array() == 0 && _logic_I2wrap.array() == 1 && _logic_I1contI2.array() == 0 && _logic_I2contI1.array() == 0 && _logic_I2overlapI1.array() == 0).select(1, _logic_I2overlapI1);
    _logic_I1overlapI2       = (_m_lb_int1.array() < _m_lb_int2.array() && _m_ub_int1.array() > _m_lb_int2.array() && _m_lb_int1.array() > _m_ub_int2.array() && _logic_I1wrap.array() == 0 && _logic_I2wrap.array() == 1 && _logic_I1contI2.array() == 0 && _logic_I2contI1.array() == 0 && _logic_I1overlapI2.array() == 0).select(1, _logic_I1overlapI2);
    _logic_I2doubleoverlapI1 = (_m_lb_int1.array() < _m_ub_int2.array() && _m_lb_int2.array() < _m_ub_int1.array() && _logic_I1wrap.array() == 0 && _logic_I2wrap.array() == 1 && _logic_I1contI2.array() == 0 && _logic_I2contI1.array() == 0 && _logic_I2overlapI1.array() == 0).select(1, _logic_I2doubleoverlapI1);
    //-->check wrapping case, int1 and int2 wrapping:
    _logic_I1overlapI2       = (_m_lb_int1.array() < _m_lb_int2.array() && _m_lb_int1.array() > _m_ub_int2.array() && _logic_I1wrap.array() == 1 && _logic_I2wrap.array() == 1 && _logic_I1contI2.array() == 0 && _logic_I2contI1.array() == 0 &&  _logic_I1overlapI2.array() == 0).select(1, _logic_I1overlapI2);
    _logic_I2overlapI1       = (_m_lb_int2.array() < _m_lb_int1.array() && _m_lb_int2.array() > _m_ub_int1.array() && _logic_I1wrap.array() == 1 && _logic_I2wrap.array() == 1 && _logic_I1contI2.array() == 0 && _logic_I2contI1.array() == 0 &&  _logic_I2overlapI1.array() == 0).select(1, _logic_I2overlapI1);
    _logic_I1doubleoverlapI2 = (_m_lb_int1.array() < _m_lb_int2.array() && _m_lb_int1.array() <= _m_ub_int2.array() && _logic_I1wrap.array() == 1 && _logic_I2wrap.array() == 1 && _logic_I1contI2.array() == 0 && _logic_I2contI1.array() == 0 && _logic_I1doubleoverlapI2.array() == 0).select(1, _logic_I1doubleoverlapI2);
    _logic_I2doubleoverlapI1 = (_m_lb_int2.array() < _m_lb_int1.array() && _m_lb_int2.array() <= _m_ub_int1.array() && _logic_I1wrap.array() == 1 && _logic_I2wrap.array() == 1 && _logic_I1contI2.array() == 0 && _logic_I2contI1.array() == 0 && _logic_I2doubleoverlapI1.array() == 0).select(1, _logic_I2doubleoverlapI1);

    //check if intervals are disjoint
    _logic_I1disjointI2 = VectorXd::Zero(_n_rows);
    _logic_I1disjointI2 = (_logic_I1contI2.array() == 0 && _logic_I2contI1.array() == 0 && _logic_I1overlapI2.array() == 0 && _logic_I2overlapI1.array() == 0 && _logic_I1doubleoverlapI2.array() == 0).select(1, _logic_I1disjointI2);

    //if int1 contains int2, keep int1 and set int2 to zero
    _lb_int1_new = (_logic_I1contI2.array() == 1).select(_m_lb_int1, _lb_int1_new);
    _ub_int1_new = (_logic_I1contI2.array() == 1).select(_m_ub_int1, _ub_int1_new);
    _lb_int2_new = (_logic_I1contI2.array() == 1).select(0, _lb_int2_new);
    _ub_int2_new = (_logic_I1contI2.array() == 1).select(0, _ub_int2_new);

    //if int2 contains int1, write int2 into int1 and set int2 to zero
    _lb_int1_new = (_logic_I2contI1.array() == 1).select(_m_lb_int2, _lb_int1_new);
    _ub_int1_new = (_logic_I2contI1.array() == 1).select(_m_ub_int2, _ub_int1_new);
    _lb_int2_new = (_logic_I2contI1.array() == 1).select(0, _lb_int2_new);
    _ub_int2_new = (_logic_I2contI1.array() == 1).select(0, _ub_int2_new); 

    //if disjoint, leave them as they are
    _lb_int1_new = (_logic_I1disjointI2.array() == 1).select(_m_lb_int1, _lb_int1_new);
    _ub_int1_new = (_logic_I1disjointI2.array() == 1).select(_m_ub_int1, _ub_int1_new);
    _lb_int2_new = (_logic_I1disjointI2.array() == 1).select(_m_lb_int2, _lb_int2_new);
    _ub_int2_new = (_logic_I1disjointI2.array() == 1).select(_m_ub_int2, _ub_int2_new);

    //if int1 and int2 overlap, combine into int1 and set int2 to zero
    _lb_int1_new = (_logic_I1overlapI2.array() == 1).select(_m_lb_int1, _lb_int1_new);
    _ub_int1_new = (_logic_I1overlapI2.array() == 1).select(_m_ub_int2, _ub_int1_new);
    _lb_int2_new = (_logic_I1overlapI2.array() == 1).select(0, _lb_int2_new);
    _ub_int2_new = (_logic_I1overlapI2.array() == 1).select(0, _ub_int2_new);

    //if int2 and int1 overlap, combine into int1 and set int2 to zero
    _lb_int1_new = (_logic_I2overlapI1.array() == 1).select(_m_lb_int2, _lb_int1_new);
    _ub_int1_new = (_logic_I2overlapI1.array() == 1).select(_m_ub_int1, _ub_int1_new);
    _lb_int2_new = (_logic_I2overlapI1.array() == 1).select(0, _lb_int2_new);
    _ub_int2_new = (_logic_I2overlapI1.array() == 1).select(0, _ub_int2_new);

    //if double overlap, combine into int1 and set int2 to zero
    _lb_int1_new = (_logic_I1doubleoverlapI2.array() == 1).select(_m_lb_int1, _lb_int1_new);
    _ub_int1_new = (_logic_I1doubleoverlapI2.array() == 1).select(_m_ub_int2.array() + 2.0*M_PI, _ub_int1_new);
    _lb_int2_new = (_logic_I1doubleoverlapI2.array() == 1).select(0, _lb_int2_new);
    _ub_int2_new = (_logic_I1doubleoverlapI2.array() == 1).select(0, _ub_int2_new);

    //if double overlap, combine into int1 and set int2 to zero
    _lb_int1_new = (_logic_I2doubleoverlapI1.array() == 1).select(_m_lb_int2, _lb_int1_new);
    _ub_int1_new = (_logic_I2doubleoverlapI1.array() == 1).select(_m_ub_int1.array() + 2.0*M_PI, _ub_int1_new);
    _lb_int2_new = (_logic_I2doubleoverlapI1.array() == 1).select(0, _lb_int2_new);
    _ub_int2_new = (_logic_I2doubleoverlapI1.array() == 1).select(0, _ub_int2_new);

    //if an interval is negative in measure, add 2pi to the value of the upper-bound
    _ub_int1_new = (_ub_int1_new.array() - _lb_int1_new.array() < 0).select(_ub_int1_new.array() + 2.0*M_PI, _ub_int1_new);
    _ub_int2_new = (_ub_int2_new.array() - _lb_int2_new.array() < 0).select(_ub_int2_new.array() + 2.0*M_PI, _ub_int2_new);

    //if an interval length exceeds 2.0*M_PI, set to 2.0*M_PI
    _lb_int1_new = (_ub_int1_new.array() - _lb_int1_new.array() > 2.0*M_PI).select(0, _lb_int1_new);
    _ub_int1_new = (_ub_int1_new.array() - _lb_int1_new.array() > 2.0*M_PI).select(2.0*M_PI, _ub_int1_new);

    //create new interval matrices
    (*Int1_new) = MatrixXd::Zero(_n_rows, 2);
    (*Int2_new) = MatrixXd::Zero(_n_rows, 2);

    (*Int1_new).col(0) = _lb_int1_new;
    (*Int1_new).col(1) = _ub_int1_new;
    (*Int2_new).col(0) = _lb_int2_new;
    (*Int2_new).col(1) = _ub_int2_new;

    #if DEBUG
        MatrixXd test_in = MatrixXd::Zero(_n_rows, 4);
        MatrixXd m_test_in = MatrixXd::Zero(_n_rows, 4);
        MatrixXd test_out = MatrixXd::Zero(_n_rows, 4);
        int row = 7;
        test_in   << _lb_int1, _ub_int1, _lb_int2, _ub_int2;
        m_test_in << _m_lb_int1, _m_ub_int1, _m_lb_int2, _m_ub_int2;
        test_out  << _lb_int1_new, _ub_int1_new, _lb_int2_new, _ub_int2_new;
        cout << "test_in: "                 << test_in.row(row)                 << "\n";
        cout << "m_test_in: "               << m_test_in.row(row)               << "\n";
        cout << "logic_I1disjointI2: "      << _logic_I1disjointI2.row(row)      << "\n";
        cout << "logic_I1contI2: "          << _logic_I1contI2.row(row)          << "\n";
        cout << "logic_I2contI1: "          << _logic_I2contI1.row(row)          << "\n";
        cout << "logic_I1overlapI2: "       << _logic_I1overlapI2.row(row)       << "\n";
        cout << "logic_I2overlapI1: "       << _logic_I2overlapI1.row(row)       << "\n";    
        cout << "logic_I1doubleoverlapI2: " << _logic_I1doubleoverlapI2.row(row) << "\n";
        cout << "logic_I2doubleoverlapI1: " << _logic_I2doubleoverlapI1.row(row) << "\n";
        cout << "test_out: "                << test_out.row(row)                << "\n";
    #endif
   
}
