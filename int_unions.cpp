#include <eigen3/Eigen/Dense>
#include <utility>
#include <iostream>
#include "include/poc_function.hpp"

using namespace Eigen;
using namespace casadi;
using namespace std;


void POCFunction::_int_unions(const MatrixXd* Int1, VectorXd* size_int1, const MatrixXd* Int2, VectorXd* size_int2, 
                              MatrixXd *Int1_new, MatrixXd *Int2_new) {
    
    //extract lower and upper bounds of int1 and int2
    VectorXd lb_int1 = (*Int1).col(0);
    VectorXd ub_int1 = (*Int1).col(1);
    VectorXd lb_int2 = (*Int2).col(0);
    VectorXd ub_int2 = (*Int2).col(1);
    VectorXd m_lb_int1 = _modified_mod(lb_int1);
    VectorXd m_ub_int1 = _modified_mod(ub_int1);
    VectorXd m_lb_int2 = _modified_mod(lb_int2);
    VectorXd m_ub_int2 = _modified_mod(ub_int2);

    //initialize new variables for int1 and int2
    VectorXd lb_int1_new = lb_int1;
    VectorXd ub_int1_new = ub_int1;
    VectorXd lb_int2_new = lb_int2;
    VectorXd ub_int2_new = ub_int2;
    VectorXd temp_int_size_1 = VectorXd::Zero(_n_rows);
    VectorXd temp_int_size_2 = VectorXd::Zero(_n_rows);
    VectorXd temp_int_size_3 = VectorXd::Zero(_n_rows);

    //check if interval is of larger size than 0
    VectorXd logic_I1greater0 = VectorXd::Zero(_n_rows);
    VectorXd logic_I2greater0 = VectorXd::Zero(_n_rows);
    logic_I1greater0 = (abs((*size_int1).array()) > 0).select(1, logic_I1greater0);
    logic_I2greater0 = (abs((*size_int2).array()) > 0).select(1, logic_I2greater0);

    //check wheather both intervals are equal in mod(x, 2*pi)
    VectorXd logic_I1eqI2 = VectorXd::Zero(_n_rows);
    logic_I1eqI2 = (m_lb_int1.array() == m_lb_int2.array() && m_ub_int1.array() == m_ub_int2.array() && logic_I1greater0.array() == 1).select(1, logic_I1eqI2);

    //check wheather one interval contains the other
    VectorXd logic_I1contI2 = VectorXd::Zero(_n_rows);
    VectorXd logic_I2contI1 = VectorXd::Zero(_n_rows);
    logic_I1contI2 = (m_lb_int1.array() <= m_lb_int2.array() && m_ub_int1.array() >= m_ub_int2.array() && m_ub_int2.array() >= m_lb_int1.array() && logic_I1eqI2.array() == 0 && logic_I2greater0.array() == 1).select(1, logic_I1contI2);
    logic_I2contI1 = (m_lb_int2.array() <= m_lb_int1.array() && m_ub_int2.array() >= m_ub_int1.array() && m_ub_int1.array() >= m_lb_int2.array() && logic_I1eqI2.array() == 0 && logic_I1greater0.array() == 1).select(1, logic_I2contI1);
    
    //check wheather there is an overlap among both intervals
    VectorXd logic_I1overlapI2 = VectorXd::Zero(_n_rows);
    VectorXd logic_I2overlapI1 = VectorXd::Zero(_n_rows);
    logic_I1overlapI2 = (m_ub_int1.array() >= m_lb_int2.array() && m_lb_int1.array() <= m_lb_int2.array() && logic_I1contI2.array() == 0 && logic_I2greater0.array() == 1).select(1, logic_I1overlapI2);
    logic_I2overlapI1 = (m_ub_int2.array() >= m_lb_int1.array() && m_lb_int2.array() <= m_lb_int1.array() && logic_I2contI1.array() == 0 && logic_I1greater0.array() == 1).select(1, logic_I2overlapI1);
     
    //check if intervals are disjoint
    VectorXd logic_I1disjointI2= VectorXd::Zero(_n_rows);
    logic_I1disjointI2 = (logic_I1eqI2.array() == 0 && logic_I1contI2.array() == 0 && logic_I2contI1.array() == 0 && logic_I1overlapI2.array() == 0 && logic_I2overlapI1.array() == 0).select(1, logic_I1disjointI2);

    //if int1 contains int2, keep int1 and set int2 to zero
    lb_int1_new = (logic_I1contI2.array() == 1).select(lb_int1, lb_int1_new);
    ub_int1_new = (logic_I1contI2.array() == 1).select(ub_int1, ub_int1_new);
    lb_int2_new = (logic_I1contI2.array() == 1).select(0, lb_int2_new);
    ub_int2_new = (logic_I1contI2.array() == 1).select(0, ub_int2_new);

    //if int2 contains int1, write int2 into int1 and set int2 to zero
    lb_int1_new = (logic_I2contI1.array() == 1).select(lb_int2, lb_int1_new);
    ub_int1_new = (logic_I2contI1.array() == 1).select(ub_int2, ub_int1_new);
    lb_int2_new = (logic_I2contI1.array() == 1).select(0, lb_int2_new);
    ub_int2_new = (logic_I2contI1.array() == 1).select(0, ub_int2_new); 

    //if both contain each other, i.e., both intervals are identical, set I2 to zero
    lb_int2_new = (logic_I1eqI2.array() == 1).select(0, lb_int2_new);
    ub_int2_new = (logic_I1eqI2.array() == 1).select(0, ub_int2_new);

    //if disjoint, leave them as they are
    lb_int1_new = (logic_I1disjointI2.array() == 1).select(lb_int1, lb_int1_new);
    ub_int1_new = (logic_I1disjointI2.array() == 1).select(ub_int1, ub_int1_new);
    lb_int2_new = (logic_I1disjointI2.array() == 1).select(lb_int2, lb_int2_new);
    ub_int2_new = (logic_I1disjointI2.array() == 1).select(ub_int2, ub_int2_new);

    //if int1 and int2 overlap, combine into int1 and set int2 to zero
    lb_int1_new = (logic_I1overlapI2.array() == 1).select(lb_int1, lb_int1_new);
    ub_int1_new = (logic_I1overlapI2.array() == 1).select(ub_int2, ub_int1_new);
    lb_int2_new = (logic_I1overlapI2.array() == 1).select(0, lb_int2_new);
    ub_int2_new = (logic_I1overlapI2.array() == 1).select(0, ub_int2_new);

    //if int2 and int1 overlap, combine into int1 and set int2 to zero
    lb_int1_new = (logic_I2overlapI1.array() == 1).select(lb_int2, lb_int1_new);
    ub_int1_new = (logic_I2overlapI1.array() == 1).select(ub_int1, ub_int1_new);
    lb_int2_new = (logic_I2overlapI1.array() == 1).select(0, lb_int2_new);
    ub_int2_new = (logic_I2overlapI1.array() == 1).select(0, ub_int2_new);

    //create new interval matrices
    (*Int1_new) = MatrixXd::Zero(_n_rows, 2);
    (*Int2_new) = MatrixXd::Zero(_n_rows, 2);

    (*Int1_new).col(0) = lb_int1_new;
    (*Int1_new).col(1) = ub_int1_new;
    (*Int2_new).col(0) = lb_int2_new;
    (*Int2_new).col(1) = ub_int2_new;

    SX int_size_1 = _eigenToCasadiSX(*size_int1);
    SX int_size_2 = _eigenToCasadiSX(*size_int2);
   
}
