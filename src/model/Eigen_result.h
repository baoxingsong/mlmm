//
// Created by baoxing on 2/23/18.
//

#ifndef MLMM_CPP_EIGEN_RESULT_H
#define MLMM_CPP_EIGEN_RESULT_H

#include <vector>
#include "My_matrix.h"
#include "My_Vector.h"
#include "../util/util.h"
class Eigen_result {
    private:
        My_Vector<double> eigen_values;
        My_matrix<double> eigen_vectors; //each column is a eigen vector
    public:
        Eigen_result( const My_Vector<double> & _eigen_values, const My_matrix<double> & _eigen_vectors);
        const My_Vector<double> &get_eigen_values() const;
        const My_matrix<double> &get_eigen_vectors() const;

};

#endif //MLMM_CPP_EIGEN_RESULT_H
