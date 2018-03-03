//
// Created by baoxing on 2/23/18.
//

#include "Eigen_result.h"

Eigen_result::Eigen_result( const My_Vector<double> & _eigen_values, const My_matrix<double> & _eigen_vectors){
    this->eigen_values=_eigen_values;
    this->eigen_vectors=_eigen_vectors;
}

const My_Vector<double> &Eigen_result::get_eigen_values() const {
    return eigen_values;
}

const My_matrix<double> &Eigen_result::get_eigen_vectors() const {
    return eigen_vectors;
}

//void Eigen_result::first_several(const unsigned long & size){
//    this->eigen_values.set_length(size);
//    this->eigen_vectors.set_num_row(size);
//}
