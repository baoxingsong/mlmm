//
// Created by baoxing on 2/17/18.
//

#ifndef MLMM_CPP_KINSHIP_MATRIX_IMPL_H
#define MLMM_CPP_KINSHIP_MATRIX_IMPL_H
#include "../util/util.h"
#include "../model/model.h"
#include "My_linear_algebra_impl.h"
#include <fstream>
#include <iostream>


class Kinship_matrix_impl {
    private:
        My_matrix<double> kin; //each column is a eigen vector
    public:
        Kinship_matrix_impl( Genotype& genotype );
        void scale_k();
        Kinship_matrix_impl(const std::string & kin_file_path);
        const My_matrix<double> &getKin() const;
        void setKin(const My_matrix<double> &kin);
};

#endif //MLMM_CPP_KINSHIP_MATRIX_IMPL_H
