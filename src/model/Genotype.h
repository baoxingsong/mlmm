//
// Created by baoxing on 2/17/18.
//

#ifndef MLMM_CPP_GENOTYPE_H
#define MLMM_CPP_GENOTYPE_H

#include <string>
#include <vector>
#include "Variant.h"
#include "Individual.h"
#include "My_matrix.h"
#include "My_Vector.h"
#include <iostream>
class Genotype {
    private:
        My_matrix<int8_t > genotype_matrix;
        size_t number_of_individual;
        size_t number_of_variant;
        Variant * variant_Vector;
        // I could not use the My_Vector here.
        // Under Mac it is ok to use My_Vector.
        // While it seems the generics technology under linux does not support custom class
        Individual * individual_Vector;
    public:
        Genotype();
        Genotype(const size_t & number_of_individuals, const size_t & number_of_variants);
        My_matrix<int8_t> &  get_genotype_matrix();
        Variant * get_variant_Vector();
        Individual * get_individual_Vector();
        void onlyKeepThoseIndividuls(const std::vector <std::string> & individual_ids);
        size_t get_number_of_individual();
        size_t get_number_of_variant();
};
#endif //MLMM_CPP_GENOTYPE_H
