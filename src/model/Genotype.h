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
#include <iostream>
class Genotype {
    private:
        My_matrix<unsigned int > genotype_matrix;
        unsigned long number_of_individual;
        unsigned long number_of_variant;
        std::vector<Variant> variant_Vector;
        std::vector<Individual> individual_Vector;
    public:
        Genotype();
        Genotype(const unsigned long& _number_of_individul,
            const unsigned long& _number_of_variant, const std::vector<Variant>& _variant_Vector,
            const std::vector<Individual>& _individual_Vector, const std::vector<int*>& genotypes, const std::string& format);
        My_matrix<unsigned int> &  get_genotype_matrix();
        const unsigned long& get_number_of_individual() const;
        const unsigned long& get_number_of_variant() const;
        const std::vector<Variant>& get_variant_Vector() const;
        const std::vector<Individual>& get_individual_Vector() const;
        void onlyKeepThoseIndividuls(const std::vector <std::string> & individual_ids);
};

#endif //MLMM_CPP_GENOTYPE_H
