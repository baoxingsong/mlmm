//
// Created by Baoxing song on 06.06.18.
//

#ifndef MLMM_CPP_PHENOTYPE_MUL_IMPL_H
#define MLMM_CPP_PHENOTYPE_MUL_IMPL_H

#include "../model/model.h"
#include <vector>
#include <fstream>
#include <cstdlib>
#include <cstdio>

class phenotype_mul_impl {
    private:
        std::vector <std::string> individual_ids;
        std::vector <std::string> family_ids;
        My_matrix <double> phenotypes;
    public:
        phenotype_mul_impl();
        phenotype_mul_impl(const std::string & file_path, const std::string & format);
        const std::vector<std::string> &getIndividual_ids() const;
        void setIndividual_ids(const std::vector<std::string> &individual_ids);
        const std::vector<std::string> &getFamily_ids() const;
        void setFamily_ids(const std::vector<std::string> &family_ids);
        const My_matrix<double> &getPhenotypes() const;
        void setPhenotypes(const My_matrix<double> &phenotypes);
};

#endif //MLMM_CPP_PHENOTYPE_MUL_IMPL_H
