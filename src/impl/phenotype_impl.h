//
// Created by baoxing on 2/25/18.
//

#ifndef MLMM_CPP_PHENOTYPE_IMPL_H
#define MLMM_CPP_PHENOTYPE_IMPL_H

#include "../model/model.h"
#include <vector>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cstdio>
class phenotype_impl {
    private:
        std::vector <std::string> individual_ids;
        std::vector <std::string> family_ids;
        My_Vector <double> phenotypes;
    public:
        phenotype_impl(const std::string & file_path, const std::string & format);
        const std::vector<std::string> &getIndividual_ids() const;
        void setIndividual_ids(const std::vector<std::string> &individual_ids);
        const std::vector<std::string> &getFamily_ids() const;
        void setFamily_ids(const std::vector<std::string> &family_ids);
        const My_Vector<double> &getPhenotypes() const;
        void setPhenotypes(const My_Vector<double> &phenotypes);
};


#endif //MLMM_CPP_PHENOTYPE_IMPL_H
