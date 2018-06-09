//
// Created by song on 5/29/18.
//

#ifndef MLMM_CPP_GEMMA_SERVICE_H
#define MLMM_CPP_GEMMA_SERVICE_H


#include "../impl/impl.h"
#include "../model/model.h"
#include <map>

void gemma_test ( phenotype_impl & pi, Kinship_matrix_impl & k_i, Genotype & genotype, const double & man_l, const double & man_u );
void gemma_test ( const std::string & phenotype_path, const std::string & genotype_path, const std::string & kinship_file, const double & maf);
#endif //MLMM_CPP_GEMMA_SERVICE_H
