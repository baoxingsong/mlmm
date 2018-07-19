//
// Created by Baoxing song on 14.06.18.
//

#ifndef MLMM_CPP_MT_GEMMA_SERVICE_H
#define MLMM_CPP_MT_GEMMA_SERVICE_H

#include "../impl/impl.h"
#include "../model/model.h"
#include <map>

void mt_gemma_test ( const std::string & phenotype_path, const std::string & genotype_path,
                     const std::string & kinship_file, const double & maf);
void mt_gemma_test_MultiAllic ( const std::string & phenotype_path, const std::string & genotype_path,
                                const std::string & kinship_file, const double & maf);
#endif //MLMM_CPP_MT_GEMMA_SERVICE_H
