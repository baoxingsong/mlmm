//
// Created by baoxing on 2/25/18.
//

#ifndef MLMM_CPP_EMMAX_SERVICE_H
#define MLMM_CPP_EMMAX_SERVICE_H

#include "../impl/impl.h"
#include "../model/model.h"
#include <map>
void emma_test ( phenotype_impl & pi, Kinship_matrix_impl & k_i, Genotype & genotype,  const double & man_l, const double & man_u );
void emmax_test( phenotype_impl & pi, Kinship_matrix_impl & k_i, Genotype & genotype,  const double & man_l, const double & man_u );
void emmax_test_multi_allic_single_test( phenotype_impl & pi, Kinship_matrix_impl & k_i, Genotype & genotype, const double & man_l, const double & man_u );
void emmax_test_multi_allic_multi_test_null_model( phenotype_impl & pi, Kinship_matrix_impl & k_i, Genotype & genotype, const double & man_l, const double & man_u );
void emmax_test_multi_allic_multi_test_full_model( phenotype_impl & pi, Kinship_matrix_impl & k_i, Genotype & genotype, const double & man_l, const double & man_u );

void emma_test ( const std::string & phenotype_path, const std::string & genotype_path, const std::string & kinship_file, const double & maf);
void emmax_test( const std::string & phenotype_path, const std::string & genotype_path, const std::string & kinship_file, const double & maf);
void emmax_test_multi_allic_single_test( const std::string & phenotype_path, const std::string & genotype_path, const std::string & kinship_file, const double & maf);
void emmax_test_multi_allic_multi_test_null_model( const std::string & phenotype_path, const std::string & genotype_path, const std::string & kinship_file, const double & maf);
void emmax_test_multi_allic_multi_test_full_model( const std::string & phenotype_path, const std::string & genotype_path, const std::string & kinship_file, const double & maf);

#endif //MLMM_CPP_EMMAX_SERVICE_H
