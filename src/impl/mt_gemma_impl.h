//
// Created by Baoxing song on 06.06.18.
//

#ifndef MLMM_CPP_MT_GEMMA_IMPL_H
#define MLMM_CPP_MT_GEMMA_IMPL_H

#include "../util/util.h"
#include "../model/model.h"
#include "My_linear_algebra_impl.h"
#include "./gemma_impl.h"
#include "./emmax_impl.h"

#include "phenotype_impl.h"
#include "Kinship_matrix_impl.h"


void AnalyzePlink(const Eigen_result & eigen_r, const My_matrix<double> & UtW, const My_matrix<double> & UtY,
                  const char & method, Kinship_matrix_impl & k_i,
                  Genotype & genotype, const My_matrix<double> & phenotype);



#endif //MLMM_CPP_MT_GEMMA_IMPL_H
