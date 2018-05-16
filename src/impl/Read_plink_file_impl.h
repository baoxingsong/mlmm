//
// Created by baoxing on 2/17/18.
//

#ifndef MLMM_CPP_READ_PLINK_FILE_IMPL_H
#define MLMM_CPP_READ_PLINK_FILE_IMPL_H
#include <string>
#include <vector>
#include <fstream>
#include <map>
#include <regex>
#include <sstream>
#include <assert.h>

#include "../model/model.h"
#include "../util/util.h"
void Read_ped_file(const std::string & mapFile, const std::string & pedFile, Genotype & genotype);
void Read_tped_file(const std::string & tfamFile, const std::string & tpedFile, Genotype & genotype);
#endif //MLMM_CPP_READ_PLINK_FILE_IMPL_H
